from qgis.processing import alg
from qgis.core import (
    QgsProject,
    QgsCoordinateTransform,
    QgsCoordinateReferenceSystem,
)
import math

# ---------------------------
# CONFIG
# ---------------------------
STEP = 0.01  # m: kleiner = glattere Kanten, mehr Punkte
COLOR = (8, 243, 8)  # gruen
MATCH_PERP_TOL = 2.5 * STEP  # Matching-Toleranz fuer Intervall-Tracking


# ---------------------------
# AOG LocalPlane helpers
# ---------------------------
def meters_per_degree_lat(lat):
    r = math.radians(lat)
    return (
        111132.92
        - 559.82 * math.cos(2 * r)
        + 1.175 * math.cos(4 * r)
        - 0.0023 * math.cos(6 * r)
    )

def meters_per_degree_lon(lat):
    r = math.radians(lat)
    return (
        111412.84 * math.cos(r)
        - 93.5 * math.cos(3 * r)
        + 0.118 * math.cos(5 * r)
    )

def to_local(lat, lon, o_lat, o_lon, mpl):
    n = (lat - o_lat) * mpl
    e = (lon - o_lon) * meters_per_degree_lon(lat)
    return e, n

def read_startfix(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        s = l.strip()
        if s.startswith("StartFix"):
            if "," not in s:
                s = lines[i + 1].strip()
            else:
                s = s.replace("StartFix", "").strip()
            a, b = s.split(",")
            return float(a), float(b)
    raise RuntimeError("StartFix not found")

def cross(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


# ---------------------------
# PCA fallback direction
# ---------------------------
def pca_direction(points):
    cx = sum(p[0] for p in points) / len(points)
    cy = sum(p[1] for p in points) / len(points)
    xx = yy = xy = 0.0
    for x, y in points:
        dx = x - cx
        dy = y - cy
        xx += dx * dx
        yy += dy * dy
        xy += dx * dy
    ang = 0.5 * math.atan2(2 * xy, xx - yy)
    return math.cos(ang), math.sin(ang)


# ---------------------------
# Line direction helper (optional)
# ---------------------------
def direction_from_line(line_layer, xform_line, o_lat, o_lon, mpl):
    feat = next(line_layer.getFeatures(), None)
    if feat is None:
        return None

    geom = feat.geometry()

    # MultiLineString -> laengsten Teil nehmen
    if geom.isMultipart():
        parts = geom.asMultiPolyline()
        if not parts:
            return None
        line = max(parts, key=len)
    else:
        line = geom.asPolyline()

    if len(line) < 2:
        return None

    p0 = line[0]
    p1 = line[-1]

    ll0 = xform_line.transform(p0.x(), p0.y())
    ll1 = xform_line.transform(p1.x(), p1.y())

    e0, n0 = to_local(ll0.y(), ll0.x(), o_lat, o_lon, mpl)
    e1, n1 = to_local(ll1.y(), ll1.x(), o_lat, o_lon, mpl)

    dx = e1 - e0
    dy = n1 - n0
    L = math.hypot(dx, dy)
    if L == 0:
        return None
    return dx / L, dy / L


# ---------------------------
# Slicing helpers (holes supported)
# ---------------------------
def project_ring(ring_pts_local, dir_x, dir_y, perp_x, perp_y):
    # ring_pts_local: [(e,n), ...] ohne closing point
    return [(e, n,
             e * dir_x + n * dir_y,
             e * perp_x + n * perp_y)
            for (e, n) in ring_pts_local]

def hits_for_ring(proj_ring, a, perp_x, perp_y):
    # returns list of tuples: (perp_scalar, e, n)
    hits = []
    n = len(proj_ring)
    for i in range(n):
        p1 = proj_ring[i]
        p2 = proj_ring[(i + 1) % n]
        a1 = p1[2]
        a2 = p2[2]
        if a1 == a2:
            continue
        # crossing or touching
        if (a1 - a) * (a2 - a) <= 0:
            t = (a - a1) / (a2 - a1)
            if 0.0 <= t <= 1.0:
                e = p1[0] + t * (p2[0] - p1[0])
                n_ = p1[1] + t * (p2[1] - p1[1])
                perp = e * perp_x + n_ * perp_y
                hits.append((perp, e, n_))
    hits.sort(key=lambda x: x[0])
    return hits

def intervals_from_hits(hits):
    # hits sorted by perp; pairs define inside-intervals (even-odd fill)
    # return list of intervals [(Lhit, Rhit)] each hit=(perp,e,n)
    if len(hits) < 2:
        return []
    # ensure even
    if len(hits) % 2 == 1:
        # drop the most suspicious one (closest pair issue) -> conservative
        hits = hits[:-1]
    intervals = []
    for i in range(0, len(hits), 2):
        L = hits[i]
        R = hits[i + 1]
        if R[0] > L[0]:
            intervals.append((L, R))
    return intervals

def subtract_intervals(base_intervals, cut_intervals):
    # both lists are (Lhit,Rhit) with perp scalars.
    # subtract each cut interval from base intervals.
    out = base_intervals
    for cL, cR in cut_intervals:
        new_out = []
        for bL, bR in out:
            # no overlap
            if cR[0] <= bL[0] or cL[0] >= bR[0]:
                new_out.append((bL, bR))
                continue
            # left remainder
            if bL[0] < cL[0]:
                new_out.append((bL, cL))
            # right remainder
            if cR[0] < bR[0]:
                new_out.append((cR, bR))
        out = new_out
        if not out:
            break
    return out


# ---------------------------
# Patch writer
# ---------------------------
def write_patch(f, left_pts, right_pts, color):
    # left_pts/right_pts lists of (e,n)
    if not left_pts or not right_pts:
        return
    if len(left_pts) != len(right_pts):
        return

    # need at least 2 stations; if only 1, duplicate
    if len(left_pts) == 1:
        left_pts = [left_pts[0], left_pts[0]]
        right_pts = [right_pts[0], right_pts[0]]

    # CCW enforce
    if cross(left_pts[0], right_pts[0], left_pts[1]) < 0:
        left_pts, right_pts = right_pts, left_pts

    count = 1 + 2 * len(left_pts)
    f.write(f"{count}\n")
    f.write(f"{color[0]},{color[1]},{color[2]}\n")
    for l, r in zip(left_pts, right_pts):
        f.write(f"{l[0]:.3f},{l[1]:.3f},0\n")
        f.write(f"{r[0]:.3f},{r[1]:.3f},0\n")


# ---------------------------
# Processing Algorithm
# ---------------------------
@alg(
    name="aog_sections_from_polygons_with_holes",
    label="AOG Sections from Polygons (holes + optional direction line=",
    group="AOG",
    group_label="AOG"
)
@alg.input(type=alg.VECTOR_LAYER, name="POLY", label="Polygon layer")
@alg.input(type=alg.VECTOR_LAYER, name="LINE", label="Direction line (optional)", optional=True)
@alg.input(type=alg.FILE, name="FIELD", label="Field.txt")
@alg.input(type=alg.FILE_DEST, name="OUT", label="Sections.txt")
def run(instance, parameters, context, feedback, inputs):
    """Exports polygon silhouettes to AOG Sections.txt, supports inner holes.
    Uses direction line if provided, otherwise PCA per polygon.
    """

    poly = instance.parameterAsVectorLayer(parameters, "POLY", context)
    line = instance.parameterAsVectorLayer(parameters, "LINE", context)
    field = instance.parameterAsFile(parameters, "FIELD", context)
    out = instance.parameterAsFileOutput(parameters, "OUT", context)

    if poly is None or poly.geometryType() != 2:
        raise RuntimeError("Polygon layer required")

    # origin
    o_lat, o_lon = read_startfix(field)
    mpl = meters_per_degree_lat(o_lat)

    crs_wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
    xform_poly = QgsCoordinateTransform(poly.crs(), crs_wgs84, QgsProject.instance())

    # optional global direction from line
    dir_from_line = None
    if line is not None:
        xform_line = QgsCoordinateTransform(line.crs(), crs_wgs84, QgsProject.instance())
        dir_from_line = direction_from_line(line, xform_line, o_lat, o_lon, mpl)

    with open(out, "w", encoding="ascii") as f:
        for feat in poly.getFeatures():
            geom = feat.geometry()
            if geom.isEmpty():
                continue

            polygons = geom.asMultiPolygon() if geom.isMultipart() else [geom.asPolygon()]

            for poly_rings in polygons:
                if not poly_rings or len(poly_rings[0]) < 4:
                    continue

                outer_ring = poly_rings[0]
                hole_rings = poly_rings[1:] if len(poly_rings) > 1 else []

                # --- convert rings to LocalPlane (each ring without closing point) ---
                def ring_to_local(ring):
                    pts_local = []
                    for p in ring[:-1]:
                        ll = xform_poly.transform(p.x(), p.y())
                        pts_local.append(to_local(ll.y(), ll.x(), o_lat, o_lon, mpl))
                    return pts_local

                outer_local = ring_to_local(outer_ring)
                if len(outer_local) < 3:
                    continue
                holes_local = [ring_to_local(hr) for hr in hole_rings if len(hr) >= 4]
                holes_local = [h for h in holes_local if len(h) >= 3]

                # --- choose direction ---
                if dir_from_line is not None:
                    dir_x, dir_y = dir_from_line
                else:
                    dir_x, dir_y = pca_direction(outer_local)

                perp_x, perp_y = -dir_y, dir_x

                # project rings
                outer_proj = project_ring(outer_local, dir_x, dir_y, perp_x, perp_y)
                holes_proj = [project_ring(h, dir_x, dir_y, perp_x, perp_y) for h in holes_local]

                along_vals = [p[2] for p in outer_proj]
                min_a, max_a = min(along_vals), max(along_vals)
                if max_a - min_a < 1e-6:
                    continue

                # --- build per-slice intervals (holes subtracted) ---
                slice_data = []  # list of (a, [ (Lpt(e,n), Rpt(e,n), mid_perp) , ... ])
                a = min_a
                while a <= max_a + 1e-6:
                    outer_hits = hits_for_ring(outer_proj, a, perp_x, perp_y)
                    outer_intervals = intervals_from_hits(outer_hits)
                    if not outer_intervals:
                        a += STEP
                        continue

                    # subtract all hole intervals
                    intervals = outer_intervals
                    for hp in holes_proj:
                        h_hits = hits_for_ring(hp, a, perp_x, perp_y)
                        h_intervals = intervals_from_hits(h_hits)
                        if h_intervals:
                            intervals = subtract_intervals(intervals, h_intervals)
                            if not intervals:
                                break

                    if intervals:
                        segs = []
                        for Lh, Rh in intervals:
                            # Lh/Rh are hits (perp,e,n)
                            Lpt = (Lh[1], Lh[2])
                            Rpt = (Rh[1], Rh[2])
                            mid = 0.5 * (Lh[0] + Rh[0])
                            segs.append((Lpt, Rpt, mid))
                        # sort segs by mid perp (stable ordering)
                        segs.sort(key=lambda s: s[2])
                        slice_data.append((a, segs))

                    a += STEP

                if len(slice_data) < 2:
                    continue

                # --- Track segments across slices into continuous strips (patches) ---
                # Each active strip: id -> {left:[...], right:[...], last_mid:float, last_a:float}
                next_id = 1
                active = {}

                def close_strip(strip):
                    write_patch(f, strip["left"], strip["right"], COLOR)

                for (a, segs) in slice_data:
                    # match segs to active strips by nearest mid
                    used_active = set()
                    used_segs = set()

                    # build match candidates
                    seg_mids = [s[2] for s in segs]
                    # greedy matching by nearest mid
                    for sid, strip in list(active.items()):
                        best_j = None
                        best_d = None
                        for j, mid in enumerate(seg_mids):
                            if j in used_segs:
                                continue
                            d = abs(mid - strip["last_mid"])
                            if best_d is None or d < best_d:
                                best_d = d
                                best_j = j
                        if best_j is not None and best_d is not None and best_d <= MATCH_PERP_TOL:
                            # match
                            Lpt, Rpt, mid = segs[best_j]
                            strip["left"].append(Lpt)
                            strip["right"].append(Rpt)
                            strip["last_mid"] = mid
                            strip["last_a"] = a
                            used_active.add(sid)
                            used_segs.add(best_j)

                    # close strips not used this slice
                    for sid in list(active.keys()):
                        if sid not in used_active:
                            close_strip(active[sid])
                            del active[sid]

                    # start new strips for unmatched segs
                    for j, (Lpt, Rpt, mid) in enumerate(segs):
                        if j in used_segs:
                            continue
                        active[next_id] = {
                            "left": [Lpt],
                            "right": [Rpt],
                            "last_mid": mid,
                            "last_a": a
                        }
                        next_id += 1

                # close remaining
                for sid in list(active.keys()):
                    close_strip(active[sid])
                    del active[sid]

    return {}
