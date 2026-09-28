"""Whether a bond passes through a ring.

One question asked in three places: the ring cure asks it forward (would this new ring
close around a bond?) and backward (did dragging these monomers push a bond through a
ring already there?), and the analysis asks it of a finished network.  The geometry is
the same each time and lives here so those three cannot drift apart.

A ring is treated as a fan of triangles from its centroid rather than as a plane,
because the thing being tested is not always flat: at candidate time in a ring cure it
is three reactive groups a search radius apart, and even a real ring puckers.

**Image before you average.**  Every ring atom is brought into the same periodic image as
an anchor atom before the centroid is taken.  A centroid of raw wrapped coordinates puts
the centre of any boundary-straddling ring in the middle of the box and its triangle fan
across the whole cell, which reports a piercing for nearly everything nearby.  About 15%
of rings straddle a boundary in a typical crosslinked system, so this is the normal case
rather than a corner one: one detector written without it reported 212 piercings in an
uncured melt, where the true answer is zero by construction.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import numpy as np


def mic(r, box):
    """Minimum-image displacement."""
    box = np.asarray(box, dtype=float)
    return r - np.round(r / box) * box


def segment_pierces_triangle(p0, p1, a, b, c, eps=1e-12):
    """Moller-Trumbore, restricted to a segment rather than an infinite ray.

    Returns:
        float or None: the fraction along p0->p1 at which it crosses, or None
    """
    e1, e2, d = b - a, c - a, p1 - p0
    h = np.cross(d, e2)
    det = float(np.dot(e1, h))
    if abs(det) < eps:
        return None
    inv = 1.0 / det
    s = p0 - a
    u = inv * float(np.dot(s, h))
    if u < 0.0 or u > 1.0:
        return None
    q = np.cross(s, e1)
    v = inv * float(np.dot(d, q))
    if v < 0.0 or u + v > 1.0:
        return None
    t = inv * float(np.dot(e2, q))
    return t if 0.0 < t < 1.0 else None


def ring_frame(ring, positions, box):
    """The ring's own coordinate frame, imaged so it does not straddle a boundary.

    Args:
        ring (list): atom indices in cyclic order
        positions (dict): atom index -> position, in nm
        box (array-like): box diagonal, in nm

    Returns:
        tuple: (points, centroid, unit normal, radius) or None if any atom is missing
    """
    if any(i not in positions for i in ring):
        return None
    anchor = np.asarray(positions[ring[0]], dtype=float)
    P = np.array([anchor + mic(np.asarray(positions[i], dtype=float) - anchor, box)
                  for i in ring])
    O = P.mean(axis=0)
    # normal from the averaged cross products of successive centre-to-vertex vectors;
    # sign is arbitrary, which is fine because every use here takes an absolute angle
    V = P - O
    n = np.zeros(3)
    for k in range(len(V)):
        n = n + np.cross(V[k], V[(k + 1) % len(V)])
    norm = np.linalg.norm(n)
    n = n / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])
    return P, O, n, float(np.linalg.norm(V, axis=1).max())


def pierces_ring(ring, bond, positions, box):
    """Whether `bond` passes through `ring`, and where.

    Args:
        ring (list): the ring's atom indices in cyclic order
        bond (tuple): (ai, aj); a bond sharing an atom with the ring never counts
        positions (dict): atom index -> position, in nm
        box (array-like): box diagonal, in nm

    Returns:
        dict or None: ``point``, ``offset`` from the ring centre in nm, ``offset_frac``
            of the ring radius, ``angle`` of the bond to the ring normal in degrees, and
            ``length`` of the bond in nm
    """
    ai, aj = int(bond[0]), int(bond[1])
    if ai in ring or aj in ring:
        return None
    frame = ring_frame([int(x) for x in ring], positions, box)
    if frame is None or ai not in positions or aj not in positions:
        return None
    P, O, n, radius = frame
    anchor = np.asarray(positions[int(ring[0])], dtype=float)
    p0 = anchor + mic(np.asarray(positions[ai], dtype=float) - anchor, box)
    # the far end is imaged against the near end, so the bond stays one segment
    d = mic(np.asarray(positions[aj], dtype=float) - np.asarray(positions[ai], dtype=float), box)
    p1 = p0 + d
    for k in range(len(P)):
        t = segment_pierces_triangle(p0, p1, O, P[k], P[(k + 1) % len(P)])
        if t is not None:
            point = p0 + t * (p1 - p0)
            offset = float(np.linalg.norm(point - O))
            length = float(np.linalg.norm(d))
            cos = abs(float(np.dot(d / length, n))) if length > 0 else 0.0
            return {'point': point, 'offset': offset,
                    'offset_frac': offset / radius if radius > 0 else float('nan'),
                    'angle_from_normal': float(np.degrees(np.arccos(min(1.0, cos)))),
                    'length': length}
    return None
