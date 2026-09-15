"""Triazine-to-cyanate-cap postcure repair.

Dismantles every incomplete triazine crosslinker (one with fewer than the
configured ``full_bond_count`` BPA-O bonds) into three independent -C#N
fragments, each then attached either in place (to the BPA-O it was already
bonded to during cure) or to a free, unreacted BPA-OH within
``cap_search_radius``.  Atom conservation is exact: across the whole
system, the count of fragments freed from incomplete rings equals the
count of unreacted BPA-OHs, so every free fragment finds a home.

The surgery is staged so that atom-index references stay valid until the
final batched deletion of sacrificial H atoms; intermediate phases edit
bonds, atom attributes, and residue tags only.
"""
import logging
from collections import defaultdict

import numpy as np
import pandas as pd

from . import topology_surgery as ts

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Planning helpers
# ---------------------------------------------------------------------------


def _find_crosslinker_residues(TC, residue_name):
    """Return {resNum: [globalIdx, ...]} for every residue with this name."""
    A = TC.Coordinates.A
    sub = A[A['resName'] == residue_name][['globalIdx', 'resNum', 'atomName']]
    out = {}
    for rn, g in sub.groupby('resNum'):
        out[int(rn)] = {row.atomName: int(row.globalIdx) for row in g.itertuples()}
    return out


def _bonded_bridge_o(TC, c_idx, bridge_residue, bridge_oxygens):
    """If atom c_idx is bonded to a (bridge_residue).(one of bridge_oxygens),
    return the bridge-O globalIdx; else return None."""
    A = TC.Coordinates.A
    for nbr in TC.Topology.bondlist.partners_of(int(c_idx)):
        row = A.loc[A['globalIdx'] == int(nbr)].iloc[0]
        if row['resName'] == bridge_residue and row['atomName'] in bridge_oxygens:
            return int(nbr)
    return None


def _ring_matching(ring_carbons, ring_nitrogens, bondlist):
    """Choose a perfect C-N matching in the 6-cycle ring.

    For a 1,3,5-triazine ring with carbons C1, C2, C3 and nitrogens
    N1, N2, N3 (in ring-traversal order ``C1-N1-C2-N2-C3-N3-C1``), both
    possible matchings are valid; we pick the one that pairs each C with
    its predecessor N in traversal order, i.e. ``Ci`` with the N that
    comes immediately before it.

    Args:
        ring_carbons (list[int]): C atom indices in ring traversal order.
        ring_nitrogens (list[int]): N atom indices in ring traversal order,
            interleaved with the carbons (N1 is between C1 and C2, etc.).
        bondlist: TC.Topology.bondlist for sanity-checking adjacency.

    Returns:
        list[tuple[int, int]]: list of (c_idx, n_idx) matched pairs.
    """
    matching = []
    for i, c in enumerate(ring_carbons):
        n = ring_nitrogens[(i - 1) % len(ring_nitrogens)]
        assert bondlist.are_bonded(c, n), (
            f'ring_matching: expected C{i+1}={c} bonded to N={n} '
            f'but they are not adjacent'
        )
        matching.append((int(c), int(n)))
    return matching


def _severed_ring_bonds(ring_carbons, ring_nitrogens, matching):
    """The 3 ring N-C bonds that are NOT in the matching."""
    kept = {(min(c, n), max(c, n)) for c, n in matching}
    severed = []
    # Ring edges in traversal: C1-N1, N1-C2, C2-N2, N2-C3, C3-N3, N3-C1
    n_ring = len(ring_carbons)
    for i in range(n_ring):
        c = ring_carbons[i]
        # the two N's adjacent to this C
        n_prev = ring_nitrogens[(i - 1) % n_ring]
        n_next = ring_nitrogens[i]
        for n in (n_prev, n_next):
            pair = (min(c, n), max(c, n))
            if pair not in kept and pair not in {(a, b) for a, b in severed}:
                severed.append((c, n))
    return severed


# ---------------------------------------------------------------------------
# Geometry placement
# ---------------------------------------------------------------------------


def _place_cyn_along(bpa_o_xyz, bpa_o_h_xyz, oc_len=0.136, cn_len=0.116):
    """Compute initial placement for a transferred -C#N group along the
    direction from BPA-O toward the (about-to-be-deleted) H, so the new
    cap sticks out along the original O-H bond.

    This is the unconditioned placement: it is where the cap goes if that
    direction happens to be empty, and it is blind to whether it is.  See
    :func:`_choose_cap_placement`, which uses it as its first candidate.

    Lengths in nm: standard -O-C ester ~1.36 A, -C#N ~1.16 A.
    Returns (C_xyz, N_xyz).
    """
    o = np.asarray(bpa_o_xyz, dtype=float)
    h = np.asarray(bpa_o_h_xyz, dtype=float)
    direction = h - o
    norm = np.linalg.norm(direction)
    if norm < 1e-8:
        direction = np.array([1.0, 0.0, 0.0])
        norm = 1.0
    u = direction / norm
    c_xyz = o + u * oc_len
    n_xyz = c_xyz + u * cn_len
    return c_xyz, n_xyz


_PLACEMENT_LOCALITY = 0.8
"""nm: an atom farther than this from the target oxygen cannot constrain a cap
that reaches only ~0.25 nm from it, so it is not scored against."""


def _report_placements(placements, target):
    """Say how well the transferred caps landed, and complain about the tight ones.

    The transfer count is the quantity that predicts whether this stage will
    survive, so it is stated plainly rather than left to be inferred, and a
    placement that could not reach the target clearance is named while it is
    still cheap to act on -- the alternative is reading it back out of a
    Lennard-Jones term at step 0.

    Args:
        placements (list): one dict per transferred cap, with 'o' and 'clearance'
        target (float): the clearance in nm that was being aimed for

    Returns:
        dict: n_transferred, n_direction_searched, n_below_target,
        n_preferred_out_of_angle, min_clearance_nm and the medians
    """
    if not placements:
        return {'n_transferred': 0, 'n_direction_searched': 0,
                'n_below_target': 0, 'n_preferred_out_of_angle': 0,
                'min_clearance_nm': None}
    gaps = np.array([p['clearance'] for p in placements], dtype=float)
    searched = sum(1 for p in placements if not p['kept_o_h_direction'])
    # why the preferred direction was abandoned matters as much as how often.
    # A metric that rejects every preferred direction for a reason unrelated
    # to crowding is what 2.6.0 shipped, and it read as a crowded box; this
    # separates "the site was occupied" from "the O-H vector was not a
    # placeable C-O-C angle", which would point at the angle window instead.
    out_of_window = sum(1 for p in placements if not p.get('blind_in_window', True))
    below = [p for p in placements if p['clearance'] < target]
    logger.info(
        f'triazine_to_cyanate_cap: placed {len(placements)} transferred caps '
        f'({searched} needed a direction search); clearance min/median '
        f'{gaps.min():.3f}/{float(np.median(gaps)):.3f} nm against a target of {target:.3f}'
    )
    blind = np.array([p['blind_clearance'] for p in placements
                      if p.get('blind_clearance') is not None], dtype=float)
    if blind.size:
        logger.info(
            f'triazine_to_cyanate_cap: along the O-H vector alone these same caps '
            f'would have had min/median {blind.min():.3f}/{float(np.median(blind)):.3f} nm, '
            f'with {int((blind < 0.10).sum())} under 0.10 nm -- that is what placement '
            f'did before the direction search existed, measured on this box'
        )
    if out_of_window:
        logger.info(
            f'triazine_to_cyanate_cap: {out_of_window} of {len(placements)} preferred '
            f'O-H directions were outside the placeable C-O-C angle window and were '
            f'searched away from for that reason rather than for crowding'
        )
    if below:
        worst = sorted(below, key=lambda p: p['clearance'])[:5]
        detail = ', '.join(f'O {int(p["o"])} at {p["clearance"]:.3f} nm' for p in worst)
        logger.warning(
            f'triazine_to_cyanate_cap: {len(below)} of {len(placements)} transferred '
            f'caps could not reach {target:.3f} nm of clearance in any direction '
            f'(worst: {detail}). The box is too crowded here for this many caps; '
            f'if the following minimization fails with a huge Lennard-Jones term, '
            f'this is why.'
        )
    # medians as well as minima: a single worst placement says whether this
    # build is in danger, but correlating placement quality against outcomes
    # across a series of runs needs a statistic that is not an extreme value
    return {
        'n_transferred': len(placements),
        'n_direction_searched': searched,
        'n_below_target': len(below),
        'n_preferred_out_of_angle': out_of_window,
        'min_clearance_nm': float(gaps.min()),
        'median_clearance_nm': float(np.median(gaps)),
        'blind_min_clearance_nm': float(blind.min()) if blind.size else None,
        'blind_median_clearance_nm': float(np.median(blind)) if blind.size else None,
        'n_blind_would_overlap': int((blind < 0.10).sum()) if blind.size else None,
    }


def _sphere_directions(n=48):
    """``n`` roughly uniform unit vectors on the sphere, deterministically.

    A Fibonacci spiral: no randomness, so a rebuild of the same system places
    its caps identically.

    Args:
        n (int): how many directions to generate

    Returns:
        numpy.ndarray: (n, 3) array of unit vectors
    """
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def _placement_neighbors(o_xyz, background, background_idx, box,
                         bonded=(), locality=_PLACEMENT_LOCALITY):
    """The atoms a cap hung on this oxygen actually has to stay clear of.

    Two filters.  Distance: an atom farther than ``locality`` from the oxygen
    cannot constrain a group that reaches only ~0.25 nm from it, so it is not
    worth scoring against.  Identity: an atom the cap is *bonded* through is
    not an obstacle, and leaving one in silently destroys the metric, because
    its distance to the cap is fixed by the bond geometry rather than by how
    crowded the site is.

    Both members of ``bonded`` matter and for the same reason.  The oxygen
    itself sits at exactly ``oc_len`` from the cap carbon in *every* candidate
    direction, so leaving it in pins :func:`_clearance` at 0.136 nm and the
    search cannot tell directions apart at all.  The oxygen's aryl carbon is
    softer but has the same character: with it in, the reachable clearance is
    bounded by the C-O-C geometry at 0.272 nm, so every comfortably-placed cap
    reports that constant instead of a measurement and the median saturates.
    A distance metric cannot do clash detection and bond-angle enforcement at
    once; the angle is enforced separately, in :func:`_choose_cap_placement`.

    Only *this* cap's bonded atoms are dropped.  Every other cap's attachment
    oxygen is a real atom in the way, as is every cap already placed -- those
    carry an index of -1 and are never masked out here.

    Args:
        o_xyz (numpy.ndarray): position of the oxygen this cap attaches to
        background (numpy.ndarray): (m, 3) candidate obstacle positions
        background_idx (numpy.ndarray): (m,) global indices, -1 for already-placed cap atoms
        bonded (iterable): global indices the cap is bonded through -- the oxygen and its aryl carbon
        box (numpy.ndarray): orthorhombic box lengths, for minimum image
        locality (float): nm; ignore anything farther than this from the oxygen

    Returns:
        numpy.ndarray: (k, 3) positions to pass to :func:`_choose_cap_placement`
    """
    if len(background) == 0:
        return background
    near = _min_image_dist(o_xyz, background, box) < locality
    drop = {int(i) for i in bonded if i is not None}
    if drop:
        near &= ~np.isin(np.asarray(background_idx, dtype=int), list(drop))
    return background[near]


def _clearance(c_xyz, n_xyz, neighbors, box):
    """Smallest distance from either cap atom to any neighbouring atom.

    Args:
        c_xyz (numpy.ndarray): position of the cap carbon
        n_xyz (numpy.ndarray): position of the cap nitrogen
        neighbors (numpy.ndarray): (m, 3) positions to stay clear of
        box (numpy.ndarray): orthorhombic box lengths, for minimum image

    Returns:
        float: the smallest such distance in nm, or inf if there are no neighbours
    """
    if neighbors is None or len(neighbors) == 0:
        return float('inf')
    return float(min(_min_image_dist(c_xyz, neighbors, box).min(),
                     _min_image_dist(n_xyz, neighbors, box).min()))


_MIN_COC_ANGLE = 90.0
_MAX_COC_ANGLE = 150.0
"""degrees: the window of C-O-C angles a cap may be placed at.  An aryl ether
sits near 120; this is deliberately generous, because placement is a starting
guess and the angle potential does the rest.  What it excludes is the two ends
that are qualitatively wrong -- a cap folded back onto the ring it hangs off,
and a linear C-O-C, which is what a pure clearance search converges to since
antiparallel maximizes distance from the rest of the molecule."""


def _angle_window(directions, o_xyz, anchor_xyz,
                  min_angle=_MIN_COC_ANGLE, max_angle=_MAX_COC_ANGLE):
    """Boolean mask of ``directions`` whose C-O-C angle is inside the window.

    Args:
        directions (numpy.ndarray): (n, 3) unit vectors from the oxygen
        o_xyz (numpy.ndarray): the oxygen's position
        anchor_xyz (numpy.ndarray): the aryl carbon's position, or None for no constraint
        min_angle (float): degrees; below this the cap is folding onto the ring
        max_angle (float): degrees; above this the ether is unphysically extended

    Returns:
        numpy.ndarray: (n,) boolean, all True when there is no anchor
    """
    d = np.atleast_2d(np.asarray(directions, dtype=float))
    if anchor_xyz is None:
        return np.ones(len(d), dtype=bool)
    a = np.asarray(anchor_xyz, dtype=float) - np.asarray(o_xyz, dtype=float)
    na = np.linalg.norm(a)
    if na < 1e-8:
        return np.ones(len(d), dtype=bool)
    cosines = d @ (a / na)
    return ((cosines <= np.cos(np.radians(min_angle))) &
            (cosines >= np.cos(np.radians(max_angle))))


def _choose_cap_placement(o_xyz, ref_xyz, neighbors, box, oc_len=0.136,
                          cn_len=0.116, target=0.15, n_directions=48,
                          anchor_xyz=None):
    """Place a transferred -C#N group in the clearest direction available.

    The old O-H vector is tried first and kept if it is clear enough, so a
    cap in open space lands exactly where it always did.  Only when that
    direction is occupied does this search the sphere -- which is the case
    that used to produce a step-0 Lennard-Jones term around 1e15 and a
    minimization that could not recover.  Bond lengths are held fixed; only
    the direction moves, so the chemistry is unchanged.

    The search is over directions inside a C-O-C angle window rather than the
    whole sphere.  Clearance alone is maximized by a linear ether, since
    antiparallel puts the cap as far as possible from the ring it hangs off,
    so a pure clearance search drifts toward a geometry no aryl ether adopts.
    The window says that explicitly instead of leaving it to a distance metric
    that cannot express it; see :data:`_MIN_COC_ANGLE`.

    Args:
        o_xyz (numpy.ndarray): position of the bridge oxygen the cap attaches to
        ref_xyz (numpy.ndarray): position defining the preferred direction, normally the O's about-to-be-deleted H
        neighbors (numpy.ndarray): (m, 3) positions the cap must avoid
        box (numpy.ndarray): orthorhombic box lengths, for minimum image
        oc_len (float): O-C bond length in nm
        cn_len (float): C#N bond length in nm
        target (float): clearance in nm that counts as good enough to stop searching; the default 0.15 is about what a steepest-descent minimization absorbs without trouble.  ``neighbors`` must already have the cap's own bonded atoms removed (see :func:`_placement_neighbors`) or this number is meaningless -- with the attachment oxygen left in, every direction scores exactly ``oc_len``
        n_directions (int): how many directions to try when the preferred one is occupied
        anchor_xyz (numpy.ndarray): position of the oxygen's aryl carbon, defining the C-O-C angle; None disables the angular guard

    Returns:
        dict: 'c' and 'n' positions, 'clearance' achieved, 'kept_o_h_direction',
        and 'blind_clearance' -- what the O-H direction alone would have given.
        The last is free to compute, since that direction is tried first, and it
        is the only way to measure on a real box what the old blind placement
        was actually doing.
    """
    o = np.asarray(o_xyz, dtype=float)
    c_xyz, n_xyz = _place_cyn_along(o, ref_xyz, oc_len=oc_len, cn_len=cn_len)
    blind = _clearance(c_xyz, n_xyz, neighbors, box)
    # the O-H direction is kept only if it is both clear enough and a sane
    # C-O-C angle; for a real hydroxyl it always is the latter, but the
    # fallback direction used when the O has no H need not be
    blind_ok = bool(_angle_window((c_xyz - o) / oc_len, o, anchor_xyz)[0])
    if blind_ok and blind >= target:
        return {'c': c_xyz, 'n': n_xyz, 'clearance': blind,
                'kept_o_h_direction': True, 'blind_clearance': blind,
                'blind_in_window': True}
    directions = _sphere_directions(n_directions)
    allowed = directions[_angle_window(directions, o, anchor_xyz)]
    # a window this wide always leaves a good fraction of the sphere, but a
    # caller narrowing it should still get a placement rather than nothing
    if len(allowed) == 0:
        allowed = directions
    best = blind if blind_ok else -1.0
    if not blind_ok:
        c_xyz, n_xyz = None, None
    for u in allowed:
        cc = o + u * oc_len
        nn = cc + u * cn_len
        gap = _clearance(cc, nn, neighbors, box)
        if gap > best:
            best, c_xyz, n_xyz = gap, cc, nn
            if best >= target:
                break
    return {'c': c_xyz, 'n': n_xyz, 'clearance': best,
            'kept_o_h_direction': False, 'blind_clearance': blind,
            'blind_in_window': blind_ok}


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def _completion_stats(residue_name, n_crosslinkers, n_incomplete, log=True, placement=None,
                      bond_histogram=None):
    """Summarize crosslinker completion and log it.

    A crosslinker survives this repair only if every one of its sites is
    filled; the incomplete ones are dismantled back into caps.  So the
    fraction of crosslinkers that survive intact -- not the bond conversion
    the cure iterates against -- is the quantity an experiment measures.
    For a cyanate ester that is the FTIR -OCN conversion, because each
    complete triazine consumes exactly three -OCN groups.

    Args:
        residue_name (str): crosslinker residue name, used only in the message
        n_crosslinkers (int): number of crosslinker residues in the system
        n_incomplete (int): how many of them are being dismantled
        log (bool): emit the summary message; False when there is nothing to summarize
        placement (dict): cap-placement summary from :func:`_report_placements`, merged into the result
        bond_histogram (dict): how many crosslinkers carried 0, 1, ... bonds
            *before* repair dismantled any of them, keyed by bond count.  This
            is the distribution ``n_complete`` is the top bin of, and it is the
            only record of it: repair rewrites the topology, and the final
            structure does not say which cap came from which ring.  Without it
            the independence assumption behind
            ``crosslinker conversion = bond conversion ** f`` can only be
            tested through its ``n_complete`` consequence.  Zero-filled
            through the full bond count, so the top bin is always present
            and always equals ``n_complete``; that is checked here.

    Returns:
        dict: n_crosslinkers, n_complete, n_dismantled, crosslinker_conversion,
        and prerepair_bond_counts when a histogram was supplied
    """
    n_complete = n_crosslinkers - n_incomplete
    chi = n_complete / n_crosslinkers if n_crosslinkers else 0.0
    if log:
        logger.info(
            f'triazine_to_cyanate_cap: {n_complete} of {n_crosslinkers} '
            f'{residue_name} residues are complete and survive repair; '
            f'crosslinker conversion {chi:.3f}'
        )
    out = {
        'residue': residue_name,
        'n_crosslinkers': n_crosslinkers,
        'n_complete': n_complete,
        'n_dismantled': n_incomplete,
        'crosslinker_conversion': chi,
    }
    if bond_histogram is not None:
        out['prerepair_bond_counts'] = {int(k): int(v) for k, v in sorted(bond_histogram.items())}
        # the histogram and the completion count are computed independently, so
        # they are each other's check: the top bin IS the surviving population
        binned = sum(out['prerepair_bond_counts'].values())
        top = out['prerepair_bond_counts'].get(max(out['prerepair_bond_counts'], default=0), 0)
        if binned != n_crosslinkers or top != n_complete:
            logger.warning(
                f'triazine_to_cyanate_cap: pre-repair histogram disagrees with the '
                f'completion count ({binned} binned vs {n_crosslinkers} residues, '
                f'{top} in the top bin vs {n_complete} complete); one of the two '
                f'is wrong and the reported conversion cannot be trusted'
            )
        if log:
            spread = ', '.join(f'{k}:{v}' for k, v in out['prerepair_bond_counts'].items())
            logger.info(
                f'triazine_to_cyanate_cap: pre-repair bonds per {residue_name} -- {spread}'
            )
    if placement:
        out.update(placement)
    return out


def triazine_to_cyanate_cap(TC, moldict, spec, reactions):
    """Execute the triazine-to-cyanate-cap repair.

    Args:
        TC: TopoCoord (modified in place).
        moldict: MoleculeDict including the linked-product template that
            describes a fully repaired BPA-O-C#N cap (looked up by the
            ``cap_template`` name in spec, with symmetry siblings auto-tried).
        spec: dict from postcure_repair config (see example 6 YAML).
        reactions: ReactionList (currently unused; reserved for cross-checks).

    Returns:
        dict: completion statistics as returned by :func:`_completion_stats`.
    """
    cl = spec['crosslinker']
    br = spec['bridge']
    cap_residue_name = spec['cap_residue']
    cap_template = spec['cap_template']
    full_bond_count = int(cl.get('full_bond_count', 3))
    ring_c_names = cl['ring_carbon_atoms']      # e.g. ['C1', 'C2', 'C3']
    ring_n_names = cl['ring_nitrogen_atoms']    # e.g. ['N1', 'N2', 'N3']
    bridge_residue = br['residue']
    bridge_o_names = br['reactive_oxygen_atoms']  # e.g. ['O1', 'O2']
    search_radius = float(spec.get('cap_search_radius', 0.6))

    bondlist = TC.Topology.bondlist

    # ---- Phase 1: find incomplete crosslinker residues and plan caps ----
    cl_residues = _find_crosslinker_residues(TC, cl['residue'])
    if not cl_residues:
        logger.info(f'triazine_to_cyanate_cap: no {cl["residue"]} residues found; nothing to do')
        return _completion_stats(cl['residue'], 0, 0, log=False,
                                 bond_histogram={n: 0 for n in range(full_bond_count + 1)})

    incomplete_plans = []     # list of dicts, one per incomplete triazine
    h_to_delete = set()       # global H indices, batched for final deletion
    # how many bonds each crosslinker carries before anything is dismantled;
    # repair destroys this distribution, so it has to be counted here or not
    # at all.  Zero-filled so the reported shape does not depend on the box.
    bond_histogram = {n: 0 for n in range(full_bond_count + 1)}

    for resnum, atoms_by_name in cl_residues.items():
        ring_c = [atoms_by_name[n] for n in ring_c_names]
        ring_n = [atoms_by_name[n] for n in ring_n_names]
        bonded_o = [_bonded_bridge_o(TC, c, bridge_residue, bridge_o_names) for c in ring_c]
        k = sum(1 for o in bonded_o if o is not None)
        bond_histogram[k] = bond_histogram.get(k, 0) + 1
        if k >= full_bond_count:
            continue
        # plan dismantle
        matching = _ring_matching(ring_c, ring_n, bondlist)
        severed = _severed_ring_bonds(ring_c, ring_n, matching)
        # plan H deletions: each dangling ring C has one C-H to remove
        for c_idx, o_partner in zip(ring_c, bonded_o):
            if o_partner is None:
                for nbr in bondlist.partners_of(c_idx):
                    name = ts.get_attr(TC, nbr, 'atomName')
                    if name.upper().startswith('H'):
                        h_to_delete.add(int(nbr))
        # one cap per (c, n) pair: tag in-place vs free
        caps = []
        c_to_o = dict(zip(ring_c, bonded_o))
        for c, n in matching:
            caps.append({
                'c': c,
                'n': n,
                'bonded_o': c_to_o[c],   # None for free caps
                'taz_resnum': resnum,
            })
        incomplete_plans.append({'resnum': resnum, 'caps': caps, 'severed': severed})

    if not incomplete_plans:
        logger.info(
            f'triazine_to_cyanate_cap: all {len(cl_residues)} {cl["residue"]} '
            f'residues are fully bonded; nothing to do'
        )
        return _completion_stats(cl['residue'], len(cl_residues), 0,
                                 bond_histogram=bond_histogram)

    total_caps = sum(len(p['caps']) for p in incomplete_plans)
    free_caps_planned = sum(1 for p in incomplete_plans for c in p['caps'] if c['bonded_o'] is None)
    logger.info(
        f'triazine_to_cyanate_cap: {len(incomplete_plans)} incomplete '
        f'{cl["residue"]} residues identified ({total_caps} caps total, '
        f'{free_caps_planned} free fragments to donate)'
    )

    # ---- Phase 2: greedy-match free caps to unreacted bridge-O's ----
    free_caps = [c for p in incomplete_plans for c in p['caps'] if c['bonded_o'] is None]
    unreacted_o = _find_unreacted_bridge_oxygens(TC, bridge_residue, bridge_o_names)
    if len(unreacted_o) != len(free_caps):
        logger.warning(
            f'triazine_to_cyanate_cap: atom-conservation mismatch — '
            f'{len(free_caps)} free fragments vs {len(unreacted_o)} unreacted '
            f'{bridge_residue}-O sites.  Excess fragments will be discarded '
            f'(their atoms will be deleted).'
        )
    donations = _greedy_match(TC, free_caps, unreacted_o, search_radius)

    # ---- Phase 3: execute surgery (atom indices stable through phase 5) ----
    # 3a. delete all severed ring N-C bonds (cascades angles/dihedrals)
    all_severed = [b for p in incomplete_plans for b in p['severed']]
    ts.delete_bonds(TC, all_severed)

    # 3b. for in-place caps: also delete the pre-existing BPA-O-C bond so
    # add_bonds_with_template can re-form it with CYN-context parameters.
    in_place_bonds_to_redo = []
    for p in incomplete_plans:
        for cap in p['caps']:
            if cap['bonded_o'] is not None:
                in_place_bonds_to_redo.append((cap['bonded_o'], cap['c']))
    ts.delete_bonds(TC, in_place_bonds_to_redo)

    # 3c. re-tag each (c, n) pair into its own CYN residue
    new_resnum = ts.next_free_resnum(TC)
    cap_records = []  # list of (cap_dict, bridge_o_idx, new_resnum)
    for p in incomplete_plans:
        for cap in p['caps']:
            bridge_o = donations.get(id(cap))  # may be None for the in-place case
            if cap['bonded_o'] is not None:
                target_o = cap['bonded_o']
            else:
                target_o = bridge_o  # may still be None if no donation found
            if target_o is None:
                # no home for this fragment — schedule its atoms for deletion
                h_to_delete.add(int(cap['c']))
                h_to_delete.add(int(cap['n']))
                continue
            ts.reassign_residue(TC, [cap['c'], cap['n']], cap_residue_name, new_resnum)
            ts.set_atom_attributes(TC, cap['c'], atomName='C1')
            ts.set_atom_attributes(TC, cap['n'], atomName='N1')
            cap_records.append({'c': cap['c'], 'n': cap['n'], 'o': target_o,
                                'resnum': new_resnum, 'free': cap['bonded_o'] is None})
            new_resnum += 1

    # 3d. for free-cap donations: relocate the CYN atoms next to their target O
    #
    # This is the phase that fails catastrophically when it fails at all.  A
    # cap dropped blindly along the old O-H vector can land on top of a
    # neighbour, and the resulting step-0 Lennard-Jones term -- of order 1e15
    # kJ/mol -- is not something the following minimization recovers from.
    # The count of caps to place here is an identity, (total reactive sites -
    # bonds formed), so it rises as conversion falls: a low-conversion build
    # places hundreds and is where this has been seen to kill runs.
    #
    # So each cap is placed against its local neighbourhood, and every cap
    # already placed is part of the neighbourhood the next one sees.
    free_records = [rec for rec in cap_records if rec['free']]
    placements = []
    placement_summary = None
    if free_records:
        box = TC.Coordinates.box.diagonal()
        target_clearance = float(spec.get('cap_min_clearance', 0.15))
        # the O's H is consumed in forming the cap; collect them all up front
        # so that a cap placed early does not dodge an atom about to vanish
        for rec in free_records:
            h = _find_bonded_hydrogen(TC, rec['o'])
            rec['o_h'] = h
            # the aryl carbon on the far side of the O: excluded from the
            # clearance metric, since its distance is set by bond geometry,
            # and used instead as the anchor for the C-O-C angle
            rec['o_anchor'] = _find_bonded_heavy(TC, rec['o'])
            if h is not None:
                h_to_delete.add(int(h))
        moving = {int(rec['c']) for rec in free_records} | {int(rec['n']) for rec in free_records}
        A = TC.Coordinates.A
        going = moving | {int(x) for x in h_to_delete}
        staying = A.loc[~A['globalIdx'].isin(going),
                        ['globalIdx', 'posX', 'posY', 'posZ']]
        # indices travel alongside the positions so each cap can drop the
        # atoms it is bonded through; placed caps get -1 and are never dropped
        background_idx = staying['globalIdx'].to_numpy(dtype=int)
        background = staying[['posX', 'posY', 'posZ']].to_numpy(dtype=float)
        for rec in free_records:
            o_xyz = ts.positions(TC, [rec['o']])[0]
            if rec['o_h'] is None:
                # no H to take a direction from (shouldn't happen for a
                # still-reactive O); any direction will do as a starting guess
                ref_xyz = o_xyz + np.array([1.0, 0.0, 0.0])
            else:
                ref_xyz = ts.positions(TC, [rec['o_h']])[0]
            local = _placement_neighbors(o_xyz, background, background_idx, box,
                                         bonded=(rec['o'], rec['o_anchor']))
            anchor_xyz = (ts.positions(TC, [rec['o_anchor']])[0]
                          if rec['o_anchor'] is not None else None)
            placed = _choose_cap_placement(o_xyz, ref_xyz, local, box,
                                           target=target_clearance,
                                           anchor_xyz=anchor_xyz)
            c_xyz, n_xyz = placed['c'], placed['n']
            ts.set_atom_attributes(TC, rec['c'],
                                   posX=float(c_xyz[0]), posY=float(c_xyz[1]), posZ=float(c_xyz[2]))
            ts.set_atom_attributes(TC, rec['n'],
                                   posX=float(n_xyz[0]), posY=float(n_xyz[1]), posZ=float(n_xyz[2]))
            background = np.vstack([background, c_xyz, n_xyz])
            background_idx = np.concatenate([background_idx, [-1, -1]])
            placements.append({'o': rec['o'], **{k: v for k, v in placed.items() if k not in ('c', 'n')}})
        placement_summary = _report_placements(placements, target_clearance)

    # ---- Phase 4: splice template params (forms BPA-O-C bond with full
    # angle/dihedral/pair entries from the cap_template molecule) ----
    pairs_to_add = [(int(rec['o']), int(rec['c']), 1) for rec in cap_records]
    # Charge is settled once, per molecule, after phase 5's deletions.  The
    # splice's own system-wide adjustment would run while the sacrificial H's
    # still exist, and phase 5 would then have to undo it on other atoms.
    touched = set()
    if pairs_to_add:
        touched |= {int(a) for a in ts.add_bonds_with_template(TC, pairs_to_add, moldict, cap_template,
                                                               adjust_charges=False)}

    # Refresh the C-N triple-bond parameters: map_from_templates updated the
    # CYN atom types (CA -> c1, NB -> n1) but only re-resolved the bond it
    # was actively mapping (BPA.O - CYN.C).  The C-N bond carries stale
    # aromatic CA-NB override parameters until we reset it.
    ts.refresh_bond_params(TC, [(rec['c'], rec['n']) for rec in cap_records])

    # zero out reactivity on every atom touched by the repair
    for rec in cap_records:
        for idx in (rec['o'], rec['c'], rec['n']):
            TC.set_gro_attribute_by_attributes('z', 0, {'globalIdx': int(idx)})
            TC.set_gro_attribute_by_attributes('nreactions', 1, {'globalIdx': int(idx)})

    # ---- Phase 5: batched H + orphan-atom deletion (reindexes) ----
    if h_to_delete:
        logger.debug(f'triazine_to_cyanate_cap: deleting {len(h_to_delete)} sacrificial/orphan atoms')
        # the heavy atoms that lose a bonded H are touched too; collect them
        # before delete_atoms reindexes
        for d_idx in h_to_delete:
            for nbr in TC.Topology.bondlist.partners_of(int(d_idx)):
                touched.add(int(nbr))
        touched -= {int(x) for x in h_to_delete}
        idx_mapper = TC.delete_atoms(sorted(h_to_delete))
        touched = {idx_mapper[a] for a in touched if a in idx_mapper}
    # Each repaired molecule must end neutral, not just the system.  A single
    # system-wide rebalance left every molecule of a zero-conversion melt
    # carrying +0.21 or -0.27 e, depending on which of its caps were formed
    # in place and which were transferred.
    neutral = ts.neutralize_touched_fragments(TC, touched)
    if neutral['n_fragments']:
        logger.info(f'triazine_to_cyanate_cap: neutralized {neutral["n_fragments"]} repaired molecule(s); '
                    f'largest excess before correction {neutral["max_excess"]:.4f} e')

    return _completion_stats(cl['residue'], len(cl_residues), len(incomplete_plans),
                             placement=placement_summary,
                             bond_histogram=bond_histogram)


# ---------------------------------------------------------------------------
# Helpers for matching free caps to unreacted bridge-O's
# ---------------------------------------------------------------------------


def _find_unreacted_bridge_oxygens(TC, bridge_residue, bridge_o_names):
    """Return list of (globalIdx, np.ndarray xyz) for every bridge-O atom
    that still carries a phenolic H — i.e. did not bond to any
    crosslinker during cure."""
    A = TC.Coordinates.A
    sub = A[(A['resName'] == bridge_residue) & (A['atomName'].isin(bridge_o_names))]
    unreacted = []
    for r in sub.itertuples():
        if _find_bonded_hydrogen(TC, int(r.globalIdx)) is not None:
            unreacted.append((int(r.globalIdx),
                              np.array([r.posX, r.posY, r.posZ])))
    return unreacted


def _find_bonded_hydrogen(TC, idx):
    """Return globalIdx of an H bonded to ``idx``, or None."""
    for nbr in TC.Topology.bondlist.partners_of(int(idx)):
        name = ts.get_attr(TC, nbr, 'atomName')
        if str(name).upper().startswith('H'):
            return int(nbr)
    return None


def _find_bonded_heavy(TC, idx):
    """Return globalIdx of the first non-hydrogen atom bonded to ``idx``, or None.

    For a bridge oxygen about to take a cap this is its aryl carbon, which is
    the anchor the C-O-C angle is measured against.
    """
    for nbr in TC.Topology.bondlist.partners_of(int(idx)):
        name = ts.get_attr(TC, nbr, 'atomName')
        if not str(name).upper().startswith('H'):
            return int(nbr)
    return None


def _greedy_match(TC, free_caps, unreacted_o, search_radius):
    """Match each free cap to its nearest available unreacted bridge-O.

    Iterates free caps in order, each grabs the closest still-unmatched
    bridge-O within ``search_radius``; on a miss, the radius is expanded
    by a factor of 2 each retry up to 10x; final fallback is the globally
    nearest still-unmatched bridge-O.

    Returns dict keyed by id(cap dict) -> bridge-O globalIdx.
    """
    if not free_caps or not unreacted_o:
        return {}
    box = TC.Coordinates.box.diagonal()
    cap_positions = ts.positions(TC, [c['c'] for c in free_caps])
    o_positions = np.array([xyz for _, xyz in unreacted_o])
    o_indices = [i for i, _ in unreacted_o]
    o_available = np.ones(len(o_indices), dtype=bool)

    out = {}
    for cap, c_xyz in zip(free_caps, cap_positions):
        d = _min_image_dist(c_xyz, o_positions, box)
        d = np.where(o_available, d, np.inf)
        # progressive radius expansion
        chosen = None
        radius = search_radius
        for _ in range(10):
            idx = int(np.argmin(d))
            if np.isfinite(d[idx]) and d[idx] <= radius:
                chosen = idx
                break
            radius *= 2.0
        if chosen is None:
            # global fallback
            chosen = int(np.argmin(d))
            if not np.isfinite(d[chosen]):
                logger.warning(
                    f'_greedy_match: cap c={cap["c"]} has no available bridge-O '
                    f'remaining; will be discarded'
                )
                continue
        o_available[chosen] = False
        out[id(cap)] = int(o_indices[chosen])
    return out


def _min_image_dist(p, qs, box):
    """Minimum-image distance from p to each row of qs in an orthorhombic box."""
    delta = qs - p
    delta -= box * np.round(delta / box)
    return np.linalg.norm(delta, axis=1)
