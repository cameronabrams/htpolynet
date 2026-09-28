"""Finds bonds threaded through rings in a finished network.

htpolynet can refuse to *create* a threading during a ring cure, and until now it could
not *tell you whether one is there*.  That asymmetry is the reason this exists: a
preventive check whose effect cannot be measured with the same tool is hard to trust,
and the numbers that justified those filters all came from scripts outside the package.

Nothing here is specific to cyanate esters.  Any network built by dragging monomers
together can thread a bond through a ring, and htpolynet builds those by construction.

The geometry lives in :mod:`htpolynet.geometry.piercing`, shared with the ring cure's own
filters so a measurement and the thing it measures cannot drift apart.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

import networkx as nx
import numpy as np

from scipy.spatial import cKDTree

from ..core.topology import _element_of
from ..geometry.piercing import pierces_ring, ring_frame

logger = logging.getLogger(__name__)


def detect_rings(bonds, max_size=8):
    """Smallest ring through each bond, for a topology with no ring list of its own.

    A build made by htpolynet carries its rings already; this is the fallback for a bare
    top/gro pair.  Enumerating all chordless cycles is not an option --- on a percolated
    network it does not finish in reasonable time --- so this takes the shortest
    alternative path between each bonded pair, which gives the smallest ring through
    that bond and nothing larger.

    Args:
        bonds (pandas.DataFrame): with ai and aj columns
        max_size (int): largest ring to look for

    Returns:
        list: rings, each a list of atom indices in cyclic order
    """
    G = nx.Graph()
    G.add_edges_from((int(r.ai), int(r.aj)) for r in bonds.itertuples())
    seen, out = set(), []
    for u, v in list(G.edges()):
        G.remove_edge(u, v)
        try:
            path = nx.shortest_path(G, u, v)
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            path = None
        G.add_edge(u, v)
        if path is None or len(path) > max_size:
            continue
        key = frozenset(path)
        if key in seen:
            continue
        seen.add(key)
        out.append([int(x) for x in path])
    return out


def ring_composition(ring, element_of):
    """A ring's composition as a formula string, heaviest-count first.

    Deliberately general: a triazine reads ``C3N3`` and a benzene ``C6`` without this
    module knowing either name.  Classify on the whole ring --- keying on residue count
    or on a truncated printout are both ways this has been got wrong.
    """
    counts = {}
    for i in ring:
        e = element_of.get(int(i), '?')
        counts[e] = counts.get(e, 0) + 1
    return ''.join(f'{e}{counts[e]}' for e in sorted(counts, key=lambda x: (-counts[x], x)))


def find_piercings(TC, max_ring=8, margin=0.3):
    """Every bond that passes through a ring.

    Args:
        TC (TopoCoord): the system, with coordinates loaded
        max_ring (int): largest ring to consider when detecting rings
        margin (float): how far beyond a ring's radius to look for bonds, in nm

    Returns:
        list: one dict per piercing, with the ring, its composition, the bond, whether
            the bond belongs to a residue of the ring, and the geometry of the crossing
    """
    A = TC.Coordinates.A
    positions = {int(r.globalIdx): np.array([r.posX, r.posY, r.posZ]) for r in A.itertuples()}
    box = np.asarray(TC.Coordinates.box.diagonal(), dtype=float)
    at = TC.Topology.D['atoms']
    element_of = {int(r.nr): _element_of(r.atom, r.type) for r in at.itertuples()}
    resnr_of = {int(r.nr): int(r.resnr) for r in at.itertuples()}
    bonds = TC.Topology.D['bonds']

    rings = [[int(x) for x in r.idx] for r in (getattr(TC.Topology, 'rings', None) or [])]
    rings = [r for r in rings if 3 <= len(r) <= max_ring]
    if not rings:
        logger.info('No ring list on this topology; detecting rings from the bond graph')
        rings = detect_rings(bonds, max_size=max_ring)
    logger.info(f'Scanning {len(rings)} ring(s) against {bonds.shape[0]} bond(s)')

    ids = np.fromiter(positions.keys(), dtype=int, count=len(positions))
    pts = np.array([positions[i] for i in ids], dtype=float)
    tree = cKDTree(np.mod(pts, box), boxsize=box)
    bonds_of = {}
    for r in bonds.itertuples():
        ai, aj = int(r.ai), int(r.aj)
        bonds_of.setdefault(ai, []).append((ai, aj))
        bonds_of.setdefault(aj, []).append((ai, aj))

    out = []
    for ring in rings:
        frame = ring_frame(ring, positions, box)
        if frame is None:
            continue
        _, O, _, radius = frame
        near = tree.query_ball_point(np.mod(O, box), radius + margin)
        ring_set, seen = set(ring), set()
        for j in near:
            a = int(ids[j])
            if a in ring_set:
                continue
            for bond in bonds_of.get(a, ()):
                key = (min(bond), max(bond))
                if key in seen:
                    continue
                seen.add(key)
                hit = pierces_ring(ring, bond, positions, box)
                if hit is None:
                    continue
                ring_res = {resnr_of.get(i) for i in ring}
                out.append({'ring': ring, 'composition': ring_composition(ring, element_of),
                            'bond': key, 'self': resnr_of.get(key[0]) in ring_res
                                              or resnr_of.get(key[1]) in ring_res,
                            **{k: v for k, v in hit.items() if k != 'point'}})
    return out


def format_report(piercings, overlong=0.2):
    """A human-readable summary, correlating piercings with over-long bonds."""
    if not piercings:
        return 'No bond passes through a ring in this system.'
    lines = [f'{len(piercings)} bond(s) pass through a ring:', '']
    by_comp = {}
    for p in piercings:
        by_comp[p['composition']] = by_comp.get(p['composition'], 0) + 1
    lines.append('  by ring composition: '
                 + ', '.join(f'{k} {v}' for k, v in sorted(by_comp.items())))
    n_self = sum(1 for p in piercings if p['self'])
    lines.append(f'  {n_self} by a bond of a residue in the ring, '
                 f'{len(piercings) - n_self} by a foreign residue')
    n_long = sum(1 for p in piercings if p['length'] > overlong)
    lines.append(f'  {n_long} of {len(piercings)} piercing bond(s) are longer than '
                 f'{overlong} nm, which is longer than any real bond')
    lines.append('')
    lines.append('  ring          composition  bond            offset  frac  angle  length')
    for p in sorted(piercings, key=lambda x: -x['length']):
        r = ','.join(str(x) for x in p['ring'][:3]) + '...'
        lines.append(f'  {r:<13} {p["composition"]:<12} '
                     f'{p["bond"][0]}-{p["bond"][1]:<10} '
                     f'{p["offset"] * 10:6.2f} {p["offset_frac"]:5.2f} '
                     f'{p["angle"]:5.1f} {p["length"] * 10:6.2f}   (A, deg)')
    return '\n'.join(lines)


def piercings(args):
    """CLI handler for ``htpolynet piercings``.

    Args:
        args (argparse.Namespace): parsed arguments with ``.top``, ``.gro``,
            ``.max_ring``, and optional ``.json``
    """
    import json
    import os

    from ..core.topocoord import TopoCoord

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    for f in (args.top, args.gro):
        if not os.path.exists(f):
            logger.error(f'{f} not found')
            raise SystemExit(1)
    try:
        TC = TopoCoord(topfilename=args.top, grofilename=args.gro)
        found = find_piercings(TC, max_ring=args.max_ring)
    except Exception as e:
        logger.error(f'Piercing analysis failed: {e}')
        raise SystemExit(1)
    print(format_report(found))
    if getattr(args, 'json', None):
        with open(args.json, 'w') as f:
            json.dump([{k: (list(v) if isinstance(v, tuple) else v)
                        for k, v in p.items()} for p in found], f, indent=2)
        logger.info(f'Wrote {args.json}')
    raise SystemExit(1 if found else 0)
