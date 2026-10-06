"""Non-covalent structure around rings in a finished network.

Three measurements on one top/gro pair, each normalised so that "more than chance" means
something:

- :func:`face_on_contacts` --- how often a chosen atom sits over the *face* of a ring,
  close and on the ring's normal, relative to a uniform random placement.  This is the
  geometry of a lone-pair/pi contact.
- :func:`centroid_rdf` --- orientation-averaged radial distribution from chosen atoms (or
  bridge groups) to ring centroids, excluding rings the selection is covalently part of.
- :func:`ring_flip_clearance` --- for each ring, the cylinder it must sweep to flip about
  its para axis and how many foreign heavy atoms sit inside it: the static steric part of
  a ring-flip barrier.

Rings are selected by composition (``C3N3`` for a triazine, ``C6`` for a phenylene) and
atoms by force-field type, so nothing here knows any chemistry by name.  Ported from
htpolynet-study's bisphenol-bridge analysis (2026-10-06), where sulfone oxygens were
found over triazine faces about twice as often as chance.

Single structure only: averaging over a trajectory means running this per frame.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging
import math
from collections import defaultdict

import numpy as np

from ..core.topology import _element_of
from ..geometry.piercing import ring_frame
from .piercings import detect_rings, ring_composition

logger = logging.getLogger(__name__)


class RingSystem:
    """Positions, box, per-atom data, bond graph and rings of one structure.

    Args:
        TC (TopoCoord): the system, with coordinates loaded
        max_ring (int): largest ring to detect when the topology carries no ring list
    """
    def __init__(self, TC, max_ring=6):
        A = TC.Coordinates.A
        self.box = np.asarray(TC.Coordinates.box.diagonal(), dtype=float)
        self.pos = {int(r.globalIdx): np.array([r.posX, r.posY, r.posZ], dtype=float)
                    for r in A.itertuples()}
        at = TC.Topology.D['atoms']
        self.type = {int(r.nr): str(r.type) for r in at.itertuples()}
        self.resnr = {int(r.nr): int(r.resnr) for r in at.itertuples()}
        self.mass = {int(r.nr): float(r.mass) for r in at.itertuples()}
        self.element = {int(r.nr): _element_of(r.atom, r.type) for r in at.itertuples()}
        self.adj = defaultdict(set)
        bonds = TC.Topology.D['bonds']
        for r in bonds.itertuples():
            ai, aj = int(r.ai), int(r.aj)
            self.adj[ai].add(aj)
            self.adj[aj].add(ai)
        self.rings = detect_rings(bonds, max_size=max_ring)
        self.composition = [ring_composition(r, self.element) for r in self.rings]
        self.volume = float(np.prod(self.box))

    def rings_of(self, composition):
        """Rings (atom lists in cyclic order) of the given composition, e.g. ``C3N3``."""
        return [r for r, c in zip(self.rings, self.composition) if c == composition]

    def atoms_of_types(self, types):
        """Atom indices whose force-field type is in `types`."""
        types = set(types)
        return sorted(i for i, t in self.type.items() if t in types)

    def heavy(self, i):
        return self.element.get(i, '?') != 'H'

    def mic(self, d):
        """Minimum-image of displacement(s) `d`."""
        return d - self.box * np.round(d / self.box)

    def unwrap(self, idx, anchor=None):
        """Positions of `idx` imaged next to `anchor` (default: the first of them)."""
        p = np.array([self.pos[i] for i in idx])
        a = p[0] if anchor is None else anchor
        return a + self.mic(p - a)

    def com(self, idx):
        p = self.unwrap(idx)
        m = np.array([self.mass[i] for i in idx])
        return (p * m[:, None]).sum(0) / m.sum()

    def within_bonds(self, seeds, maxd):
        """Atoms within `maxd` bonds of any seed, seeds included."""
        seen = set(seeds)
        frontier = list(seeds)
        for _ in range(maxd):
            nxt = []
            for a in frontier:
                for b in self.adj[a]:
                    if b not in seen:
                        seen.add(b)
                        nxt.append(b)
            frontier = nxt
        return seen


def face_on_contacts(S, donors, rings, rmax=0.37, max_angle=35.0, exclude_bonds=4,
                     rmin=0.0, exclude='bonds'):
    """Face-on donor-over-ring contacts, against a uniform random expectation.

    A contact is a donor between `rmin` and `rmax` of a ring centroid with the
    centroid-to-donor vector within `max_angle` of the ring normal (either face).  The
    expectation is the number of eligible donor-ring pairs times the fraction of the box
    the two acceptance cones occupy, ``(1 - cos max_angle) * (4/3) pi (rmax^3 - rmin^3) / V``.

    Rings covalently close to the donor are not eligible: with ``exclude='bonds'``, rings
    within `exclude_bonds` bonds; with ``exclude='residue'``, rings containing an atom of
    the donor's own residue.  Without an exclusion the result is meaningless: a donor
    scores contacts with its own substituents, which are two or three bonds away.

    Args:
        S (RingSystem): the structure
        donors (list): donor atom indices
        rings (list): rings, as atom lists in cyclic order
        rmax (float): contact distance, nm
        max_angle (float): largest angle from the ring normal, degrees
        exclude_bonds (int): covalent exclusion depth for ``exclude='bonds'``
        rmin (float): inner contact distance, nm; nonzero makes the region a shell
        exclude (str): ``'bonds'`` or ``'residue'``

    Returns:
        dict: observed, expected, enhancement (observed/expected), pairs (eligible
            donor-ring pairs) and contacts, an (n,2) array of (distance nm, angle deg)
            for every eligible pair inside 0.6 nm, for histogramming
    """
    frames = [ring_frame(r, S.pos, S.box) for r in rings]
    keep = [k for k, f in enumerate(frames) if f is not None]
    if not donors or not keep:
        return {'observed': 0, 'expected': 0.0, 'enhancement': float('nan'), 'pairs': 0,
                'contacts': np.empty((0, 2))}
    cen = np.array([frames[k][1] for k in keep])
    nor = np.array([frames[k][2] for k in keep])
    if exclude not in ('bonds', 'residue'):
        raise ValueError(f"exclude must be 'bonds' or 'residue', not {exclude!r}")
    ringsets = [set(rings[k]) for k in keep]
    ringres = [{S.resnr[i] for i in rings[k]} for k in keep]
    near_r = max(rmax, 0.6)
    out, npairs, obs = [], 0, 0
    for o in donors:
        if exclude == 'bonds':
            near = S.within_bonds([o], exclude_bonds)
            ok = np.array([not (near & rs) for rs in ringsets])
        else:
            ok = np.array([S.resnr[o] not in rr for rr in ringres])
        npairs += int(ok.sum())
        d = S.mic(cen[ok] - S.pos[o])
        r = np.linalg.norm(d, axis=1)
        m = r < near_r
        if not m.any():
            continue
        cosang = np.abs(np.einsum('ij,ij->i', d[m] / r[m][:, None], nor[ok][m]))
        ang = np.degrees(np.arccos(np.clip(cosang, 0.0, 1.0)))
        obs += int(((r[m] >= rmin) & (r[m] < rmax) & (ang < max_angle)).sum())
        out.extend(zip(r[m], ang))
    vcone = (1 - math.cos(math.radians(max_angle))) * (4.0 / 3.0) * math.pi * (rmax ** 3 - rmin ** 3)
    exp = npairs * vcone / S.volume
    return {'observed': obs, 'expected': exp,
            'enhancement': obs / exp if exp > 0 else float('nan'), 'pairs': npairs,
            'contacts': np.array(out) if out else np.empty((0, 2))}


def bridge_groups(S, composition):
    """Bridge groups between rings of one composition.

    A bridge centre is a heavy atom outside the rings bonded to atoms of two *different*
    rings of `composition`; its group is the centre plus its neighbours outside those
    rings (a sulfone's two oxygens, an isopropylidene's methyls).  For ``C6`` this picks
    -O-, -S-, -CH2-, -C(CH3)2-, -C(CF3)2- and -SO2- alike, and rejects an aryl ether
    oxygen that links one phenylene to a triazine.

    Returns:
        list: groups, each a list of atom indices with the centre first
    """
    ring_of = {}
    for k, r in enumerate(S.rings_of(composition)):
        for i in r:
            ring_of[i] = k
    groups = []
    for i in sorted(S.type):
        if i in ring_of or not S.heavy(i):
            continue
        if len({ring_of[j] for j in S.adj[i] if j in ring_of}) >= 2:
            groups.append([i] + sorted(j for j in S.adj[i] if j not in ring_of))
    return groups


def centroid_rdf(S, groups, rings, rmax=1.5, dr=0.02):
    """Radial distribution from group centres of mass to ring centres of mass.

    A ring containing any atom of the group's own residue is excluded: a monomer's
    bridge is covalently a fixed distance from rings it is part of, and counting those
    measures molecular geometry, not packing.

    Args:
        S (RingSystem): the structure
        groups (list): atom groups (a single atom is a group of one)
        rings (list): rings, as atom lists
        rmax (float): largest distance, nm
        dr (float): bin width, nm

    Returns:
        tuple: (bin centres in nm, g(r)); g is normalised by the eligible pair count
    """
    nbin = int(round(rmax / dr))
    edges = np.linspace(0.0, rmax, nbin + 1)
    if not groups or not rings:
        return 0.5 * (edges[1:] + edges[:-1]), np.zeros(nbin)
    gcom = np.array([S.com(g) for g in groups])
    gres = [{S.resnr[i] for i in g} for g in groups]
    rcom = np.array([S.com(r) for r in rings])
    rres = [{S.resnr[i] for i in r} for r in rings]
    hist = np.zeros(nbin)
    npair = 0
    for k in range(len(groups)):
        ok = np.array([not (gres[k] & s) for s in rres])
        npair += int(ok.sum())
        r = np.linalg.norm(S.mic(rcom[ok] - gcom[k]), axis=1)
        hist += np.histogram(r, bins=nbin, range=(0.0, rmax))[0]
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / (npair * shell / S.volume) if npair else hist
    return 0.5 * (edges[1:] + edges[:-1]), g


def flip_axis(S, ring):
    """The two ring atoms a ring flips about: the para pair carrying exocyclic heavy atoms.

    Falls back to the farthest-apart pair when a ring does not have exactly one para pair
    of substituted atoms.  The farthest pair alone is not a flip axis --- a hexagon has
    three para pairs at nearly equal separation --- so the fallback is logged.
    """
    n = len(ring)
    rs = set(ring)
    subst = [any(S.heavy(j) and j not in rs for j in S.adj[i]) for i in ring]
    if n % 2 == 0:
        pairs = [(k, k + n // 2) for k in range(n // 2) if subst[k] and subst[k + n // 2]]
        if len(pairs) == 1:
            return ring[pairs[0][0]], ring[pairs[0][1]], True
    p = S.unwrap(ring)
    dmat = np.linalg.norm(p[:, None] - p[None, :], axis=-1)
    i, j = np.unravel_index(np.argmax(dmat), dmat.shape)
    return ring[i], ring[j], False


def ring_flip_clearance(S, rings):
    """Room each ring needs to flip about its para axis, and what is in the way.

    The ring's own group is its atoms, their exocyclic neighbours that are not in another
    ring, and the hydrogens on those neighbours --- an ortho methyl's hydrogens are what
    actually sweep.  ``R_sweep`` is the largest distance of that group from the axis;
    ``blocking`` counts foreign heavy atoms inside the cylinder of that radius over the
    group's own axial extent.  This is the static steric part of a flip barrier, not the
    barrier: no electrostatics and no relaxation of the neighbours as the ring turns.

    Returns:
        dict: arrays R_sweep (nm) and blocking, one entry per ring, and the number of
            rings whose axis fell back to the farthest-apart pair
    """
    all_rings = set().union(*(set(r) for r in S.rings)) if S.rings else set()
    heavy = np.array(sorted(i for i in S.pos if S.heavy(i)))
    hx = np.array([S.pos[i] for i in heavy])
    R, N, fallback = [], [], 0
    for ring in rings:
        ia, ib, para = flip_axis(S, ring)
        if not para:
            fallback += 1
        a0 = S.pos[ia]
        a1 = a0 + S.mic(S.pos[ib] - a0)
        ax = (a1 - a0) / np.linalg.norm(a1 - a0)
        rs = set(ring)
        own = rs | {k for i in ring for k in S.adj[i] if k not in all_rings}
        own |= {k for i in own - rs for k in S.adj[i] if S.element.get(k) == 'H'}
        q = S.unwrap(sorted(own), anchor=a0)
        rad = np.linalg.norm(np.cross(q - a0, ax), axis=1)
        rsweep = rad.max()
        t_own = (q - a0) @ ax
        v = S.mic(hx - a0)
        t = v @ ax
        perp = np.linalg.norm(v - t[:, None] * ax, axis=1)
        inside = (perp < rsweep) & (t > t_own.min()) & (t < t_own.max())
        R.append(rsweep)
        N.append(int(sum(1 for k in heavy[inside] if k not in own)))
    return {'R_sweep': np.array(R), 'blocking': np.array(N), 'axis_fallback': fallback}


def _donor_groups(S, args):
    """(label, atom list) per donor selection: each ``-donors`` entry is a type or a
    comma-joined set of types, and ``-bridge`` adds the bridges' heavy atoms."""
    out = [(d, S.atoms_of_types(d.split(','))) for d in (args.donors or [])]
    if args.bridge:
        heavy = [i for g in bridge_groups(S, args.bridge) for i in g if S.heavy(i)]
        out.append((f'bridge({args.bridge}) heavy atoms', heavy))
    return out


def ring_contacts(args):
    """CLI handler for ``htpolynet ring-contacts``.

    Args:
        args (argparse.Namespace): parsed arguments; see ``htpolynet ring-contacts -h``
    """
    import os

    from ..core.topocoord import TopoCoord

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    for f in (args.top, args.gro):
        if not os.path.exists(f):
            logger.error(f'{f} not found')
            raise SystemExit(1)
    S = RingSystem(TopoCoord(topfilename=args.top, grofilename=args.gro), max_ring=args.max_ring)
    found = {}
    for c in S.composition:
        found[c] = found.get(c, 0) + 1
    print('rings: ' + (', '.join(f'{n} {c}' for c, n in sorted(found.items())) or 'none'))
    comps = args.rings or (['C6'] if args.analysis == 'flip' else ['C3N3'])
    missing = [c for c in comps if c not in found]
    if missing:
        logger.error(f'no ring of composition {", ".join(missing)} in this system')
        raise SystemExit(1)

    if args.analysis == 'flip':
        print(f'{"rings":>8}{"n":>7}{"R_sweep (nm)":>18}{"blocking heavy atoms":>24}')
        for c in comps:
            res = ring_flip_clearance(S, S.rings_of(c))
            R, N = res['R_sweep'], res['blocking']
            print(f'{c:>8}{len(R):7d}{R.mean():11.3f} +/- {R.std():.3f}'
                  f'{N.mean():15.2f} +/- {N.std() / np.sqrt(len(N)):.2f}')
            if res['axis_fallback']:
                logger.warning(f'{res["axis_fallback"]} {c} ring(s) had no unique substituted '
                               'para pair; their axis is the farthest-apart pair')
        return

    groups = _donor_groups(S, args)
    if not groups:
        logger.error('select atoms with -donors <type> [...] and/or -bridge <composition>')
        raise SystemExit(1)
    if args.analysis == 'contacts':
        region = (f'{args.rmin * 10:.1f}-{args.rmax * 10:.1f}' if args.rmin
                  else f'< {args.rmax * 10:.1f}')
        print(f'face-on: {region} A from the centroid, < {args.angle:g} deg from the normal; '
              + (f'rings within {args.exclude_bonds} bonds excluded' if args.exclude == 'bonds'
                 else "rings of the donor's own residue excluded"))
        print(f'{"donors":>28}{"n":>7}{"rings":>8}{"observed":>10}{"expected":>10}{"enhancement":>13}')
        for lab, don in groups:
            for c in comps:
                r = face_on_contacts(S, don, S.rings_of(c), rmax=args.rmax, max_angle=args.angle,
                                     exclude_bonds=args.exclude_bonds, rmin=args.rmin,
                                     exclude=args.exclude)
                print(f'{lab:>28}{len(don):7d}{c:>8}{r["observed"]:10d}{r["expected"]:10.1f}'
                      f'{r["enhancement"]:12.2f}x')
        return

    # rdf: one group per atom for -donors, the whole bridge group for -bridge
    sel = [(lab, [[i] for i in don]) for lab, don in groups if not lab.startswith('bridge(')]
    if args.bridge:
        sel.append((f'bridge({args.bridge})', bridge_groups(S, args.bridge)))
    cols, names = [], []
    for lab, grp in sel:
        for c in comps:
            x, g = centroid_rdf(S, grp, S.rings_of(c), rmax=args.rdf_max, dr=args.dr)
            w = (x > 0.3) & (x < 1.2)
            pk = x[w][np.argmax(g[w])] if w.any() else float('nan')
            print(f'{lab} -> {c}: {len(grp)} group(s); first peak in 3-12 A at '
                  f'{pk * 10:.2f} A, g = {g[w].max() if w.any() else float("nan"):.2f}')
            cols.append(g)
            names.append(f'{lab}->{c}'.replace(' ', '_'))
    if args.o:
        np.savetxt(args.o, np.column_stack([x] + cols), fmt='%.4f',
                   header='r_nm ' + ' '.join(names))
        logger.info(f'Wrote {args.o}')
