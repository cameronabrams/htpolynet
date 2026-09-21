"""Finds the triples of reactive groups that one ring-closing reaction can join.

The pairwise bond search asks, for each reaction bond, which A atom is near which B
atom.  A cyclotrimerization does not decompose that way: three cyanate groups have to
be mutually close, and choosing two of them well can still leave the third far off.
So candidates here are *triangles*, scored as a whole, and each group may join only
one of them.

Kept free of TopoCoord and MD so it can be exercised on coordinates alone.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

import numpy as np
import pandas as pd

from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)


def reactive_sites(adf, resname, groups=(('N1', 'C1'),), require_z=True):
    """Finds the unreacted reactive groups a ring-closing reaction can use.

    A group is one (donor, acceptor) pair of atoms within a residue: for a cyanate
    ester, the nitrogen and carbon of one -O-C#N, with the ring alternating them.  A
    residue may carry several -- a bisphenol dicyanate has two, which is what makes it
    a crosslinker -- so a group is identified by residue *and* atom names, never by
    residue alone.

    Args:
        adf (pandas.DataFrame): the system's atoms, with globalIdx, resNum, resName,
            atomName, molecule and z
        resname (str): residue name carrying the groups
        groups (iterable): (donor name, acceptor name) pairs, one per group a residue
            of this kind carries
        require_z (bool): keep only groups whose atoms both still have z > 0

    Returns:
        pandas.DataFrame: one row per group, with resNum, molecule, group (the donor
            and acceptor names), and the two atom indices
    """
    sel = adf[adf['resName'] == resname]
    rows = []
    for donor, acceptor in groups:
        d = sel[sel['atomName'] == donor].set_index('resNum')
        a = sel[sel['atomName'] == acceptor].set_index('resNum')
        common = d.index.intersection(a.index)
        if common.empty:
            continue
        part = pd.DataFrame({
            'resNum': common,
            'molecule': d.loc[common, 'molecule'].to_numpy(),
            'group': [f'{donor}-{acceptor}'] * len(common),
            'donor': d.loc[common, 'globalIdx'].to_numpy(),
            'acceptor': a.loc[common, 'globalIdx'].to_numpy(),
        })
        if require_z and 'z' in sel.columns:
            keep = (d.loc[common, 'z'].to_numpy() > 0) & (a.loc[common, 'z'].to_numpy() > 0)
            part = part[keep]
        rows.append(part)
    if not rows:
        return pd.DataFrame(columns=['resNum', 'molecule', 'group', 'donor', 'acceptor'])
    return pd.concat(rows, ignore_index=True).sort_values(['resNum', 'group']).reset_index(drop=True)


def _wrapped(positions, box):
    """Positions folded into [0, box), as a periodic KD-tree requires."""
    return np.mod(np.asarray(positions, dtype=float), np.asarray(box, dtype=float))


def _mic(d, box):
    """Minimum-image displacement."""
    box = np.asarray(box, dtype=float)
    return d - np.round(d / box) * box


def candidate_triples(sites, positions, radius, box, same_molecule=False, same_residue=False):
    """Scores every triangle of groups whose three new bonds would all be short enough.

    A triangle of groups admits two ring orientations -- which group donates to which
    -- and they are not equivalent, so both are scored and the better kept.

    Args:
        sites (pandas.DataFrame): from :func:`reactive_sites`
        positions (dict): global atom index -> position, in nm
        radius (float): longest new bond allowed, in nm
        box (array-like): box diagonal, in nm
        same_molecule (bool): allow two groups of one molecule in a ring; defaults to
            False, which is what the published method does and what the geometry of a
            para-para bisphenol makes impossible anyway
        same_residue (bool): allow two groups of one residue in a ring; defaults to
            False.  Both arms of one dicyanate closing one triazine is a 14-membered
            ring, which for these monomers the O...O span cannot reach, and the
            pairwise bondtest would hit an assertion on such a pair rather than
            decline it

    Returns:
        list: (score, (i, j, k), [(ai, aj), ...]) tuples, ascending by score, where the
            indices are rows of ``sites`` in ring order and the score is the total
            length of the three new bonds
    """
    if sites.shape[0] < 3:
        return []
    box = np.asarray(box, dtype=float)
    # a group's neighbourhood is judged on its acceptor, which is what the previous
    # group's donor has to reach
    acc = np.array([positions[int(x)] for x in sites['acceptor']])
    don = np.array([positions[int(x)] for x in sites['donor']])
    tree = cKDTree(_wrapped(acc, box), boxsize=box)
    # a ring needs donor-to-acceptor reach; 2*radius on acceptor separation is a
    # generous envelope that cannot exclude a feasible ring
    pairs = tree.query_pairs(2.0 * radius)
    neighbors = {}
    for i, j in pairs:
        neighbors.setdefault(i, set()).add(j)
        neighbors.setdefault(j, set()).add(i)
    resnum = sites['resNum'].to_numpy()
    molecule = sites['molecule'].to_numpy()

    def bond_length(a, b):
        return float(np.linalg.norm(_mic(don[a] - acc[b], box)))

    out = []
    seen = set()
    for i, j in pairs:
        for k in neighbors.get(i, set()) & neighbors.get(j, set()):
            tri = tuple(sorted((i, j, k)))
            if tri in seen:
                continue
            seen.add(tri)
            if not same_residue and len(set(resnum[list(tri)])) < 3:
                continue
            if not same_molecule and len(set(molecule[list(tri)])) < 3:
                continue
            best = None
            for order in ({(tri[0], tri[1], tri[2]), (tri[0], tri[2], tri[1])}):
                a, b, c = order
                bonds = [(a, b), (b, c), (c, a)]
                lengths = [bond_length(x, y) for x, y in bonds]
                if max(lengths) > radius:
                    continue
                score = float(sum(lengths))
                if best is None or score < best[0]:
                    best = (score, order, bonds)
            if best is not None:
                out.append(best)
    out.sort(key=lambda x: x[0])
    logger.debug(f'{len(out)} candidate triple(s) from {sites.shape[0]} group(s) within {radius} nm')
    return out


def select_disjoint(candidates, max_triples=0):
    """Chooses triples so that no group is used twice, shortest total bonds first.

    Greedy, which is what a cutoff-driven search does with pairs.  A better packing
    exists in principle -- the same combinatorial problem the single-step annealer
    solves -- and this is where that would go.

    Args:
        candidates (list): from :func:`candidate_triples`, ascending by score
        max_triples (int): stop after this many; 0 for no limit

    Returns:
        list: the accepted candidates, in acceptance order
    """
    used, chosen = set(), []
    for score, order, bonds in candidates:
        if used.intersection(order):
            continue
        used.update(order)
        chosen.append((score, order, bonds))
        if max_triples and len(chosen) >= max_triples:
            break
    return chosen


def triple_bonds_dataframe(chosen, sites, product, order=1):
    """Turns chosen triples into the bond table the topology update consumes.

    Args:
        chosen (list): from :func:`select_disjoint`
        sites (pandas.DataFrame): from :func:`reactive_sites`
        product (str): the reaction's product name, recorded per bond as elsewhere
        order (int): bond order

    Returns:
        pandas.DataFrame: three rows per triple, with ai, aj, ri, rj, order,
            reactantName and a triple column grouping them
    """
    rows = []
    for n, (score, ring, bonds) in enumerate(chosen):
        for a, b in bonds:
            rows.append({'ai': int(sites.at[a, 'donor']), 'aj': int(sites.at[b, 'acceptor']),
                         'ri': int(sites.at[a, 'resNum']), 'rj': int(sites.at[b, 'resNum']),
                         'order': order, 'reactantName': product, 'triple': n})
    return pd.DataFrame(rows, columns=['ai', 'aj', 'ri', 'rj', 'order', 'reactantName', 'triple'])


def residue_map(chosen, sites, template_resids):
    """Maps each triple's template residues onto its instance residues.

    ``productsplice.map_product_from_template`` needs this and will not guess it.

    Args:
        chosen (list): from :func:`select_disjoint`
        sites (pandas.DataFrame): from :func:`reactive_sites`
        template_resids (list): the template's residue numbers, in the order the
            reaction lists its reactants

    Returns:
        list: one {template resnr: instance resNum} dict per triple
    """
    out = []
    for score, ring, bonds in chosen:
        if len(template_resids) != len(ring):
            raise ValueError(f'template has {len(template_resids)} residue(s) for a '
                             f'{len(ring)}-group ring')
        out.append({int(t): int(sites.at[s, 'resNum']) for t, s in zip(template_resids, ring)})
    return out
