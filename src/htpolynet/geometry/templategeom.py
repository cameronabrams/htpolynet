"""Makes a freshly assembled product template a plausible molecule before it is parameterized.

A template built by bonding reactants together inherits its geometry from how those
reactants were placed, not from the chemistry of the product.  For a condensation that
is close enough: one new bond, formed between atoms that were already adjacent, and the
minimization after parameterization tidies the rest.

A ring closure is not close enough.  Three cyanate groups pulled toward each other close
to about 0.22 nm and no further -- the restraint ladder is fighting each group's own
C#N, whose sp parameters hold it linear, and no force constant wins that.  The bonds are
recorded, so the connectivity says triazine, but the coordinates still say three
cyanates 2.2 A apart.

That matters because charges are computed from the geometry.  antechamber repairs the
connectivity it is given -- it perceives the ring as aromatic and assigns `ca`/`nb`
correctly -- but AM1-BCC takes Mulliken charges from a semiempirical calculation on the
coordinates, and on those coordinates the calculation sees three separate cyanates.
Measured on a methyl cyanate trimer: each ring carbon comes out at +0.699 against
+0.909 for the same molecule properly embedded, 0.21 e low on every ring carbon, with
the deficit smeared over the rest of the molecule so the net charge still checks out.
htpolynet-study measured 0.2113 e on a bisphenol A dicyanate trimer independently.

Correcting the bond orders alone does not fix it -- that arm of the experiment was
bit-identical to the unfixed one, because the QM step does not read bond orders.  The
geometry has to be right before antechamber sees it, which is what this module does.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)

_ORDER_ = {'1': 1.0, '2': 2.0, '3': 3.0, 'ar': 1.5, 'am': 1.0}
"""mol2 bond order strings as numbers."""

_MAX_PLAUSIBLE_ = {1.0: 0.20, 1.5: 0.18, 2.0: 0.18, 3.0: 0.16}
"""Longest believable length, in nm, for a bond of each order between light elements.

Generous: the longest ordinary single bond here is C-S at about 0.182 nm, and the point
is to catch a bond that is not a bond at all, not to police a strained one.
"""


def _has_mol2_bonds(TC):
    """True when a TopoCoord carries the bond orders these checks read.

    A monomer reaches parameterization with its structure still in a file rather than
    loaded, and a system-scale TopoCoord never carries mol2 bonds at all, so both simply
    have nothing to check here.
    """
    top = getattr(TC, 'Topology', None)
    D = getattr(top, 'D', None) if top is not None else None
    if not D or 'atoms' not in D:
        return False
    mb = D.get('mol2_bonds')
    if mb is None or mb.empty:
        return False
    A = getattr(getattr(TC, 'Coordinates', None), 'A', None)
    return A is not None and not A.empty


def _rdkit_mol(TC):
    """Builds an RDKit molecule from a TopoCoord's elements, mol2 bonds and coordinates.

    Returns:
        tuple: (mol, index_map) where index_map takes an htpolynet atom number to an
            RDKit atom index, or (None, None) if the molecule cannot be sanitized
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem  # noqa: F401  (registers conformer machinery)
    from rdkit.Geometry import Point3D

    from ..core.topology import _element_of

    if not _has_mol2_bonds(TC):
        return None, None
    at = TC.Topology.D['atoms']
    mb = TC.Topology.D['mol2_bonds']
    A = TC.Coordinates.A.set_index('globalIdx')

    bt = {1.0: Chem.BondType.SINGLE, 1.5: Chem.BondType.AROMATIC,
          2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE}
    rw = Chem.RWMol()
    index_map = {}
    for r in at.itertuples():
        e = _element_of(r.atom, r.type)
        e = e.capitalize() if len(e) == 2 else e
        index_map[int(r.nr)] = rw.AddAtom(Chem.Atom(e))
    for r in mb.itertuples():
        rw.AddBond(index_map[int(r.ai)], index_map[int(r.aj)],
                   bt[_ORDER_.get(str(r.order), 1.0)])
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        logger.debug(f'RDKit could not sanitize the template: {e}')
        return None, None

    conf = Chem.Conformer(mol.GetNumAtoms())
    for nr, idx in index_map.items():
        p = A.loc[nr]
        # htpolynet holds coordinates in nm; RDKit force fields want Angstrom
        conf.SetAtomPosition(idx, Point3D(float(p.posX) * 10.0, float(p.posY) * 10.0,
                                          float(p.posZ) * 10.0))
    mol.AddConformer(conf, assignId=True)
    return mol, index_map


def _minimize(mol, maxIters):
    """Runs MMFF if it has parameters for this molecule, UFF otherwise.

    Returns:
        tuple: (succeeded, energy before, energy after)
    """
    from rdkit.Chem import AllChem
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        ff = (AllChem.MMFFGetMoleculeForceField(mol, props) if props is not None
              else AllChem.UFFGetMoleculeForceField(mol))
        if ff is None:
            return False, 0.0, 0.0
        before = ff.CalcEnergy()
        ff.Minimize(maxIts=maxIters)
        return True, before, ff.CalcEnergy()
    except Exception as e:
        logger.debug(f'force-field minimization failed: {e}')
        return False, 0.0, 0.0


def relax_geometry(TC, name='', maxIters=5000):
    """Relaxes a template's coordinates against its own connectivity, in place.

    Tries first from the coordinates the template already has, so a product keeps the
    arrangement its reactants were placed in and only the bonds that were never closed
    actually move.  If that minimization cannot get going -- a bond recorded across 2.2 A
    can leave the starting geometry strained enough that the line search fails on its
    first step -- it falls back to embedding the molecule from scratch.  A fresh embed
    discards the reactants' arrangement, which is a real loss but a much smaller one than
    parameterizing an open ring: the template's geometry is a starting point that gets
    minimized again under GAFF, while its charges are computed once and kept.

    Args:
        TC (TopoCoord): the template
        name (str): the molecule's name, for logging
        maxIters (int): force-field iteration cap

    Returns:
        bool: True if the coordinates were updated
    """
    from rdkit.Chem import AllChem

    mol, index_map = _rdkit_mol(TC)
    if mol is None:
        logger.warning(f'{name}: template geometry left as built; RDKit could not read it. '
                       f'If this product closes a ring, its charges would be computed on an '
                       f'open ring.')
        return False

    ok, before, after = _minimize(mol, maxIters)
    how = 'relaxed against its own bonds'
    if not ok:
        logger.debug(f'{name}: minimizing from the built geometry failed; re-embedding')
        mol.RemoveAllConformers()
        if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
            logger.warning(f'{name}: template geometry left as built; RDKit could neither '
                           f'minimize nor embed it.')
            return False
        ok, before, after = _minimize(mol, maxIters)
        how = 're-embedded and relaxed'
        if not ok:
            logger.warning(f'{name}: template geometry left as built; no force field would '
                           f'run on it.')
            return False

    conf = mol.GetConformer()
    A = TC.Coordinates.A
    pos = {nr: conf.GetAtomPosition(idx) for nr, idx in index_map.items()}
    for col, comp in (('posX', 'x'), ('posY', 'y'), ('posZ', 'z')):
        A[col] = [getattr(pos[int(g)], comp) / 10.0 if int(g) in pos else v
                  for g, v in zip(A['globalIdx'], A[col])]
    logger.info(f'{name}: template geometry {how} ({before:.1f} -> {after:.1f} kcal/mol)')
    return True


def geometry_complaints(TC):
    """Finds what would make a template a bad thing to hand to a charge calculation.

    Two checks, because they fail independently: an atom whose recorded bond orders
    exceed what its element can carry, and a bond whose length is not credible for its
    order.  The second is what catches a ring that never closed -- an aromatic C-N
    recorded at 2.2 A is unmistakable -- and it catches it whether or not the orders
    were ever repaired.

    Args:
        TC (TopoCoord): the template

    Returns:
        list: human-readable complaints, empty when the template is fit to parameterize
    """
    from ..core.topology import _element_of

    out = []
    if not _has_mol2_bonds(TC):
        return out
    mb = TC.Topology.D['mol2_bonds']
    at = TC.Topology.D['atoms']
    element = {int(r.nr): _element_of(r.atom, r.type) for r in at.itertuples()}
    name = {int(r.nr): r.atom for r in at.itertuples()}
    limit = {'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'CL': 1, 'BR': 1, 'I': 1}

    val = {}
    for r in mb.itertuples():
        o = _ORDER_.get(str(r.order), 1.0)
        val[int(r.ai)] = val.get(int(r.ai), 0.0) + o
        val[int(r.aj)] = val.get(int(r.aj), 0.0) + o
    for a, v in sorted(val.items()):
        e = element.get(a)
        if e in limit and v > limit[e] + 0.01:
            out.append(f'atom {a} ({name[a]}, {e}) has bond orders summing to {v:g}, '
                       f'more than {e} can carry')

    A = TC.Coordinates.A.set_index('globalIdx')
    for r in mb.itertuples():
        i, j = int(r.ai), int(r.aj)
        if i not in A.index or j not in A.index:
            continue
        d = float(np.linalg.norm(A.loc[i, ['posX', 'posY', 'posZ']].to_numpy(dtype=float)
                                 - A.loc[j, ['posX', 'posY', 'posZ']].to_numpy(dtype=float)))
        o = _ORDER_.get(str(r.order), 1.0)
        if d > _MAX_PLAUSIBLE_.get(o, 0.20):
            out.append(f'bond {name[i]}({i})-{name[j]}({j}) is recorded at order '
                       f'{r.order} but is {d * 10:.2f} A long')
    return out


def check_ring_closed(TC, name):
    """Refuses to parameterize a ring-closing product whose ring is still open.

    Scope matters here.  EVERY htpolynet template reaches antechamber with its new bonds
    about 2.2 A long -- a styrene dimer as much as a triazine -- because the bond is
    recorded before anything relaxes it.  For a single new bond that has always been
    fine in practice, and this check deliberately says nothing about it: changing how
    those are parameterized would move the charges of every template htpolynet has ever
    built, which is not a thing to do blind.  See ROADMAP.md.

    What is not fine is a ring, because its charges were measured wrong: 0.21 e low on
    every ring carbon, twice, on two different monomers.  So after
    :func:`relax_geometry` has had its turn, a ring that is still open means the repair
    did not work and the charges would be silently wrong again.

    Note that antechamber echoes its input coordinates into its output, so this must be
    checked on what we hand it.  Its output geometry is not evidence of anything.

    Raises:
        RuntimeError: listing every complaint, because they usually share one cause
    """
    complaints = geometry_complaints(TC)
    if not complaints:
        return
    raise RuntimeError(
        f'{name} closes a ring, but after relaxation it is still not a plausible '
        f'molecule, so its charges would be wrong and it will not be parameterized:\n  '
        + '\n  '.join(complaints) +
        f'\nCharges are computed from the geometry, not from the bond orders, and a ring '
        f'left open is charged as the separate groups it still looks like.')
