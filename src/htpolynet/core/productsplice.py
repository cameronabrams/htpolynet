"""Splices a whole product template onto the residues of one reaction event.

``TopoCoord.map_from_templates`` works one bond at a time: it identifies each new
bond by its local context, then copies the interactions that contain that bond from
whichever template matches.  That is right for a reaction that forms one bond between
two reactants, and it cannot form a ring.

A cyclotrimerization bonds three cyanate groups 1-2, 2-3 and 3-1.  In the trimer
template all three ring bonds exist; in the system, partway through the batch, only
some do.  So the bond-by-bond context of an instance bond never matches the
template's, in either direction, and the per-bond path cannot be made to work by
ordering the bonds differently.

This module takes the other route: given the participating residues, make the mapped
region of the system match the template exactly -- atom types, charges, and every
bonded interaction among those residues.  The caller says which residues correspond
to which template residues, so no context matching is needed at all.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

import pandas as pd

logger = logging.getLogger(__name__)

_SECTIONS_ = ('bonds', 'angles', 'dihedrals', 'pairs')
"""Interaction sections copied from the template, in the order they are applied."""

_ATOM_COLUMNS_ = {'bonds': ('ai', 'aj'), 'angles': ('ai', 'aj', 'ak'),
                  'dihedrals': ('ai', 'aj', 'ak', 'al'), 'pairs': ('ai', 'aj')}


def _key(section, row):
    """Returns a direction-independent key for one interaction row.

    Bonds and pairs are unordered; an angle is the same read either way, so its end
    atoms are sorted with the apex held; a dihedral is the same read either way.
    """
    idx = tuple(int(row[c]) for c in _ATOM_COLUMNS_[section])
    if section == 'angles':
        i, j, k = idx
        return (min(i, k), j, max(i, k))
    return min(idx, idx[::-1])


def _keys_of(section, df):
    """Returns the set of keys present in an interaction dataframe."""
    if df is None or df.empty:
        return set()
    return {_key(section, r) for _, r in df.iterrows()}


def atom_map(instance_atoms, template_atoms, resid_map):
    """Maps template atoms to instance atoms, residue by residue, on atom name.

    Args:
        instance_atoms (pandas.DataFrame): the system's [ atoms ], with nr, atom, resnr
        template_atoms (pandas.DataFrame): the template's [ atoms ], same columns
        resid_map (dict): {template resnr: instance resnr} for every participating residue

    Returns:
        tuple: (temp2inst, unmapped) -- the mapping, and the template atom numbers that
            have no instance counterpart
    """
    temp2inst = {}
    for t_res, i_res in resid_map.items():
        t = template_atoms[template_atoms['resnr'] == t_res][['nr', 'atom']]
        i = instance_atoms[instance_atoms['resnr'] == i_res][['nr', 'atom']]
        merged = t.merge(i, on='atom', how='left', suffixes=('_t', '_i'))
        for _, r in merged.iterrows():
            if pd.notna(r['nr_i']):
                temp2inst[int(r['nr_t'])] = int(r['nr_i'])
    mapped_t = set(temp2inst)
    unmapped = [int(x) for x in template_atoms[template_atoms['resnr'].isin(resid_map)]['nr']
                if int(x) not in mapped_t]
    return temp2inst, unmapped


def map_product_from_template(TC, template, resid_map, new_bonds=(), strict=True):
    """Makes the mapped residues of TC match a product template exactly.

    Copies atom types and charges for every mapped atom, and every bond, angle,
    dihedral and 1-4 pair of the template whose atoms are all mapped -- replacing the
    instance's own parameters where they exist and adding them where they do not.
    That is what a ring closure needs: the interactions it changes include ones
    wholly inside a residue (a cyanate C#N becoming a ring bond), which a per-bond
    splice never touches.

    The caller is responsible for creating the new bonds themselves (and for any
    sacrificial-atom deletion) before calling this, exactly as with
    ``map_from_templates``, and for settling the system charge afterwards -- the
    touched atoms are returned for that purpose.

    Args:
        TC (TopoCoord): the system
        template (Molecule): the product template
        resid_map (dict): {template resnr: instance resnr} for every participating residue
        new_bonds (iterable): (ai, aj) instance pairs formed for this event, checked against
            the template
        strict (bool): raise when the template and instance disagree about which bonds
            exist among the mapped residues; defaults to True

    Returns:
        dict: ``atoms`` (instance atom numbers touched), ``types`` and ``charges``
            (how many changed), ``added`` and ``replaced`` (per section), and
            ``unmapped`` (template atoms with no counterpart)

    Raises:
        ValueError: if a new bond has no counterpart in the template, or, under
            strict, if a template bond among mapped residues is missing from the
            instance
    """
    inst_top, temp_top = TC.Topology, template.TopoCoord.Topology
    temp2inst, unmapped = atom_map(inst_top.D['atoms'], temp_top.D['atoms'], resid_map)
    if unmapped:
        logger.debug(f'{len(unmapped)} template atom(s) have no instance counterpart and are skipped')

    # Every bond the caller formed must be one the template has, or the residue map is
    # wrong and everything below would be quietly misapplied.
    inst2temp = {v: k for k, v in temp2inst.items()}
    temp_bond_keys = _keys_of('bonds', temp_top.D.get('bonds'))
    for ai, aj in new_bonds:
        t = (inst2temp.get(int(ai)), inst2temp.get(int(aj)))
        if None in t or min(t, t[::-1]) not in temp_bond_keys:
            raise ValueError(f'new bond {ai}-{aj} maps to {t}, which is not a bond of template '
                             f'{template.name}; the residue map is wrong')

    stats = {'atoms': sorted(temp2inst.values()), 'types': 0, 'charges': 0,
             'added': {}, 'replaced': {}, 'unmapped': unmapped}

    # types and charges, atom by atom
    iat, tat = inst_top.D['atoms'], temp_top.D['atoms']
    t_by_nr = tat.set_index('nr')
    for t_nr, i_nr in temp2inst.items():
        t_row = t_by_nr.loc[t_nr]
        mask = iat['nr'] == i_nr
        if iat.loc[mask, 'type'].iloc[0] != t_row['type']:
            iat.loc[mask, 'type'] = t_row['type']
            stats['types'] += 1
        if iat.loc[mask, 'charge'].iloc[0] != t_row['charge']:
            iat.loc[mask, 'charge'] = float(t_row['charge'])
            stats['charges'] += 1

    for section in _SECTIONS_:
        t_df = temp_top.D.get(section)
        if t_df is None or t_df.empty:
            continue
        cols = _ATOM_COLUMNS_[section]
        keep = t_df[[all(int(r[c]) in temp2inst for c in cols) for _, r in t_df.iterrows()]]
        if keep.empty:
            continue
        mapped = keep.copy()
        for c in cols:
            mapped[c] = [temp2inst[int(x)] for x in keep[c]]
        i_df = inst_top.D.get(section)
        existing = _keys_of(section, i_df)
        new_keys = {_key(section, r) for _, r in mapped.iterrows()}
        if section == 'bonds':
            missing = new_keys - existing
            if missing and strict:
                raise ValueError(f'template {template.name} has {len(missing)} bond(s) among the mapped '
                                 f'residues that the instance does not: {sorted(missing)[:4]}. Form the '
                                 f'reaction bonds before splicing, or fix the residue map.')
        # replace the instance's own version of every interaction the template covers
        if i_df is not None and not i_df.empty:
            drop = [_key(section, r) in new_keys for _, r in i_df.iterrows()]
            stats['replaced'][section] = int(sum(drop)) if any(drop) else 0
            inst_top.D[section] = pd.concat([i_df[[not d for d in drop]], mapped], ignore_index=True)
        else:
            stats['replaced'][section] = 0
            inst_top.D[section] = mapped.reset_index(drop=True)
        stats['added'][section] = len(new_keys) - stats['replaced'][section]

    # pandas upcasts an int column to float wherever a concat had to align frames, and
    # the topology writer and delete_atoms both assert int atom indices
    for section in _SECTIONS_:
        df = inst_top.D.get(section)
        if df is None or df.empty:
            continue
        for c in _ATOM_COLUMNS_[section] + ('funct',):
            if c in df.columns:
                df[c] = df[c].astype(int)
    if 'bonds' in inst_top.D:
        from ..geometry.bondlist import Bondlist
        inst_top.bondlist = Bondlist.fromDataFrame(inst_top.D['bonds'])
    logger.debug(f'product splice from {template.name}: {stats["types"]} type(s) and '
                 f'{stats["charges"]} charge(s) changed, added {stats["added"]}, replaced {stats["replaced"]}')
    return stats
