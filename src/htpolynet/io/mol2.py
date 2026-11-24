"""Sybyl MOL2 file reader and writer.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

from io import StringIO

import numpy as np
import pandas as pd

from ..geometry.bondlist import Bondlist

logger = logging.getLogger(__name__)

MOL2_ATOM_ATTRIBUTES = [
    'globalIdx', 'atomName', 'posX', 'posY', 'posZ',
    'type', 'resNum', 'resName', 'charge',
]
MOL2_BOND_ATTRIBUTES = ['bondIdx', 'ai', 'aj', 'order']
MOL2_BOND_TYPES = {k: v for k, v in zip(MOL2_BOND_ATTRIBUTES, [int, int, int, str])}


def read(filename):
    """Parse a Sybyl MOL2 file.

    Only MOLECULE, ATOM, and BOND sections are read.
    All lengths are converted from Ångström to nm.

    Args:
        filename (str): path to the .mol2 file

    Returns:
        dict: parsed data with keys:
            - 'name'         (str)
            - 'N'            (int): atom count
            - 'metadat'      (dict): nBonds, nSubs, nFeatures, nSets, mol2type, mol2chargetype
            - 'A'            (pd.DataFrame): atom data frame (positions in nm)
            - 'mol2_bonds'   (pd.DataFrame)
            - 'mol2_bondlist' (Bondlist)
    """
    with open(filename, 'r') as f:
        rawsections = f.read().split('@<TRIPOS>')[1:]

    sections = {}
    for rs in rawsections:
        s = rs.split('\n')
        key = s[0].strip().lower()
        val = [a.strip() for a in s[1:] if len(a) > 0]
        if key == 'atom' or key == 'bond':
            val = StringIO('\n'.join(val))
        sections[key] = val

    name = sections['molecule'][0]
    imetadat = list(map(int, sections['molecule'][1].strip().split()))
    metadat = {
        'N':            imetadat[0],
        'nBonds':       imetadat[1],
        'nSubs':        imetadat[2],
        'nFeatures':    imetadat[3],
        'nSets':        imetadat[4],
        'mol2type':     sections['molecule'][2],
        'mol2chargetype': sections['molecule'][3],
    }
    N = metadat['N']

    A = pd.read_csv(sections['atom'], sep=r'\s+', names=MOL2_ATOM_ATTRIBUTES)
    A[['posX', 'posY', 'posZ']] *= [0.1, 0.1, 0.1]   # Å → nm
    N = A.shape[0]

    mol2_bonds = pd.read_csv(
        sections['bond'], sep=r'\s+',
        names=MOL2_BOND_ATTRIBUTES, dtype=MOL2_BOND_TYPES,
    )
    # Ensure ai < aj in every bond record
    for i, r in mol2_bonds.iterrows():
        ai, aj = r['ai'], r['aj']
        if aj < ai:
            logger.debug(f'mol2 bonds swapping {ai} and {aj}')
            mol2_bonds.iloc[i, mol2_bonds.columns == 'ai'] = aj
            mol2_bonds.iloc[i, mol2_bonds.columns == 'aj'] = ai

    mol2_bondlist = Bondlist.fromDataFrame(mol2_bonds)

    return {
        'name': name,
        'N': N,
        'metadat': metadat,
        'A': A,
        'mol2_bonds': mol2_bonds,
        'mol2_bondlist': mol2_bondlist,
    }


def write(filename, A, mol2_bonds, metadat, name,
          bondsDF=None, molname='', other_attributes=None):
    """Write a Sybyl MOL2 file.

    All coordinates are written in Ångström (positions in A are assumed to be in nm
    and are multiplied by 10 after centering).

    Args:
        filename (str): path to write
        A (pd.DataFrame): atom data frame
        mol2_bonds (pd.DataFrame): internal mol2 bonds (may be empty)
        metadat (dict): molecule metadata (nFeatures, nSets, mol2type, mol2chargetype …)
        name (str): molecule name used when molname is empty
        bondsDF (pd.DataFrame or None): caller-supplied bonds, overrides mol2_bonds when present
        molname (str): explicit molecule name for the MOLECULE section, defaults to ''
        other_attributes (pd.DataFrame or None): additional per-atom columns to merge into A
    """
    if bondsDF is None:
        bondsDF = pd.DataFrame()
    if other_attributes is None:
        other_attributes = pd.DataFrame()

    acopy = A.copy()

    # Determine which bond table to use
    if bondsDF.empty and (mol2_bonds is None or mol2_bonds.empty):
        logger.warning(f'Cannot write any bonds to MOL2 file {filename}')
        bdf = pd.DataFrame()
    elif (mol2_bonds is not None and not mol2_bonds.empty) and bondsDF.empty:
        bdf = mol2_bonds
    elif (not bondsDF.empty) and (mol2_bonds is None or mol2_bonds.empty):
        bdf = bondsDF
    else:
        logger.info('write_mol2 provided with both a bondsDF parameter and a mol2_bonds attribute')
        logger.info('Using the parameter')
        bdf = bondsDF

    for col in other_attributes.columns:
        acopy[col] = other_attributes[col]

    # Geometric centre for centering before writing
    com = np.array([acopy.posX.mean(), acopy.posY.mean(), acopy.posZ.mean()])

    atomformatters = [
        lambda x: f'{x:>7d}',
        lambda x: f'{x:<8s}',
        lambda x: f'{x:>9.4f}',
        lambda x: f'{x:>9.4f}',
        lambda x: f'{x:>9.4f}',
        lambda x: f'{x:<5s}',
        lambda x: f'{x:>3d}',
        lambda x: f' {x:<7s}',
        lambda x: f'{x:>9.4f}',
    ]
    bondformatters = [
        lambda x: f'{x:>6d}',
        lambda x: f'{x:>5d}',
        lambda x: f'{x:>5d}',
        lambda x: f'{str(x):>4s}',
    ]
    substructureformatters = [
        lambda x: f'{x:>6d}',
        lambda x: f'{x:<7s}',
        lambda x: f'{x:>6d}',
        lambda x: f'{x:<7s}',
    ]

    rdf = acopy[['resNum', 'resName']].copy().drop_duplicates()
    rdf['rootatom'] = [1] * len(rdf)
    rdf['residue'] = ['RESIDUE'] * len(rdf)

    N = acopy.shape[0]
    nBonds = bdf.shape[0]
    nSubs = len(rdf)

    with open(filename, 'w') as f:
        f.write('@<TRIPOS>MOLECULE\n')
        f.write(f'{molname if molname else name}\n')
        # each field is a space plus the old width, so the line is byte-identical
        # whenever a value fits; the old fixed widths ran 100+ residues into the
        # bond count ('  3600  3600200') and the header could not be read back
        f.write(' {:>5d} {:>5d} {:>2d} {:>2d} {:>2d}\n'.format(
            N, nBonds, nSubs,
            metadat.get('nFeatures', 0),
            metadat.get('nSets', 0),
        ))
        f.write(f"{metadat.get('mol2type', 'SMALL')}\n")
        f.write(f"{metadat.get('mol2chargetype', 'GASTEIGER')}\n")
        f.write('\n')
        f.write('@<TRIPOS>ATOM\n')
        # Convert nm → Å and centre
        pos = (acopy.loc[:, ['posX', 'posY', 'posZ']] - com) * 10.0
        acopy.loc[:, ['posX', 'posY', 'posZ']] = pos
        f.write(acopy.to_string(
            columns=MOL2_ATOM_ATTRIBUTES, header=False, index=False,
            formatters=atomformatters,
        ))
        f.write('\n')
        f.write('@<TRIPOS>BOND\n')
        if not bondsDF.empty:
            logger.debug('Mol2 bonds from outside')
            out_bdf = bondsDF[['bondIdx', 'ai', 'aj', 'order']].copy()
            out_bdf['bondIdx'] = out_bdf['bondIdx'].astype(int)
            out_bdf['ai'] = out_bdf['ai'].astype(int)
            out_bdf['aj'] = out_bdf['aj'].astype(int)
            f.write(out_bdf.to_string(
                columns=MOL2_BOND_ATTRIBUTES, header=False, index=False,
                formatters=bondformatters,
            ))
        elif mol2_bonds is not None and not mol2_bonds.empty:
            logger.debug(f'write_mol2 ({filename}): Mol2 bonds from mol2_bonds attribute')
            f.write(mol2_bonds.to_string(
                columns=MOL2_BOND_ATTRIBUTES, header=False, index=False,
                formatters=bondformatters,
            ))
        f.write('\n')
        f.write('@<TRIPOS>SUBSTRUCTURE\n')
        f.write(rdf.to_string(header=False, index=False, formatters=substructureformatters))
        f.write('\n')
