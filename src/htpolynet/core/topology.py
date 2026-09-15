"""Class for managing gromacs .top file data.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import json
import logging
import os

from copy import deepcopy
from itertools import product

import networkx as nx
import numpy as np
import pandas as pd

from networkx.readwrite import json_graph
from scipy.constants import physical_constants

from ..geometry.bondlist import Bondlist
from ..geometry.ring import Ring, RingList

logger = logging.getLogger(__name__)

_PAD_ = -99.99

def typeorder(a):
    """Correctly orders the tuple of atom types for particular interaction types.

    Args:
        a (tuple): atom indices/types from a [ bond ], [ pair ], [ angle ], or [ dihedral ] record

    Returns:
        tuple: same atom indices/types correctly ordered to allow for easy searching/sorting
    """
    assert type(a) == tuple, 'error: typeorder() requires a tuple argument'
    if len(a) == 2: # bond
        return a if a[0] < a[1] else a[::-1]
    elif len(a) == 3: # angle
        return a if a[0] < a[2] else a[::-1]
    elif len(a) == 4: # dihedral
        if a[1] == a[2]:
            return a if a[0] < a[3] else a[::-1]
        return a if a[1] < a[2] else a[::-1]
idxorder = typeorder  # same syntax to order global atom indices in an interaction index

def repeat_check(t, msg=''):
    """Checks for repeated index tuples.

    Args:
        t (list): list of index tuples
        msg (str): optional message, defaults to ''
    """
    for i in range(len(t)):
        for j in range(i + 1, len(t)):
            assert t[i] != t[j], f'Error: repeated index in {len(t)}-tuple {t}: t({i})={t[i]}\n{msg}'

def df_typeorder(df, typs):
    """Type-orders the atom type attributes in each row of dataframe df.

    Args:
        df (pandas.DataFrame): a Topology type-directive dataframe; [ atomtypes ], [ bondtypes ], etc.
        typs (list): list of type-attribute names; typically ['i','j',...]
    """
    for i in df.index:
        df.loc[i,typs] = typeorder(tuple(df.loc[i, typs]))

_GromacsIntegers_ = ('nr', 'atnum', 'resnr', 'ai', 'aj', 'ak', 'al', '#mols', 'nrexcl', 'funct', 'func', 'nbfunc', 'comb-rule')
_GromacsFloats_ = ('charge', 'mass', 'chargeB', 'massB', *tuple([f'c{i}' for i in range(5)]),
                 'b0', 'kb', 'th0', 'cth', 'rub', 'kub', 'phase', 'kd', 'pn',
                 'fudgeLJ', 'fudgeQQ')
def typedata(h, s):
    if h in _GromacsIntegers_:
        return int(s)
    if h in _GromacsFloats_:
        return float(s)
    return s

_GromacsExtensiveDirectives_ = ['atoms', 'pairs', 'bonds', 'angles', 'dihedrals']
_NonGromacsExtensiveDirectives_ = ['mol2_bonds']
_GromacsTopologyDirectiveOrder_ = ['defaults', 'atomtypes', 'bondtypes', 'angletypes', 'dihedraltypes', 'moleculetype', 'atoms', 'pairs', 'bonds', 'angles', 'dihedrals', 'system', 'molecules']
_GromacsTopologyDirectiveHeaders_ = {
    'atoms': ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass','typeB', 'chargeB', 'massB'],
    'pairs': ['ai', 'aj', 'funct', 'c0', 'c1'],
    'bonds': ['ai', 'aj', 'funct', 'c0', 'c1'],
    'angles': ['ai', 'aj', 'ak', 'funct', 'c0', 'c1'],
    'dihedrals': ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'],
    'atomtypes': ['name', 'atnum', 'mass', 'charge', 'ptype', 'sigma', 'epsilon'],
    'moleculetype': ['name', 'nrexcl'],
    'bondtypes': ['i', 'j', 'func', 'b0','kb'],
    'angletypes': ['i', 'j', 'k', 'func', 'th0', 'cth', 'rub', 'kub'],
    'dihedraltypes': ['i', 'j', 'k', 'l', 'func', 'phase', 'kd', 'pn'],
    'system': ['Name'],
    'molecules': ['Compound', '#mols'],
    'defaults': ['nbfunc', 'comb-rule', 'gen-pairs', 'fudgeLJ', 'fudgeQQ']
}

_GromacsTopologyIntegerColumns_ = { # columns that must be written without a trailing '.0'; the dataframes can carry these as float when pandas demotes during NaN handling, so we re-cast to nullable Int64 just before serialization
    'atoms': ['nr', 'resnr', 'cgnr'],
    'pairs': ['ai', 'aj', 'funct'],
    'bonds': ['ai', 'aj', 'funct'],
    'angles': ['ai', 'aj', 'ak', 'funct'],
    'dihedrals': ['ai', 'aj', 'ak', 'al', 'funct'],
    'atomtypes': ['atnum'],
    'moleculetype': ['nrexcl'],
    'bondtypes': ['func'],
    'angletypes': ['func'],
    'dihedraltypes': ['func', 'pn'],
    'molecules': ['#mols'],
    'defaults': ['nbfunc', 'comb-rule'],
}

_GromacsTopologyHashables_ = { # attributes/columns that should always have values, no NaNs; these are how each item is sorted
    'atoms': ['nr'],
    'pairs': ['ai', 'aj'],
    'bonds': ['ai', 'aj'],
    'angles': ['ai', 'aj', 'ak'],
    'dihedrals': ['ai', 'aj', 'ak', 'al'],
    'atomtypes': ['name'],
    'bondtypes': ['i', 'j'],
    'angletypes': ['i', 'j', 'k'],
    'dihedraltypes': ['i', 'j', 'k', 'l']
}

# htpolynet writes an entire build as ONE moleculetype, on purpose: a cured
# network is a single covalently connected molecule, and leftover monomers
# ride along in the same block.  The name used to default to the string
# 'None', which reads to a user like an unset field or a bug -- and users who
# hit GROMACS's "inconsistent shifts" warning under trjconv -pbc whole go
# looking at exactly this line.  The name must match between [ moleculetype ]
# and [ molecules ], so both defaults change together.
_WHOLE_SYSTEM_MOLECULETYPE_ = 'whole_system'

_GromacsTopologyDirectiveDefaults_ = {
    'system': ['A_generic_system'],
    'molecules': [_WHOLE_SYSTEM_MOLECULETYPE_, 1],
    'moleculetype': [_WHOLE_SYSTEM_MOLECULETYPE_, 3],
    'defaults': [1, 2, 'yes', 0.5, 0.83333333]
}

def conflicting_types(topologies, typenames=('bondtypes', 'angletypes', 'dihedraltypes')):
    """Finds bonded types that two groups of topologies parameterize differently.

    A merged topology keeps one table per type, shared by every molecule in
    it, so a type with different parameters in two molecules cannot be
    honored for both: drop_duplicates keeps both sets, and GROMACS applies
    them together.  A user frcmod loaded for some molecules and not others is
    the usual way to create one.

    Args:
        topologies (dict): {name: (Topology, group)}; types are compared only between topologies
            whose groups differ, e.g. a digest of the frcmod files each was built with
        typenames (tuple): type sections to compare

    Returns:
        list: (typename, type key, name of one topology, name of another) tuples, one per conflicting
            type; the key is the atom types followed by the function number
    """
    seen = {}
    conflicts = []
    for name, (T, group) in topologies.items():
        for typename in typenames:
            df = T.D.get(typename)
            if df is None or df.empty:
                continue
            hashables = _GromacsTopologyHashables_[typename]
            params = [c for c in df.columns if c not in hashables]
            rows = {}
            for r in df.itertuples(index=False):
                r = r._asdict()
                t = tuple(r[h] for h in hashables)
                # a proper and an improper can share atom types; they are different types
                key = (*min(t, t[::-1]), int(r['func']) if 'func' in r else 0)
                rows.setdefault(key, []).append(tuple(round(float(r[p]), 4) if isinstance(r[p], (int, float)) else r[p]
                                                      for p in params))
            for key, values in rows.items():
                values = tuple(sorted(values))
                first = seen.setdefault((typename, key), (name, group, values))
                if first[1] != group and first[2] != values:
                    conflicts.append((typename, key, first[0], name))
    return conflicts


def select_topology_type_option(options, typename='dihedraltypes', rule='stiffest'):
    """Selects from a list of topological interaction options of type typename using the provided rule.

    Args:
        options (list): list of parameterization options for a particular interaction
        typename (str): string designation of interaction type, defaults to 'dihedraltypes'
        rule (str): string describing the selection rule, defaults to 'stiffest'

    Returns:
        dict: the selected parameterization option (an element of options)
    """
    hashables = _GromacsTopologyHashables_[typename]
    headers = _GromacsTopologyDirectiveHeaders_[typename].copy()
    for i in hashables:
        headers.remove(i)
    if rule == 'stiffest' or rule == 'softest':
        if typename == 'bondtypes':
            parindex = headers.index('kb')
        elif typename == 'angletypes':
            parindex = headers.index('cth')
        elif typename == 'dihedraltypes':
            parindex = headers.index('kd')
    sorted_options = list(sorted(options, key=lambda x: x[parindex]))
    if rule == 'stiffest':
        return sorted_options[-1]
    elif rule == 'softest':
        return sorted_options[0]
    return []

class Topology:
    """ Class for handling gromacs topology data
    """
    def __init__(self, system_name=''):
        """Constructor for Topology class.

        Args:
            system_name (str): optional name of system, defaults to ''
        """
        ''' D: a dictionay keyed on Gromacs topology directives with values that are lists of
               one or more pandas dataframes corresponding to sections '''
        self.D: dict[str, pd.DataFrame] = {}
        for k, v in _GromacsTopologyDirectiveDefaults_.items():
            hdr = _GromacsTopologyDirectiveHeaders_[k]
            dfdict = {kk: [a] for kk, a in zip(hdr, v)}
            self.D[k] = pd.DataFrame(dfdict)
        self.D['system'] = pd.DataFrame({'name': [system_name]})
        ''' bondlist: a class that owns a dictionary keyed on atom global index with values that are lists of global atom indices bound to the key '''
        self.bondlist = Bondlist()
        self.residue_network = nx.Graph()
        self.rings = RingList([])
        self.empty = True

    @classmethod
    def read_top(cls, filename, pad=_PAD_):
        """Reads a GROMACS-style topology file and returns a Topology instance.

        Each directive section is stored as a pandas DataFrame. Duplicate 'dihedrals' and
        'dihedraltypes' sections are merged.

        Args:
            filename (str): name of GROMACS top file to read
            pad (float): padding value for missing fields, defaults to _PAD_

        Raises:
            KeyError: if an unrecognized topology directive is encountered

        Returns:
            Topology: a new Topology instance
        """
        assert os.path.exists(filename), f'Error: {filename} not found.'
        inst = cls()
        inst.filename = filename
        dirname = os.path.dirname(filename)
        inst.includes = []
        with open(filename, 'r') as f:
            data = f.read().split('[')
            stanzas = [('['+x).split('\n') for x in data][1:]
            for s in stanzas:
                directive = s[0].split()[1].strip()
                contentlines = [l for l in s[1:] if not (len(l.strip()) == 0 or l.strip().startswith(';'))]
                if not directive in _GromacsTopologyDirectiveHeaders_:
                    raise KeyError(f'unrecognized topology directive "{directive}"')
                header = _GromacsTopologyDirectiveHeaders_[directive]
                series = {k: [] for k in header}
                for line in contentlines:
                    if line.startswith('#include'):
                        tokens = [x.strip() for x in line.split()]
                        inst.includes.append(tokens[1].replace('"', ''))
                        continue
                    # ignore anything after a ';'
                    line = line.split(';')[0].strip()
                    if directive != 'system':  # no need to split line
                        tokens = [x.strip() for x in line.split()]
                    else:
                        tokens = [line]
                    padded = [typedata(k, _) for _, k in zip(tokens, header)]
                    # pad row so that it is same length as series
                    for _ in range(len(tokens), len(header)):
                        padded.append(pad)
                    # logger.debug(f'padded row: {padded}')
                    assert len(padded) == len(header), f'Error: Padding solution does not work! {directive} {len(tokens)}:{len(padded)}!={len(header)} {",".join(tokens)} {",".join(header)}'
                    for k, v in zip(header, padded):
                        series[k].append(v)
                tdf = pd.DataFrame(series)
                if directive == 'dihedraltypes':
                    if directive in inst.D:
                        # logger.info(f'Found second set of {len(tdf)} [ dihedraltypes ] in {inst.filename}; merging into set of {len(inst.D["dihedraltypes"])} types already read in...')
                        # we have already read-in a dihedraltypes section
                        # so let's append this one
                        inst.D['dihedraltypes'] = pd.concat([inst.D['dihedraltypes'], tdf], ignore_index=True).drop_duplicates()
                        # logger.info(f'    -> now there are {len(inst.D["dihedraltypes"])} dihedral types.')
                    else:
                        inst.D[directive] = tdf
                elif directive == 'dihedrals':
                    # if there is already a dihedrals section, assume
                    if directive in inst.D:
                        inst.D['dihedrals'] = pd.concat([inst.D['dihedrals'], tdf], ignore_index=True)
                    else:
                        inst.D[directive] = tdf
                else:
                    inst.D[directive] = tdf
            # we must assume the 'atoms' are sorted by global index; however, all other
            # sections need not be sorted.  For convenience, we will keep them sorted by
            # atom indices or atom type name, where appropriate.
            inst.null_check(msg=f'read from {filename}')
            if 'atomtypes' in inst.D:
                inst.D['atomtypes'] = inst.D['atomtypes'].sort_values(by='name').reset_index(drop=True)
            if 'bonds' in inst.D:
                inst.D['bonds'] = inst.D['bonds'].sort_values(by=['ai','aj']).reset_index(drop=True)
                inst.bondlist = Bondlist.fromDataFrame(inst.D['bonds'])
            if 'bondtypes' in inst.D:
                df_typeorder(inst.D['bondtypes'], typs=['i', 'j'])
                inst.D['bondtypes'] = inst.D['bondtypes'].sort_values(by=['i', 'j']).reset_index(drop=True)
            if 'pairs' in inst.D:
                inst.D['pairs'] = inst.D['pairs'].sort_values(by=['ai', 'aj']).reset_index(drop=True)
            if 'pairtypes' in inst.D:
                df_typeorder(inst.D['pairtypes'], typs=['i', 'j'])
                inst.D['pairtypes'] = inst.D['pairtypes'].sort_values(by=['i', 'j']).reset_index(drop=True)
            if 'angles' in inst.D:
                inst.D['angles'] = inst.D['angles'].sort_values(by=['aj', 'ai', 'ak']).reset_index(drop=True)
            if 'angletypes' in inst.D:
                df_typeorder(inst.D['angletypes'], typs=['i', 'j', 'k'])
                inst.D['angletypes'] = inst.D['angletypes'].sort_values(by=['j', 'i', 'k']).reset_index(drop=True)
            if 'dihedrals' in inst.D:
                inst.D['dihedrals'] = inst.D['dihedrals'].sort_values(by=['aj', 'ak', 'ai', 'al']).reset_index(drop=True)
            if 'dihedraltypes' in inst.D:
                df_typeorder(inst.D['dihedraltypes'], typs=['i', 'j', 'k', 'l'])
                inst.D['dihedraltypes'] = inst.D['dihedraltypes'].sort_values(by=['j', 'k', 'i', 'l']).reset_index(drop=True)
            for f in inst.includes:
                inst.merge(Topology.read_top(os.path.join(dirname, f)))
            inst.empty = False
            return inst

    def bond_source_check(self):
        """Checks that the 'bonds' and 'mol2_bonds' dataframes contain the same bonds."""
        if 'bonds' in self.D and 'mol2_bonds' in self.D:
            # logger.debug(f'Consistency check between gromacs-top bonds and mol2-bonds requested.')
            grobonds = self.D['bonds'].sort_values(by=['ai', 'aj'])
            bmi = grobonds.set_index(['ai', 'aj']).index
            mol2bonds = self.D['mol2_bonds'].sort_values(by=['ai', 'aj'])
            mbmi = mol2bonds.set_index(['ai', 'aj']).index
            check = all([x == y for x, y in zip(bmi, mbmi)])
            # logger.debug(f'Result: {check}')
            if not check:
                logger.error(f'Gromacs/Mol2 bond inconsistency detected')
                logger.debug(f'GROMACS:')
                logger.debug(grobonds[['ai','aj']].head().to_string())
                logger.debug(f'MOL2:')
                logger.debug(mol2bonds[['ai','aj']].head().to_string())
                for x, y in zip(bmi, mbmi):
                    logger.debug(f'{x} {y} {x == y}')

    def shiftatomsidx(self, idxshift, directive, rows=[], idxlabels=[]):
        """Shifts all atom indexes in a topology directive dataframe.

        Args:
            idxshift (int): integer index shift
            directive (str): name of GROMACS topology directive ('atoms','bonds','pairs','angles','dihedrals')
            rows (list): row boundaries [start, end], defaults to []
            idxlabels (list): names of columns that contain atom indexes, defaults to []
        """
        if directive in self.D:
            cols = self.D[directive].columns.get_indexer(idxlabels)
            self.D[directive].iloc[rows[0]:rows[1], cols] += idxshift

    def detect_rings(self):
        """Detects unique rings in the topology."""
        g = self.bondlist.graph()
        self.rings = RingList([])
        for c in nx.chordless_cycles(g):
            self.rings.append(Ring(c))

    def read_tpx(self, filename):
        if not os.path.exists(filename):
            # .tpx only holds the ring list, which we can rebuild from the bond
            # graph. This path matters for "Using cached parameterization" on a
            # molecule whose library had .gro/.top/.mol2 but no .tpx checked in.
            logger.debug(f'{filename} not found; detecting rings from topology graph.')
            self.detect_rings()
            return
        with open(filename, 'r') as f:
            data = f.read().split('[')
            stanzas = [('[' + x).split('\n') for x in data][1:]
            for s in stanzas:
                directive = s[0].split()[1].strip()
                contentlines = [l for l in s[1:] if not (len(l.strip()) == 0 or l.strip().startswith(';'))]
                if directive == 'rings':
                    # each line is a space-delimited list of global atom indices that identifies one ring
                    self.rings = RingList([])
                    for line in contentlines:
                        self.rings.append(Ring([int(x) for x in line.split()]))
                else:
                    logger.debug(f'tpx directive {directive} in {filename} is not recognized.')
    
    def write_tpx(self, filename):
        with open(filename, 'w') as f:
            f.write('[ rings ]\n')
            for r in self.rings:
                f.write(' '.join([str(x) for x in r.idx]) + '\n')

    def rep_ex(self, count=0):
        """Replicates extensive topology components (atoms, pairs, bonds, angles, dihedrals).

        Args:
            count (int): number of replicas to generate, defaults to 0

        Raises:
            Exception: if self is missing an atoms dataframe
        """
        if count > 0:
            counts = {k: 0 for k in _GromacsExtensiveDirectives_}
            counts.update({k: 0 for k in _NonGromacsExtensiveDirectives_})
            for t in _GromacsExtensiveDirectives_:
                if t in self.D:
                    counts[t] = len(self.D[t])
            for t in _NonGromacsExtensiveDirectives_:
                if t in self.D:
                    counts[t] = len(self.D[t])
            try:
                idxshift = counts['atoms']
            except:
                raise Exception(f'Error: expected an "atoms" dataframe')
            # Per-copy resnr span: shift each replica by this many residues so
            # multi-residue molecules don't collide on the same resid.
            resnr_per_copy = int(self.D['atoms']['resnr'].max())
            for t in _GromacsExtensiveDirectives_:
                if t in self.D:
                    self.D[t] = pd.concat([self.D[t]] * count, ignore_index=True)
            new_rings = RingList([])
            for c in range(1, count):
                self.shiftatomsidx(idxshift * c, 'atoms', rows=[(c * counts['atoms']), ((c + 1) * counts['atoms'])], idxlabels=['nr'])
                self.shiftatomsidx(resnr_per_copy * c, 'atoms', rows=[(c * counts['atoms']), ((c + 1) * counts['atoms'])], idxlabels=['resnr'])
                self.shiftatomsidx(idxshift * c, 'bonds', rows=[(c * counts['bonds']), ((c + 1) * counts['bonds'])], idxlabels=['ai', 'aj'])
                self.shiftatomsidx(idxshift * c, 'mol2_bonds', rows=[(c * counts['mol2_bonds']), ((c + 1) * counts['mol2_bonds'])], idxlabels=['ai', 'aj'])
                self.shiftatomsidx(idxshift * c, 'pairs', rows=[(c * counts['pairs']), ((c + 1) * counts['pairs'])], idxlabels=['ai', 'aj'])
                self.shiftatomsidx(idxshift * c, 'angles', rows=[(c * counts['angles']), ((c + 1) * counts['angles'])], idxlabels=['ai', 'aj', 'ak'])
                self.shiftatomsidx(idxshift * c, 'dihedrals', rows=[(c * counts['dihedrals']), ((c + 1) * counts['dihedrals'])], idxlabels=['ai', 'aj', 'ak', 'al'])
                new_ringblock = RingList([])
                for r in self.rings:
                    new_ringblock.append(r.copy().shift(idxshift * c))
                new_rings.extend(new_ringblock)
            self.rings.extend(new_rings)
            if 'bonds' in self.D:
                self.bondlist = Bondlist.fromDataFrame(self.D['bonds'])

    @classmethod
    def from_ex(cls,other):
        """Makes a new Topology instance by copying only the extensive dataframes from an existing topology.

        Args:
            other (Topology): the other topology

        Returns:
            Topology: a new Topology built from the extensive dataframes of other
        """
        ''' '''
        inst = cls()
        for t in _GromacsExtensiveDirectives_:
            if t in other.D:
                inst.D[t] = other.D[t].copy()
        if 'bonds' in inst.D:
            inst.bondlist = Bondlist.fromDataFrame(inst.D['bonds'])
        inst.rings = deepcopy(other.rings)
        return inst

    def write_top(self, filename):
        """Writes topology to a GROMACS-format top file.

        Args:
            filename (str): name of top file to write
        """
        self.null_check(msg=f'writing {filename}')
        assert 'defaults' in self.D, 'Error: no [ defaults ] in topology?'
        with open(filename, 'w') as f:
            f.write('; Gromacs-format topology written by HTPolyNet\n')
            for k in _GromacsTopologyDirectiveOrder_:
                if k not in self.D:
                    continue
                # Replace pad sentinels with NA so to_csv writes nothing for missing fields;
                # operate on a copy so self.D[k] is never mutated.
                df = self.D[k].replace(_PAD_, pd.NA)
                f.write(f'[ {k} ]\n; ')
                if k in _GromacsTopologyHashables_:
                    df = df.sort_values(by=_GromacsTopologyHashables_[k])
                # Cast integer-valued columns to nullable Int64 so they serialize
                # as e.g. "2" rather than "2.0" (parmed's gromacs reader insists
                # on parsing pn/func/atom-index columns as int).
                for col in _GromacsTopologyIntegerColumns_.get(k, []):
                    if col in df.columns:
                        df[col] = df[col].astype('Int64')
                df.to_csv(f, sep=' ', index=False, header=True, doublequote=False)
                f.write('\n')
            f.write('; end\n')

    def null_check(self, msg=''):
        """Checks for NaNs in dataframe locations that should never have NaNs.

        Args:
            msg (str): context message, defaults to ''

        Raises:
            Exception: if a NaN is found in a required field
        """
        for k in _GromacsTopologyDirectiveOrder_:
            if k in self.D and k in _GromacsTopologyHashables_:
                for a in _GromacsTopologyHashables_[k]:
                    has_nulls = any(self.D[k][a].isnull())
                    if has_nulls:
                        logger.debug(f'{msg} null in {k} {a}\n{self.D[k].to_string()}')
                        raise ValueError(f'{msg} NaN in {k}[{a}]\n{self.D[k].to_string()}')

    def total_charge(self):
        """Computes and returns total system charge.

        Returns:
            float: total system charge
        """
        if 'atoms' in self.D:
            return self.D['atoms']['charge'].sum()
        return 0.0

    def adjust_charges(self, atoms=[], desired_charge=0.0, overcharge_threshhold=0.1, msg=''):
        """Adjusts atom partial charges so that total system charge equals the desired value.

        Args:
            atoms (list): global atom indices to adjust, defaults to []
            desired_charge (float): target system charge, defaults to 0.0
            overcharge_threshhold (float): threshold overcharge that triggers a warning message, defaults to 0.1
            msg (str): message to write if pre-adjusted system charge is too high, defaults to ''

        Returns:
            Topology: self
        """
        apparent_charge = self.total_charge()
        overcharge = apparent_charge - desired_charge
        logger.debug(f'Adjusting charges of {len(atoms)} atoms due to overcharge of {overcharge:.6f}')
        if np.abs(overcharge) > overcharge_threshhold:
            logger.debug(f'{msg}')
        cpa = -overcharge / len(atoms)
        logger.debug(f'Adjustment is {cpa:.4e} per atom')
        for i in atoms:
            self.D['atoms'].loc[self.D['atoms']['nr'] == i, 'charge'] += cpa
        logger.debug(f'New total charge after adjustment: {self.total_charge():.6f}')
        return self
        
    def molecule_charges(self):
        """Returns the size and net charge of every covalently bonded molecule.

        Returns:
            list: (number of atoms, net charge) tuples, one per connected component of the bond
                graph, largest molecule first
        """
        if 'atoms' not in self.D or self.D['atoms'].empty:
            return []
        atoms = self.D['atoms']
        charge = dict(zip(atoms['nr'].astype(int), atoms['charge'].astype(float)))
        g = nx.Graph()
        g.add_nodes_from(charge)
        if 'bonds' in self.D and not self.D['bonds'].empty:
            g.add_edges_from(zip(self.D['bonds']['ai'].astype(int), self.D['bonds']['aj'].astype(int)))
        result = [(len(c), sum(charge[a] for a in c)) for c in nx.connected_components(g)]
        return sorted(result, key=lambda x: -x[0])

    def total_mass(self, units='gromacs'):
        """Returns total mass of all atoms in the Topology.

        Args:
            units (str): unit system; 'SI' returns kg, 'gromacs' returns amu, defaults to 'gromacs'

        Returns:
            float: total mass in amu (gromacs) or kg (SI)
        """
        fac = 1.0
        if units == 'SI':
            fac = physical_constants['atomic mass constant'][0]
        if 'atoms' in self.D:
            M_amu = self.D['atoms']['mass'].sum()
            logger.debug(f'mass {M_amu} fac {fac}')
            return M_amu * fac
        return 0.0

    def atomcount(self):
        """Returns the total number of atoms.

        Returns:
            int: number of atoms
        """
        if 'atoms' in self.D:
            return self.D['atoms'].shape[0]
        return 0
    
    def add_restraints(self, pairdf, typ=6, kb=300000.):
        """Adds type-6 (non-topological) bonds to help drag atoms toward each other.

        Args:
            pairdf (pandas.DataFrame): dataframe of pairs ['ai','aj']
            typ (int): bond type, defaults to 6
            kb (float): bond spring constant (kJ/mol/nm^2), defaults to 300000
        """
        bmi = self.D['bonds'].set_index(['ai','aj']).sort_index().index
        for i, b in pairdf.iterrows():
            ai, aj = idxorder((b['ai'],b['aj']))
            b0 = b['initial_distance']
            if not (ai, aj) in bmi:
                h = _GromacsTopologyDirectiveHeaders_['bonds']
                data = [ai, aj, typ, b0, kb]  # this new bond will have override parameters
                bonddict = {k: [v] for k, v in zip(h, data)}
                bdtoadd = pd.DataFrame(bonddict)
                self.D['bonds'] = pd.concat((self.D['bonds'], bdtoadd), ignore_index=True)

    def remove_restraints(self, pairdf):
        """Removes all non-topological restraint bonds represented in pairdf.

        Deleting these bonds does not influence angles or dihedrals.

        Args:
            pairdf (pandas.DataFrame): dataframe of pairs ['ai','aj']
        """
        d = self.D['bonds']
        to_drop = []
        for i, b in pairdf.iterrows():
            ai, aj = idxorder((b['ai'], b['aj']))
            to_drop.append(d[(d.ai == ai) & (d.aj == aj)].index.values[0])
        self.D['bonds'] = self.D['bonds'].drop(to_drop)

    def add_bonds(self, pairs=[]):
        """Adds bonds indicated in list pairs to the topology.

        Args:
            pairs (list): list of (ai, aj, order) tuples, defaults to []

        Raises:
            Exception: if an existing bond is in the list of pairs
        """
        # logger.debug('begins')
        at = self.D['atoms']
        ij = self.D['bondtypes'].set_index(['i', 'j'])
        #mb=self.D['mol2_bonds']
        bmi = self.D['bonds'].set_index(['ai', 'aj']).sort_index().index
        pmi = self.D['pairs'].set_index(['ai', 'aj']).sort_index().index
        newbonds = []
        for b in pairs:
            # logger.debug(f'{b}')
            # assert type(b[0])==int
            bondtuple = (int(b[0]), int(b[1]))
            order = int(b[2])
            ai, aj = idxorder(bondtuple)
            assert type(ai) == int
            assert type(aj) == int
            '''
            if this bond is not in the topology, then add it
            '''
            if not (ai, aj) in bmi:
                newbonds.append((ai, aj))
                # logger.debug(f'asking types of {ai} and {aj}; at.shape {at.shape}')
                it = at.iloc[ai-1].type
                jt = at.iloc[aj-1].type
                idx = typeorder((it, jt))
                if idx in ij.index:
                    bt = ij.loc[idx, 'func']  # why don't i need need values[0]
                    kb = ij.loc[idx, 'kb']
                    b0 = ij.loc[idx, 'b0']
                else:
                    logger.debug(f'no bondtype {idx} found; using placeholder parameters')
                    bt = 1
                    b0 = 0.15
                    kb = 999999
                    # raise Exception(f'no bondtype {idx} found.')
                '''
                add a new bond!
                '''
                h = _GromacsTopologyDirectiveHeaders_['bonds']
                data = [ai, aj, bt, b0, kb]  # this new bond will have override parameters
                assert len(h) == len(data), 'Error: not enough data for new bond?'
                bonddict = {k: [v] for k, v in zip(h, data)}
                bdtoadd = pd.DataFrame(bonddict)
                self.D['bonds'] = pd.concat((self.D['bonds'], bdtoadd), ignore_index=True)
                # logger.info(f'add_bond:\n{bdtoadd.to_string()}')
                logger.debug(f'just added {bonddict}')
                if 'mol2_bonds' in self.D:
                    data = [len(self.D['mol2_bonds']), ai, aj, '1']  # assume single bond; order is str dtype
                    bonddict = {k: [v] for k, v in zip(['bondIdx', 'ai', 'aj', 'order'], data)}
                    self.D['mol2_bonds'] = pd.concat((self.D['mol2_bonds'], pd.DataFrame(bonddict)),    ignore_index=True)
                # remove this pair from pairs if it's in there (it won't be)
                if idx in pmi:
                    logger.debug(f'Warning: new bond {ai}-{aj} was evidently in the [ pairs ]!')
                    d = self.D['pairs']
                    indexes_to_drop = d[(d.ai.isin(idx)) & (d.aj.isin(idx))].index
                    logger.debug(f'Dropping [ pair ]:\n{d[indexes_to_drop].to_string()}')
                    indexes_to_keep = set(range(d.shape[0])) - set(indexes_to_drop)
                    self.D['pairs'] = d.take(list(indexes_to_keep)).reset_index(drop=True)
            else:
                ''' 
                if it is, do nothing; it will be templated; if mol2_bonds are present (usually
                because a Topology is part of a molecule being parameterized), update the order
                of the bond.
                '''
                if 'mol2_bonds' in self.D:
                    mb = self.D['mol2_bonds']
                    bi = (mb['ai'] == ai) & (mb['aj'] == aj)
                    mb.loc[bi, 'order'] = str(order)
        '''
        update the bondlist
        '''
        for b in newbonds:
            self.bondlist.append(b)
        logger.debug(f'Added {len(newbonds)} new bonds')


    def delete_atoms(self, idx=[], reindex=True, return_idx_of=[], **kwargs):
        """Deletes atoms from topology.

        Args:
            idx (list): atom indexes to delete, defaults to []
            reindex (bool): reindex atoms after deleting, defaults to True
            return_idx_of (list): old indices whose new indices should be returned, defaults to []

        Returns:
            dict: old-index-to-new-index mapper
        """
        idx_set = set(idx)
        paranoid_about_pairs = kwargs.get('paranoid_about_pairs', False)
        self.null_check(msg='beginning of delete atoms')
        d = self.D['atoms']
        new_idx = []
        drop_mask = d['nr'].isin(idx_set)
        total_missing_charge = d.loc[drop_mask, 'charge'].sum()
        logger.debug(f'Deleting {drop_mask.sum()} [ atoms ]; charge to make up: {total_missing_charge:.4f}')
        self.D['atoms'] = d[~drop_mask].reset_index(drop=True)
        mapper = {}
        if reindex:
            d = self.D['atoms']
            oldGI = d['nr'].copy()
            d['nr'] = d.index + 1
            mapper = {k: v for k, v in zip(oldGI, d['nr'])}
            assert not any(x in mapper for x in idx_set), 'Error: Some deleted atoms in mapper.'
            if any(np.isnan(list(mapper.keys()))):
                logger.error('null in mapper keys')
            if any(np.isnan(list(mapper.values()))):
                logger.error('null in mapper values')
            if len(return_idx_of) > 0:
                new_idx = [mapper[o] for o in return_idx_of]
        for pt in ['bonds', 'mol2_bonds', 'pairs']:
            if pt not in self.D:
                continue
            d = self.D[pt]
            drop_mask = d['ai'].isin(idx_set) | d['aj'].isin(idx_set)
            logger.debug(f'Deleting {drop_mask.sum()} [ {pt} ]')
            self.D[pt] = d[~drop_mask].reset_index(drop=True)
            if reindex:
                d = self.D[pt]
                assert not any(x in idx_set for x in d.ai), f'Error: deleted atom survived in {pt} ai'
                assert not any(x in idx_set for x in d.aj), f'Error: deleted atom survived in {pt} aj'
                assert all(x in mapper for x in d.ai), f'Error: surviving {pt} atom ai old idx not in mapper'
                assert all(x in mapper for x in d.aj), f'Error: surviving {pt} atom aj old idx not in mapper'
                # Defer remapping pairs until after dihedrals are pruned.
                if pt != 'pairs':
                    d.ai = d.ai.map(mapper)
                    d.aj = d.aj.map(mapper)
                if pt == 'bonds':
                    self.bondlist = Bondlist.fromDataFrame(d)
                if pt == 'mol2_bonds':
                    self.D[pt]['bondIdx'] = list(range(1, self.D[pt].shape[0] + 1))
        d = self.D['angles']
        drop_mask = d['ai'].isin(idx_set) | d['aj'].isin(idx_set) | d['ak'].isin(idx_set)
        logger.debug(f'Deleting {drop_mask.sum()} [ angles ]')
        self.D['angles'] = d[~drop_mask].reset_index(drop=True)
        assert d.ai.dtype == int, f'post-delete lost angle ai dtype {d.ai.dtype}'
        assert d.aj.dtype == int, f'post-delete lost angle aj dtype {d.aj.dtype}'
        assert d.ak.dtype == int, f'post-delete lost angle ak dtype {d.ak.dtype}'
        self.null_check(msg='inside delete atoms before angles reindex')
        if reindex:
            d = self.D['angles']
            assert not any(x in idx_set for x in d.ai), 'Error: deleted atom survived in angle ai'
            assert not any(x in idx_set for x in d.aj), 'Error: deleted atom survived in angle aj'
            zombie_tags = [x in idx_set for x in d.ak]
            if any(zombie_tags):
                logger.debug(f'Zombie ak angles:\n{self.D["angles"][zombie_tags].to_string()}')
            assert not any(zombie_tags), 'Error: deleted atom survived in angle ak'
            assert all(x in mapper for x in d.ai), 'Error: surviving angle atom ai old idx not in mapper'
            assert all(x in mapper for x in d.aj), 'Error: surviving angle atom aj old idx not in mapper'
            assert all(x in mapper for x in d.ak), 'Error: surviving angle atom ak old idx not in mapper'
            d.ai = d.ai.map(mapper)
            d.aj = d.aj.map(mapper)
            d.ak = d.ak.map(mapper)
        self.null_check(msg='inside delete atoms after angles reindex')
        d = self.D['dihedrals']
        drop_mask = d['ai'].isin(idx_set) | d['aj'].isin(idx_set) | d['ak'].isin(idx_set) | d['al'].isin(idx_set)
        logger.debug(f'Deleting {drop_mask.sum()} [ dihedrals ]')
        # Optionally remove pairs whose defining dihedral (ai…al) is being deleted.
        if paranoid_about_pairs:
            dp = self.D['pairs']
            ddp = d[drop_mask]
            pair_drop = pd.Index([])
            for pi, pl in zip(ddp.ai, ddp.al):
                pair_drop = pair_drop.union(
                    dp[((dp.ai == pi) & (dp.aj == pl)) | ((dp.ai == pl) & (dp.aj == pi))].index
                )
            if len(pair_drop) > 0:
                logger.debug(f'  -> and deleting {len(pair_drop)} [ pairs ] from those dihedrals')
            self.D['pairs'] = dp.drop(index=pair_drop).reset_index(drop=True)
        self.D['dihedrals'] = d[~drop_mask].reset_index(drop=True)
        if reindex:
            d = self.D['dihedrals']
            d.ai = d.ai.map(mapper)
            d.aj = d.aj.map(mapper)
            d.ak = d.ak.map(mapper)
            d.al = d.al.map(mapper)
            d = self.D['pairs']
            d.ai = d.ai.map(mapper)
            d.aj = d.aj.map(mapper)
            # Drop any duplicate 1-4 pairs that survive after remapping.
            pdrops = d.duplicated(subset=['ai', 'aj'], keep='first')
            if pdrops.any():
                logger.debug(f'Deleting {pdrops.sum()} duplicate 1-4 pair descriptors -- this is likely due to a bug somewhere')
            self.D['pairs'] = d[~pdrops].reset_index(drop=True)
        self.null_check(msg='end of delete atoms')
        logger.debug('finished.')
        if len(return_idx_of) > 0:
            return new_idx
        return mapper

    def prune_stale_14_pairs(self):
        """Drops entries from [ pairs ] whose atoms are not actually 1-4 in the
        current bond graph.

        A template-derived 1-4 pair can become stale during cure when a
        subsequently formed bond shortens the topological path between its two
        atoms to 1-3 (or 1-2).  Such an entry incorrectly re-introduces a
        scaled nonbonded interaction on top of the bond/angle exclusion.

        Returns:
            int: number of pair entries dropped.
        """
        if 'pairs' not in self.D or self.D['pairs'].empty:
            return 0
        bl = self.bondlist
        keep_mask = []
        for ai, aj in zip(self.D['pairs']['ai'], self.D['pairs']['aj']):
            one_hop = set(bl.partners_of(int(ai)))
            if int(aj) in one_hop:
                keep_mask.append(False)
                continue
            two_hop = set()
            for p in one_hop:
                two_hop.update(bl.partners_of(int(p)))
            two_hop.discard(int(ai))
            two_hop -= one_hop
            if int(aj) in two_hop:
                keep_mask.append(False)
                continue
            three_hop = set()
            for p in two_hop:
                three_hop.update(bl.partners_of(int(p)))
            three_hop.discard(int(ai))
            three_hop -= one_hop
            three_hop -= two_hop
            keep_mask.append(int(aj) in three_hop)
        dropped = len(keep_mask) - sum(keep_mask)
        if dropped:
            self.D['pairs'] = self.D['pairs'][keep_mask].reset_index(drop=True)
        return dropped

    def _myconcat(self, other, directive='', idxlabel=[], idxshift=0, drop_duplicates=False):
        if not directive in other.D:
            return
        if directive in self.D:
            # shift atom indices
            for i in idxlabel:
                other.D[directive][i] += idxshift
            if drop_duplicates:
                self.D[directive] = pd.concat((self.D[directive], other.D[directive]), ignore_index=True).drop_duplicates()
            else:
                self.D[directive] = pd.concat((self.D[directive], other.D[directive]), ignore_index=True)
        else:
            self.D[directive] = other.D[directive]

    def merge(self, other):
        """Merges another topology into self.

        Args:
            other (Topology): topology to merge in
        """
        # logger.debug('Topology.merge begins')
        # look for duplicated types between self and other.  If any are found, delete those types from other and copy their parameters into the explicit interactions they correspond to.
        self.merge_ex(other)
        self.merge_types(other)
        # logger.debug('Topology.merge ends')

    # def handle_duplicate_types(self,other,copy_directive='other_to_self',typename='',funcidx=4):
#         if not typename in self.D or not typename in other.D:
#             return
#         stdf=self.D[typename]
#         otdf=other.D[typename]
#         hashables=_GromacsTopologyHashables_[typename]
#         headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
#         logger.debug(f'{hashables} {headers}')
#         for i in hashables:
#             headers.remove(i)
#         common=[]
#         fi=headers.index('func')
#         copy_idx_pairs=[]
#         for idx,r in otdf.iterrows():
#             typidx=typeorder(tuple([r[i] for i in hashables]))
#             idata=[r[i] for i in headers]
#             for jdx,q in stdf.iterrows():
#                 typjdx=typeorder(tuple([q[i] for i in hashables]))
#                 jdata=[q[i] for i in headers]
#                 if typidx==typjdx and idata==jdata:
#                     common.append(typidx)
#             for jdx,q in stdf.iterrows():
#                 typjdx=typeorder(tuple([q[i] for i in hashables]))
#                 jdata=[q[i] for i in headers]
#                 if typidx==typjdx and idata!=jdata and not typidx in common and idata[fi]==funcidx and jdata[fi]==funcidx:
#                     logger.debug(f'duplicate {typename} {typidx}')
#                     logger.debug(f'o {idx} {idata}')
#                     logger.debug(f's {jdx} {jdata}')
#                     if not (idx,jdx) in copy_idx_pairs:
#                         copy_idx_pairs.append((idx,jdx))
#         if copy_directive=='other_to_self' and len(copy_idx_pairs)>0:
#             for o,s in copy_idx_pairs:
#                 self.D[typename].iloc[s]=other.D[typename].iloc[o]
#         elif copy_directive=='self_to_other' and len(copy_idx_pairs)>0:
#             for o,s in copy_idx_pairs:
#                 other.D[typename].iloc[o]=self.D[typename].iloc[s]
#           # logger.debug(f'{typidx}')
# #        pass
    def report_type(self, typidx_q, typename='', funcidx=-1):
        if not typename in self.D:
            return
        stdf = self.D[typename]
        hashables = _GromacsTopologyHashables_[typename]
        headers = _GromacsTopologyDirectiveHeaders_[typename].copy()
        for i in hashables:
            headers.remove(i)
        fi = headers.index('func')
        for idx, r in stdf.iterrows():
            typidx=typeorder(tuple([r[i] for i in hashables]))
            idata=[r[i] for i in headers]
            if typidx==typidx_q and idata[fi]==funcidx:
                return [r[i] for i in headers]
        return []
    '''
    'bonds':        ['ai', 'aj', 'funct', 'c0', 'c1'],
    'bondtypes':    [ 'i',  'j', 'func',  'b0', 'kb'],
    'angles':       ['ai', 'aj', 'ak', 'funct', 'c0',  'c1'],
    'angletypes':   [ 'i',  'j',  'k', 'func',  'th0', 'cth', 'rub', 'kub'],
    'dihedrals':    ['ai', 'aj', 'ak', 'al', 'funct', 'c0',    'c1', 'c2', 'c3', 'c4', 'c5'],
    'dihedraltypes':[ 'i',  'j',  'k',  'l', 'func',  'phase', 'kd', 'pn'],
    '''
    def reset_override_from_type(self, interactionname, typename, inst_idx):
        if not typename in self.D:
            return
        typ_hashables = _GromacsTopologyHashables_[typename]
        stdf = self.D[typename].set_index(typ_hashables)
        sidf = self.D[interactionname]
        typ_headers = _GromacsTopologyDirectiveHeaders_[typename].copy()
        ins_hashables = _GromacsTopologyHashables_[interactionname]
        ins_headers = _GromacsTopologyDirectiveHeaders_[interactionname].copy()
        typidx = [self.D['atoms'].iloc[x - 1]['type'] for x in inst_idx]
        typidx = typeorder(tuple(typidx))
        for i in typ_hashables:
            typ_headers.remove(i)
        for i in ins_hashables:
            ins_headers.remove(i)
        iidx = idxorder(tuple(inst_idx))
        idx = -1
        for i, r in sidf.iterrows():
            jdx = tuple([r[a] for a in ins_hashables])
            if iidx == jdx:
                idx = i
                break
        assert idx != -1
        num_data = min([len(typ_headers), len(ins_headers)])
        typ_headers = typ_headers[:num_data]
        ins_headers = ins_headers[:num_data]
        typrec = stdf.loc[typidx][typ_headers]
        cols = self.D[interactionname].columns.get_indexer(ins_headers)
        logger.debug(f'Resetting override in {interactionname} for {inst_idx} from')
        for ln in self.D[interactionname].iloc[idx,cols].to_string().split('\n'):
            logger.debug(ln)
        logger.debug('to')
        self.D[interactionname].iloc[idx,cols] = typrec
        for ln in self.D[interactionname].iloc[idx,cols].to_string().split('\n'):
            logger.debug(ln)
        

    def reset_type(self, typename, typidx_t, values):
        if not typename in self.D:
            return
        stdf = self.D[typename]
        hashables = _GromacsTopologyHashables_[typename]
        headers = _GromacsTopologyDirectiveHeaders_[typename].copy()
        for i in hashables:
            headers.remove(i)
        cols = self.D[typename].columns.get_indexer(headers)
        idxs = []
        for idx, r in stdf.iterrows():
            typidx = typeorder(tuple([r[i] for i in hashables]))
            if typidx == typidx_t:
                idxs.append(idx)
        assert len(headers) == len(values)
        logger.debug(f'Resetting {len(idxs)} entries {headers} to {values}')
        for idx in idxs:
            self.D[typename].iloc[idx,cols] = values
            for ln in self.D[typename].iloc[idx][headers].to_string().split('\n'):
                logger.debug(ln)

    def report_duplicate_types(self,other,typename='',funcidx=4):
        if not typename in self.D or not typename in other.D:
            return
        stdf = self.D[typename]
        otdf = other.D[typename]
        hashables = _GromacsTopologyHashables_[typename]
        headers = _GromacsTopologyDirectiveHeaders_[typename].copy()
        # logger.debug(f'{hashables} {headers}')
        for i in hashables:
            headers.remove(i)
        common = []
        fi = headers.index('func')
        true_duplicate_types = []
        for idx, r in otdf.iterrows():
            typidx = typeorder(tuple([r[i] for i in hashables]))
            idata = [r[i] for i in headers]
            for jdx, q in stdf.iterrows():
                typjdx = typeorder(tuple([q[i] for i in hashables]))
                jdata = [q[i] for i in headers]
                if typidx == typjdx and idata == jdata:
                    common.append(typidx)
        for idx, r in otdf.iterrows():
            typidx = typeorder(tuple([r[i] for i in hashables]))
            idata = [r[i] for i in headers]
            for jdx, q in stdf.iterrows():
                typjdx = typeorder(tuple([q[i] for i in hashables]))
                jdata = [q[i] for i in headers]
                if typidx == typjdx and idata != jdata and not typidx in common and idata[fi] == funcidx and jdata[fi] == funcidx:
                    # logger.debug(f'duplicate {typename} {typidx}')
                    true_duplicate_types.append(typidx)
        return true_duplicate_types

    def dup_check(self, die=True):
        """Checks for duplicate type-like topology records.

        Args:
            die (bool): if True, raise an exception when a duplicate is found, defaults to True

        Raises:
            Exception: if a duplicate type with different parameters is detected and die is True
        """
        L = ['atomtypes', 'bondtypes', 'angletypes']
        Not = ' not' if not die else ''
        logger.debug(f'Checking for duplicate {L}; will{Not} die if found.')
        for t in L:
            i = _GromacsTopologyHashables_[t]
            ''' checking for types with duplicate atom-type indices '''
            dups = self.D[t].duplicated(subset=i, keep=False)
            if any(dups):
                logger.error(f'Duplicate {t} with different parameters detected\n' + self.D[t][dups].to_string())
                if die:
                    raise Exception('duplicate topology types with different parameters detected')

    def merge_types(self, other):
        """Merges type-like topology dataframes from other to self.

        Args:
            other (Topology): topology containing attribute D, a dictionary of dataframes
        """
        # self.handle_duplicate_types(other,typename='dihedraltypes',funcidx=4,drop_directive='drop_from_self')
        L = ['atomtypes', 'bondtypes', 'angletypes', 'dihedraltypes']
        for t in L:
            self._myconcat(other, directive=t, drop_duplicates=True)

    def merge_ex(self, other):
        logger.debug(f'   extensive merging...')
        ''' merge EXTENSIVE quantities '''
        idxshift = 0 if 'atoms' not in self.D else len(self.D['atoms'])
        rdxshift = 0 if 'atoms' not in self.D else self.D['atoms'].iloc[-1]['resnr']
        if 'atoms' in other.D:
            other.D['atoms']['resnr'] += rdxshift
        self._myconcat(other, directive='atoms', idxlabel=['nr'], idxshift=idxshift)
        self._myconcat(other, directive='bonds', idxlabel=['ai','aj'], idxshift=idxshift)
        if 'bonds' in other.D:
            self.bondlist.update(other.D['bonds'])
        if 'mol2_bonds' in self.D and 'mol2_bonds' in other.D:
            other.D['mol2_bonds'].bondIdx += self.D['mol2_bonds'].shape[0]
        self._myconcat(other, directive='mol2_bonds', idxlabel=['ai','aj'], idxshift=idxshift)
        self._myconcat(other, directive='pairs', idxlabel=['ai','aj'], idxshift=idxshift)
        self._myconcat(other, directive='angles', idxlabel=['ai','aj','ak'], idxshift=idxshift)
        self._myconcat(other, directive='dihedrals', idxlabel=['ai','aj','ak','al'], idxshift=idxshift)
        logger.debug(f'merging {len(other.rings)} rings into base list of {len(self.rings)} with idxshift {idxshift}')
        other.rings.shift(idxshift)
        self.rings.extend(other.rings)

    def get_atom_attribute(self, idx, attribute):
        """Returns value of attribute of atom idx.

        Args:
            idx (int): global atom index
            attribute (str): atom attribute name

        Returns:
            atom attribute value
        """
        return self.D['atoms'].loc[self.D['atoms']['nr'] == idx, attribute].values[0]

    def get_atomtype(self, idx):
        """Gets atom type of atom with global index idx.

        Args:
            idx (int): atom global index

        Returns:
            str: atom type
        """
        return self.D['atoms'].loc[self.D['atoms']['nr'] == idx, 'type'].values[0]

    def build_interresidue_graph(self, G, ri):
        if ri in G:
            return
        adf = self.D['atoms']
        ri_at_set = set(adf.loc[adf['resnr'] == ri, 'nr'])
        # Collect all bonded neighbours of ri's atoms that belong to other residues.
        neighbour_atoms = set()
        for a in ri_at_set:
            neighbour_atoms.update(self.bondlist.partners_of(a))
        neighbour_atoms -= ri_at_set
        ri_partners = adf.loc[adf['nr'].isin(neighbour_atoms), 'resnr'].unique().tolist()
        G.add_node(ri)
        for rj in ri_partners:
            self.build_interresidue_graph(G, rj)

    def local_resid_cluster(self,ri):
        G = nx.DiGraph()
        # logger.debug(str(self.bondlist))
        self.build_interresidue_graph(G, ri)
        return [x for x in G]

    def make_resid_graph(self, json_file=None):
        adf = self.D['atoms']
        N = adf.shape[0]
        self.residue_network = nx.DiGraph()
        residues = adf['residue'].unique()
        for rn in residues:
            rs = adf[adf['residue'] == rn]['resnr'].unique()
            self.residue_network.add_nodes_from(rs, resName=rn)
        
        resnrs = adf['resnr'].unique()
        for i in resnrs:
            # atom global indexes in this resnr
            ats = adf[adf['resnr'] == i]['nr'].to_list()
            connectors = []
            for a in ats:
                an = self.bondlist.partners_of(a)
                natsrn = adf[adf['nr'].isin(an)]['resnr'].unique()
                if len(natsrn) > 1:
                    for n in natsrn:
                        if n != i:  # this is an inter-residue connection
                            connectors.append((a, n))
            for c in connectors:
                a, n = c
                bondtype = "cross"
                if not self.residue_network.has_edge(i, n):
                    self.residue_network.add_edge(i, n, bondtype=bondtype)
        if json_file:
            the_data = json_graph.node_link_data(self.residue_network)
            logger.debug(f'writing graph node_link_data to {json_file}')
            with open(json_file, 'w') as f:
                json.dump(the_data, f, default=str)

    def copy_bond_parameters(self,bonds):
        """Generates and returns a copy of a bonds dataframe that contains all bonds listed in bonds.

        Args:
            bonds (pandas.DataFrame): dataframe of bonds managed by runtime, 'ai','aj','reactantName'

        Returns:
            pandas.DataFrame: bonds dataframe extracted from system with all parameters
        """
        bdf = self.D['bonds']
        assert not any(bdf['c0'].isna())
        pairs = [idxorder((b.ai, b.aj)) for b in bonds.itertuples()]
        idx = pd.MultiIndex.from_tuples(pairs, names=['ai', 'aj'])
        rows = []
        for ai, aj in pairs:
            rows.append(bdf[(bdf['ai'] == ai) & (bdf['aj'] == aj)].copy())
        saveme = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=bdf.columns)
        for c in bdf.columns:
            if c not in saveme:
                saveme[c] = _PAD_
        return saveme

    def attenuate_bond_parameters(self, bondsdf, stage, max_stages, minimum_distance=0.0, init_colname='initial_distance'):
        """Alters the kb and b0 parameters for new crosslink bonds according to the values prior to relaxation, their equilibrium values, and the ratio stage/max_stages.

        Let stage/max_stages be x, and 1/max_stages <= x <= 1. The spring constant for each
        bond is multiplied by x and the distance is 1/x of the way from its maximum value to
        its equilibrium value.

        Args:
            bondsdf (pandas.DataFrame): dataframe of bonds managed by runtime, 'ai','aj','reactantName'
            stage (int): index of stage in the series of post-bond-formation relaxation
            max_stages (int): total number of relaxation stages for this iteration
            minimum_distance (float): minimum bond length allowed, overriding type-specific b0 if greater than 0
        """
        bdf = self.D['bonds']
        factor = (stage + 1) / max_stages
        logger.debug(f'Attenuating {bondsdf.shape[0]} bond{"s" if bondsdf.shape[0] > 1 else ""} in stage {stage + 1}/{max_stages}')
        jdx = list(bondsdf.columns).index(init_colname)
        for b in bondsdf.itertuples(index=False):
            ai, aj = idxorder((b.ai, b.aj))
            rij = b[jdx]
            b0, kb = self.get_bond_parameters(ai, aj)
            if minimum_distance > 0.0:
                b0 = minimum_distance
            new_b0 = rij - factor * (rij - b0)
            new_kb = kb * factor
            # logger.debug(f'bond attenuation target for {ai}-{aj}:\nb0 {b0:.5f} kb {kb:.2f}; using b0 {new_b0:.5f} kb {new_kb:.2f}')
            bdf.loc[(bdf['ai'] == ai) & (bdf['aj'] == aj), 'c0'] = new_b0
            bdf.loc[(bdf['ai'] == ai) & (bdf['aj'] == aj), 'c1'] = new_kb

    def get_bond_parameters(self, ai, aj):
        """Gets b0 and kb for bond between atoms with global indexes ai and aj.

        Args:
            ai (int): global atom index
            aj (int): global atom index

        Returns:
            tuple: b0, kb — equilibrium bond length and spring constant
        """
        ai, aj = idxorder((ai, aj))
        adf = self.D['atoms']
        bdf = self.D['bonds']
        tdf = self.D['bondtypes']
        b0 = bdf.loc[(bdf['ai'] == ai) & (bdf['aj'] == aj), 'c0'].values[0]
        kb = bdf.loc[(bdf['ai'] == ai) & (bdf['aj'] == aj), 'c1'].values[0]
        if pd.isna(b0) or pd.isna(kb):
            ''' no overrides for this bond, so take from types '''
            it = adf.loc[adf['nr'] == ai, 'type'].values[0]
            jt = adf.loc[adf['nr'] == aj, 'type'].values[0]
            it, jt = typeorder((it, jt))
            b0 = tdf.loc[(tdf['i'] == it) & (tdf['j'] == jt), 'b0'].values[0]
            kb = tdf.loc[(tdf['i'] == it) & (tdf['j'] == jt), 'kb'].values[0]
        return b0, kb

    def restore_bond_parameters(self, df):
        """Copies data from all bonds in dataframe df to global dataframe.

        Args:
            df (pandas.DataFrame): dataframe of bonds ['ai','aj','c0','c1']
        """
        bdf = self.D['bonds']
        for r in df.itertuples():
            ai, aj = r.ai, r.aj
            c0, c1 = r.c0, r.c1
            bdf.loc[(bdf['ai'] == ai) & (bdf['aj'] == aj), 'c0'] = c0
            bdf.loc[(bdf['ai'] == ai) & (bdf['aj'] == aj), 'c1'] = c1

    def attenuate_pair_parameters(self, pairsdf, stage, max_stages, draglimit_nm=0.3):
        """Alters the kb and b0 parameters for new pre-crosslink pairs according to the values prior to dragging, the desired lower limit of interatomic distance, and the ratio stage/max_stages.

        Args:
            pairsdf (pandas.DataFrame): pairs dataframe (['ai'],['aj'],['initial_distance'])
            stage (int): index of stage in the series of pre-bond-formation dragging
            max_stages (int): total number of drag stages for this iteration
            draglimit_nm (float): lower limit of interatomic distance requested from drag
        """
        pdf = self.D['pairs']
        ess = 's' if pairsdf.shape[0] > 1 else ''
        factor = (stage + 1) / max_stages
        logger.debug(f'Attenuating {pairsdf.shape[0]} pair{ess} in stage {stage+1}/{max_stages}')
        for b in pairsdf.itertuples():
            ai, aj = idxorder((b.ai, b.aj))
            b0 = b.initial_distance
            kb = pdf.loc[(pdf['ai'] == ai) & (pdf['aj'] == aj), 'c1'].values[0]
            pdf.loc[(pdf['ai'] == ai) & (pdf['aj'] == aj), 'c0'] = draglimit_nm - factor * (b0 - draglimit_nm)
            pdf.loc[(pdf['ai'] == ai) & (pdf['aj'] == aj), 'c1'] = kb * factor
