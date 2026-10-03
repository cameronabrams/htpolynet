"""

.. module:: test_triazine_to_cyanate_cap
   :synopsis: runs the triazine-to-cyanate-cap repair driver end to end on a small synthetic system

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

The system is a real TopoCoord, read by the real readers from a .top/.gro
pair written to a temporary directory: TAZ-like triazines carrying 0, 1, 2
and 3 aryl-ether bonds to BPA-like bridges, plus free bridges, sized so the
number of unreacted bridge -OH equals the number of dangling ring carbons.

One thing is stubbed: ``topology_surgery.add_bonds_with_template``.  The real
one splices angles/dihedrals/pairs, types and charges from a parameterized
cap template molecule (map_from_templates over a MoleculeDict), which cannot
be had without AmberTools.  The stand-in forms the O-C bonds through the real
``TopoCoord.make_bonds``, writes cap-like types and charges onto the O, C and
N the template would map, and returns those atoms as the touched set.
Planning, bond deletion, residue reassignment, cap placement,
``refresh_bond_params``, H deletion with reindexing, and per-molecule
neutralization all run for real.
"""
import os
import unittest
import logging
import tempfile
from unittest import mock

import networkx as nx
import numpy as np

from htpolynet.core.topocoord import TopoCoord
from htpolynet.geometry.bondlist import Bondlist
from htpolynet.repair import topology_surgery as ts
from htpolynet.repair.cyanate_cap import triazine_to_cyanate_cap

logger = logging.getLogger(__name__)

SPEC = {
    'type': 'triazine_to_cyanate_cap',
    'crosslinker': {'residue': 'TAZ', 'full_bond_count': 3,
                    'ring_carbon_atoms': ['C1', 'C2', 'C3'],
                    'ring_nitrogen_atoms': ['N1', 'N2', 'N3']},
    'bridge': {'residue': 'BPA', 'reactive_oxygen_atoms': ['O1', 'O2']},
    'cap_residue': 'CYN',
    'cap_template': 'BPA~CYN',
    'cap_search_radius': 0.6,
}

MASS = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'H': 1.008}
CHARGE = {'ca': 0.10, 'nb': -0.30, 'h4': 0.15, 'os': -0.35, 'oh': -0.50, 'ho': 0.40, 'c3': -0.05}
BONDTYPES = [('ca', 'nb', 0.134), ('ca', 'h4', 0.108), ('ca', 'os', 0.136), ('ca', 'oh', 0.136),
             ('ca', 'c3', 0.151), ('oh', 'ho', 0.097), ('c1', 'n1', 0.1155), ('c1', 'os', 0.131),
             ('ca', 'n1', 0.130), ('ca', 'c1', 0.140)]
CN_B0 = 0.1155
# template-like charges the stand-in writes on a spliced cap
CAP_Q = {'o': -0.30, 'c': 0.50, 'n': -0.45}

R_RING = 0.14
Z = np.array([0.0, 0.0, 1.0])


class _Builder:
    """Accumulates atoms and bonds; residues are numbered as they are added."""
    def __init__(self, reverse_ring_bonds=False):
        self.reverse_ring_bonds = reverse_ring_bonds
        self.atoms = []   # dicts: name, res, resnr, type, xyz
        self.bonds = []
        self.resnr = 0

    def atom(self, name, res, typ, xyz):
        self.atoms.append({'name': name, 'res': res, 'resnr': self.resnr, 'type': typ,
                           'xyz': np.asarray(xyz, dtype=float)})
        return len(self.atoms)

    def bond(self, i, j, reverse=False):
        # ai < aj is the usual convention; reverse stores (higher, lower)
        self.bonds.append((max(i, j), min(i, j)) if reverse else (min(i, j), max(i, j)))

    def bpa(self, origin, u, o1_reacted):
        """O1-C1-C3-C2-O2 along u from origin; O2 always carries an H, O1 only if unreacted."""
        self.resnr += 1
        u = np.asarray(u, dtype=float)
        o1 = self.atom('O1', 'BPA', 'os' if o1_reacted else 'oh', origin)
        c1 = self.atom('C1', 'BPA', 'ca', origin + 0.137 * u)
        c3 = self.atom('C3', 'BPA', 'c3', origin + 0.289 * u)
        c2 = self.atom('C2', 'BPA', 'ca', origin + 0.441 * u)
        o2 = self.atom('O2', 'BPA', 'oh', origin + 0.577 * u)
        for i, j in ((o1, c1), (c1, c3), (c3, c2), (c2, o2)):
            self.bond(i, j)
        # H at a C-O-H angle of 110 degrees, off the chain axis
        if not o1_reacted:
            h1 = self.atom('H1', 'BPA', 'ho', origin + 0.097 * (-0.342 * u + 0.940 * Z))
            self.bond(o1, h1)
        h2 = self.atom('H2', 'BPA', 'ho', origin + 0.577 * u + 0.097 * (0.342 * u + 0.940 * Z))
        self.bond(o2, h2)
        return o1

    def taz(self, center, k):
        """Triazine C1-N1-C2-N2-C3-N3 in the xy plane; the first k ring C's take a BPA."""
        self.resnr += 1
        resnr = self.resnr
        center = np.asarray(center, dtype=float)
        ring_c, ring_n, units = [], [], []
        for i in range(3):
            tc = np.radians(120.0 * i)
            tn = np.radians(120.0 * i + 60.0)
            u = np.array([np.cos(tc), np.sin(tc), 0.0])
            units.append(u)
            ring_c.append(self.atom(f'C{i+1}', 'TAZ', 'ca', center + R_RING * u))
            ring_n.append(self.atom(f'N{i+1}', 'TAZ', 'nb',
                                    center + R_RING * np.array([np.cos(tn), np.sin(tn), 0.0])))
        for i in range(3):
            self.bond(ring_c[i], ring_n[i], reverse=self.reverse_ring_bonds)
            self.bond(ring_n[i], ring_c[(i + 1) % 3], reverse=self.reverse_ring_bonds)
        for i in range(3):
            if i >= k:
                h = self.atom(f'H{i+1}', 'TAZ', 'h4', center + (R_RING + 0.108) * units[i])
                self.bond(ring_c[i], h)
        for i in range(k):
            o1 = self.bpa(center + (R_RING + 0.136) * units[i], units[i], o1_reacted=True)
            self.bond(ring_c[i], o1)
        return resnr

    def write(self, dirname, box):
        """Write top/gro; charges make every molecule neutral, as a real build would be."""
        n = len(self.atoms)
        q = np.array([CHARGE[a['type']] for a in self.atoms])
        g = nx.Graph()
        g.add_nodes_from(range(1, n + 1))
        g.add_edges_from(self.bonds)
        for comp in nx.connected_components(g):
            first = min(comp)
            q[first - 1] -= sum(q[i - 1] for i in comp)
        angles = set()
        for j in g.nodes:
            nb = sorted(g.neighbors(j))
            for a in range(len(nb)):
                for b in range(a + 1, len(nb)):
                    angles.add((nb[a], j, nb[b]))
        dihedrals = set()
        for j, k in g.edges:
            for i in g.neighbors(j):
                if i == k:
                    continue
                for l in g.neighbors(k):
                    if l in (j, i):
                        continue
                    dihedrals.add((i, j, k, l) if j < k else (l, k, j, i))
        pairs = {(min(d[0], d[3]), max(d[0], d[3])) for d in dihedrals}
        top = os.path.join(dirname, 'sys.top')
        gro = os.path.join(dirname, 'sys.gro')
        with open(top, 'w') as f:
            f.write('[ defaults ]\n1 2 yes 0.5 0.83333333\n\n[ atomtypes ]\n')
            for t, e in (('ca', 'C'), ('c3', 'C'), ('c1', 'C'), ('nb', 'N'), ('n1', 'N'),
                         ('os', 'O'), ('oh', 'O'), ('ho', 'H'), ('h4', 'H')):
                f.write(f'{t} {dict(C=6, N=7, O=8, H=1)[e]} {MASS[e]} 0.0 A 0.3 0.4\n')
            f.write('\n[ bondtypes ]\n')
            for i, j, b0 in BONDTYPES:
                f.write(f'{i} {j} 1 {b0} 300000.0\n')
            f.write('\n[ moleculetype ]\nwhole_system 3\n\n[ atoms ]\n')
            for i, a in enumerate(self.atoms, start=1):
                f.write(f'{i} {a["type"]} {a["resnr"]} {a["res"]} {a["name"]} {i} '
                        f'{q[i-1]:.8f} {MASS[a["name"][0]]}\n')
            f.write('\n[ pairs ]\n')
            for i, j in sorted(pairs):
                f.write(f'{i} {j} 1\n')
            f.write('\n[ bonds ]\n')
            for i, j in sorted(self.bonds):
                f.write(f'{i} {j} 1 0.14 300000.0\n')
            f.write('\n[ angles ]\n')
            for i, j, k in sorted(angles):
                f.write(f'{i} {j} {k} 1 120.0 500.0\n')
            f.write('\n[ dihedrals ]\n')
            for i, j, k, l in sorted(dihedrals):
                f.write(f'{i} {j} {k} {l} 9 180.0 10.0 2\n')
            f.write('\n[ system ]\nsynthetic\n\n[ molecules ]\nwhole_system 1\n')
        with open(gro, 'w') as f:
            f.write(f'synthetic\n{n}\n')
            for i, a in enumerate(self.atoms, start=1):
                x, y, zz = a['xyz']
                f.write(f'{a["resnr"]:5d}{a["res"]:<5s}{a["name"]:>5s}{i:5d}{x:8.3f}{y:8.3f}{zz:8.3f}\n')
            f.write(f'{box:10.5f}{box:10.5f}{box:10.5f}\n')
        return top, gro


def build_system(dirname, k_list, n_free_bpa, strip_one_oh=False, box=6.0, reverse_ring_bonds=False):
    """Return a TopoCoord with one TAZ per entry of k_list and n_free_bpa free bridges.

    With strip_one_oh, the phenolic H on the first bonded bridge's O2 is left
    off, leaving one more dangling ring C than unreacted bridge -OH.
    """
    B = _Builder(reverse_ring_bonds=reverse_ring_bonds)
    centers = [(1.0 + 2.0 * (i % 3), 1.0 + 2.0 * ((i // 3) % 3), 1.5 + 3.0 * (i // 9))
               for i in range(len(k_list))]
    for c, k in zip(centers, k_list):
        B.taz(c, k)
    for m in range(n_free_bpa):
        B.bpa(np.array([0.8 + 2.0 * m, 1.0, 4.5]), np.array([1.0, 0.0, 0.0]), o1_reacted=False)
    if strip_one_oh:
        o2 = next(i for i, a in enumerate(B.atoms, start=1) if a['res'] == 'BPA' and a['name'] == 'O2')
        h2 = next(j if i == o2 else i for i, j in B.bonds if o2 in (i, j)
                  and B.atoms[(j if i == o2 else i) - 1]['name'].startswith('H'))
        # keep the H as an atom bonded to nothing would change the H census; drop it
        B.bonds = [b for b in B.bonds if h2 not in b]
        B.bonds = [(i - (i > h2), j - (j > h2)) for i, j in B.bonds]
        del B.atoms[h2 - 1]
    top, gro = B.write(dirname, box)
    TC = TopoCoord(topfilename=top, grofilename=gro)
    # cure-stage extended attributes the driver resets on repaired atoms
    A = TC.Coordinates.A
    A['z'] = 0
    A['nreactions'] = 0
    return TC


def splice_standin(calls):
    """Stand-in for topology_surgery.add_bonds_with_template; see module docstring."""
    def standin(TC, pairs, moldict, product_name, chain_manager=None, adjust_charges=True):
        calls.append({'pairs': list(pairs), 'moldict': moldict, 'product_name': product_name,
                      'adjust_charges': adjust_charges})
        if not pairs:
            return []
        TC.make_bonds(pairs, explicit_sacH={i: [] for i in range(len(pairs))},
                      chain_manager=chain_manager)
        touched = []
        for o, c, _ in pairs:
            n = [p for p in TC.Topology.bondlist.partners_of(c)
                 if ts.get_attr(TC, p, 'resName') == 'CYN' and ts.get_attr(TC, p, 'atomName') == 'N1']
            assert len(n) == 1, f'cap C {c} has {len(n)} cap N partners'
            n = n[0]
            ts.set_atom_attributes(TC, o, type='os', charge=CAP_Q['o'])
            ts.set_atom_attributes(TC, c, type='c1', charge=CAP_Q['c'])
            ts.set_atom_attributes(TC, n, type='n1', charge=CAP_Q['n'])
            touched += [o, c, n]
        ts._fix_atom_index_dtypes(TC)
        return touched
    return standin


def census(TC):
    A = TC.Coordinates.A
    el = A['atomName'].str[0]
    return {e: int((el == e).sum()) for e in 'CNOH'}


def run_repair(TC, calls=None):
    calls = [] if calls is None else calls
    moldict = object()
    with mock.patch.object(ts, 'add_bonds_with_template', splice_standin(calls)):
        result = triazine_to_cyanate_cap(TC, moldict, SPEC, reactions=None)
    return result, moldict


def bond_set(TC):
    b = TC.Topology.D['bonds']
    return {(min(int(i), int(j)), max(int(i), int(j))) for i, j in zip(b['ai'], b['aj'])}


def residues(TC, name):
    A = TC.Coordinates.A
    sub = A[A['resName'] == name]
    return {int(rn): {r.atomName: int(r.globalIdx) for r in g.itertuples()}
            for rn, g in sub.groupby('resNum')}


class TestTriazineToCyanateCap(unittest.TestCase):
    # k distribution: dangling ring C = 2*3 + 2*2 + 1 = 11; bonded bridge O2-H = 2+2+3 = 7,
    # plus 2 free bridges with 2 OH each = 11 unreacted -OH
    K_LIST = [0, 0, 1, 1, 2, 3]
    N_FREE_BPA = 2

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        TC = build_system(cls._tmp.name, cls.K_LIST, cls.N_FREE_BPA)
        cls.before = census(TC)
        cls.n_before = len(TC.Coordinates.A)
        cls.dangling_ch = len(cls.K_LIST) * 3 - sum(cls.K_LIST)
        A = TC.Coordinates.A
        bl = TC.Topology.bondlist
        cls.unreacted_oh = sum(
            1 for r in A[(A['resName'] == 'BPA') & A['atomName'].isin(['O1', 'O2'])].itertuples()
            if any(ts.get_attr(TC, p, 'atomName').startswith('H') for p in bl.partners_of(int(r.globalIdx))))
        cls.calls = []
        cls.result, cls.moldict = run_repair(TC, cls.calls)
        cls.TC = TC

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_system_is_balanced_as_built(self):
        self.assertEqual(self.dangling_ch, 11)
        self.assertEqual(self.unreacted_oh, self.dangling_ch)

    def test_heavy_atoms_are_conserved(self):
        after = census(self.TC)
        for e in 'CNO':
            self.assertEqual(after[e], self.before[e], f'{e} count changed')

    def test_hydrogens_removed_are_ring_ch_plus_oh(self):
        after = census(self.TC)
        self.assertEqual(self.before['H'] - after['H'], self.dangling_ch + self.unreacted_oh)
        self.assertEqual(self.n_before - len(self.TC.Coordinates.A), self.dangling_ch + self.unreacted_oh)

    def test_complete_triazine_survives_intact(self):
        taz = residues(self.TC, 'TAZ')
        self.assertEqual(len(taz), self.K_LIST.count(3))
        bonds = bond_set(self.TC)
        bl = self.TC.Topology.bondlist
        for atoms in taz.values():
            self.assertEqual(sorted(atoms), ['C1', 'C2', 'C3', 'N1', 'N2', 'N3'])
            c = [atoms[f'C{i}'] for i in (1, 2, 3)]
            n = [atoms[f'N{i}'] for i in (1, 2, 3)]
            for i in range(3):
                self.assertIn((min(c[i], n[i]), max(c[i], n[i])), bonds)
                nxt = c[(i + 1) % 3]
                self.assertIn((min(n[i], nxt), max(n[i], nxt)), bonds)
                o = [p for p in bl.partners_of(c[i]) if ts.get_attr(self.TC, p, 'resName') == 'BPA']
                self.assertEqual(len(o), 1)

    def test_incomplete_triazines_become_three_caps_each(self):
        n_incomplete = sum(1 for k in self.K_LIST if k < 3)
        cyn = residues(self.TC, 'CYN')
        self.assertEqual(len(cyn), 3 * n_incomplete)
        A = self.TC.Coordinates.A
        bl = self.TC.Topology.bondlist
        bonds = self.TC.Topology.D['bonds']
        atoms = self.TC.Topology.D['atoms']
        o_seen = set()
        for atoms_by_name in cyn.values():
            self.assertEqual(sorted(atoms_by_name), ['C1', 'N1'])
            c, n = atoms_by_name['C1'], atoms_by_name['N1']
            self.assertTrue(bl.are_bonded(c, n))
            self.assertEqual(bl.partners_of(n), [c])
            o = [p for p in bl.partners_of(c) if p != n]
            self.assertEqual(len(o), 1)
            row = A[A['globalIdx'] == o[0]].iloc[0]
            self.assertEqual(row['resName'], 'BPA')
            self.assertIn(row['atomName'], ('O1', 'O2'))
            o_seen.add(o[0])
            # the C-N bond carries the triple-bond (c1-n1) parameters after refresh
            self.assertEqual(atoms.loc[atoms['nr'] == c, 'type'].iloc[0], 'c1')
            m = (((bonds['ai'] == c) & (bonds['aj'] == n)) | ((bonds['ai'] == n) & (bonds['aj'] == c)))
            self.assertEqual(int(m.sum()), 1)
            self.assertAlmostEqual(float(bonds.loc[m, 'c0'].iloc[0]), CN_B0)
        self.assertEqual(len(o_seen), len(cyn), 'two caps share one bridge O')

    def test_no_bridge_oxygen_keeps_its_h(self):
        A = self.TC.Coordinates.A
        bl = self.TC.Topology.bondlist
        for r in A[(A['resName'] == 'BPA') & A['atomName'].isin(['O1', 'O2'])].itertuples():
            names = [ts.get_attr(self.TC, p, 'atomName') for p in bl.partners_of(int(r.globalIdx))]
            self.assertFalse(any(x.startswith('H') for x in names), f'O {r.globalIdx} still has an H')
            self.assertEqual(len(names), 2)

    def test_no_ring_h_survives(self):
        A = self.TC.Coordinates.A
        self.assertEqual(int(((A['resName'] == 'TAZ') & A['atomName'].str.startswith('H')).sum()), 0)

    def test_statistics(self):
        r = self.result
        hist = {k: self.K_LIST.count(k) for k in range(4)}
        self.assertEqual(r['prerepair_bond_counts'], hist)
        self.assertEqual(r['n_crosslinkers'], len(self.K_LIST))
        self.assertEqual(r['n_dismantled'], sum(1 for k in self.K_LIST if k < 3))
        self.assertEqual(r['n_complete'], self.K_LIST.count(3))
        self.assertEqual(r['n_transferred'], self.dangling_ch)

    def test_splice_called_once_without_charge_adjustment(self):
        self.assertEqual(len(self.calls), 1)
        call = self.calls[0]
        self.assertIs(call['moldict'], self.moldict)
        self.assertEqual(call['product_name'], SPEC['cap_template'])
        self.assertFalse(call['adjust_charges'])
        self.assertEqual(len(call['pairs']), 3 * sum(1 for k in self.K_LIST if k < 3))

    def test_every_molecule_is_neutral(self):
        T = self.TC.Topology
        self.assertAlmostEqual(float(T.total_charge()), 0.0, places=6)
        atoms = T.D['atoms']
        q = dict(zip(atoms['nr'].astype(int), atoms['charge'].astype(float)))
        g = Bondlist.fromDataFrame(T.D['bonds']).graph()
        g.add_nodes_from(q)
        for comp in nx.connected_components(g):
            self.assertAlmostEqual(sum(q[i] for i in comp), 0.0, places=6)

    def test_indices_are_consistent_after_reindexing(self):
        n = len(self.TC.Coordinates.A)
        A = self.TC.Coordinates.A
        atoms = self.TC.Topology.D['atoms']
        self.assertEqual(A['globalIdx'].tolist(), list(range(1, n + 1)))
        self.assertEqual(atoms['nr'].tolist(), list(range(1, n + 1)))
        self.assertEqual(A['resName'].tolist(), atoms['residue'].tolist())
        self.assertEqual(A['atomName'].tolist(), atoms['atom'].tolist())
        bonds = bond_set(self.TC)
        for i, j in bonds:
            self.assertTrue(1 <= i <= n and 1 <= j <= n)
        # every surviving angle and dihedral runs along bonds that still exist
        for r in self.TC.Topology.D['angles'].itertuples():
            for a, b in ((r.ai, r.aj), (r.aj, r.ak)):
                self.assertIn((min(a, b), max(a, b)), bonds)
        for r in self.TC.Topology.D['dihedrals'].itertuples():
            for a, b in ((r.ai, r.aj), (r.aj, r.ak), (r.ak, r.al)):
                self.assertIn((min(a, b), max(a, b)), bonds)

    def test_repaired_atoms_are_no_longer_reactive(self):
        A = self.TC.Coordinates.A
        cyn = A[A['resName'] == 'CYN']
        self.assertTrue((cyn['z'] == 0).all())
        self.assertTrue((cyn['nreactions'] == 1).all())

    def test_caps_sit_at_bond_length_from_their_oxygen(self):
        # in-place caps keep their ring-C position, 0.136 from O by construction;
        # transferred ones are placed at oc_len=0.136 with C#N 0.116
        bl = self.TC.Topology.bondlist
        n_short_cn = 0
        for atoms_by_name in residues(self.TC, 'CYN').values():
            c, n = atoms_by_name['C1'], atoms_by_name['N1']
            o = [p for p in bl.partners_of(c) if p != n][0]
            pc, pn, po = ts.positions(self.TC, [c, n, o])
            self.assertAlmostEqual(np.linalg.norm(pc - po), 0.136, delta=0.002)
            d_cn = np.linalg.norm(pn - pc)
            self.assertTrue(abs(d_cn - 0.116) < 0.002 or abs(d_cn - R_RING) < 0.002, d_cn)
            n_short_cn += abs(d_cn - 0.116) < 0.002
        self.assertEqual(n_short_cn, self.dangling_ch)
        self.assertEqual(self.result['n_below_target'], 0)
        self.assertGreaterEqual(self.result['min_clearance_nm'], 0.15)


class TestFullyCuredIsUntouched(unittest.TestCase):
    def test_all_complete_is_a_no_op(self):
        with tempfile.TemporaryDirectory() as d:
            TC = build_system(d, [3, 3], 0)
            n = len(TC.Coordinates.A)
            calls = []
            result, _ = run_repair(TC, calls)
            self.assertEqual(len(TC.Coordinates.A), n)
            self.assertEqual(calls, [])
            self.assertEqual(result['n_dismantled'], 0)
            self.assertEqual(result['prerepair_bond_counts'], {0: 0, 1: 0, 2: 0, 3: 2})


class TestConservationMismatch(unittest.TestCase):
    def test_extra_dangling_carbon_is_reported_and_its_fragment_discarded(self):
        # balanced k = 0,1,2,3 system (6 dangling, 6 -OH) with one -OH H left off
        with tempfile.TemporaryDirectory() as d:
            TC = build_system(d, [0, 1, 2, 3], 0, strip_one_oh=True)
            before = census(TC)
            with self.assertLogs('htpolynet.repair.cyanate_cap', level='WARNING') as cm:
                run_repair(TC)
            self.assertTrue(any('atom-conservation mismatch' in m for m in cm.output))
            after = census(TC)
            # the homeless fragment's C and N are deleted, nothing else heavy
            self.assertEqual(before['C'] - after['C'], 1)
            self.assertEqual(before['N'] - after['N'], 1)
            self.assertEqual(before['O'], after['O'])
            self.assertEqual(len(residues(TC, 'CYN')), 3 * 3 - 1)


class TestReversedRingBonds(unittest.TestCase):
    def test_cap_cn_bonds_refresh_when_stored_higher_lower(self):
        # ring bonds written as (higher, lower); the refresh must still find
        # each cap's C-N bond and give it the c1-n1 parameters
        with tempfile.TemporaryDirectory() as d:
            TC = build_system(d, [0, 1, 2, 3], 0, reverse_ring_bonds=True)
            b = TC.Topology.D['bonds']
            self.assertTrue((b['ai'] > b['aj']).any())
            with self.assertNoLogs('htpolynet.repair.topology_surgery', level='WARNING'):
                run_repair(TC)
            bonds = TC.Topology.D['bonds']
            cyn = residues(TC, 'CYN')
            self.assertEqual(len(cyn), 9)
            for atoms_by_name in cyn.values():
                c, n = atoms_by_name['C1'], atoms_by_name['N1']
                m = (((bonds['ai'] == c) & (bonds['aj'] == n)) | ((bonds['ai'] == n) & (bonds['aj'] == c)))
                self.assertEqual(int(m.sum()), 1)
                self.assertAlmostEqual(float(bonds.loc[m, 'c0'].iloc[0]), CN_B0)
