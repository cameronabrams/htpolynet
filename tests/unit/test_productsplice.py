"""

.. module:: test_productsplice
   :synopsis: splicing a whole product template onto the residues of one reaction event,
              which is what a ring closure needs and a per-bond splice cannot do

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import pandas as pd
import pytest

from htpolynet.core.productsplice import atom_map, map_product_from_template
from htpolynet.core.topology import Topology


def atoms(rows):
    """rows: (nr, resnr, residue, atom, type, charge)"""
    return pd.DataFrame(rows, columns=['nr', 'resnr', 'residue', 'atom', 'type', 'charge'])


def bonds(rows):
    return pd.DataFrame(rows, columns=['ai', 'aj', 'funct', 'c0', 'c1'])


def angles(rows):
    return pd.DataFrame(rows, columns=['ai', 'aj', 'ak', 'funct', 'c0', 'c1'])


class Holder:
    """Stands in for TopoCoord or Molecule.TopoCoord: only .Topology is read."""
    def __init__(self, T):
        self.Topology = T


class Mol:
    def __init__(self, name, T):
        self.name, self.TopoCoord = name, Holder(T)


def cyanate_group(nr0, resnr, charges=(-0.24, 0.23, -0.16), types=('os', 'c1', 'n1')):
    """One -O-C#N group as three atoms O1, C1, N1, with its two internal bonds."""
    a = atoms([(nr0, resnr, 'BCY', 'O1', types[0], charges[0]),
               (nr0 + 1, resnr, 'BCY', 'C1', types[1], charges[1]),
               (nr0 + 2, resnr, 'BCY', 'N1', types[2], charges[2])])
    b = bonds([(nr0, nr0 + 1, 1, 0.136, 3.0e5),      # O-C
               (nr0 + 1, nr0 + 2, 1, 0.116, 8.0e5)])  # C#N, a triple bond
    ang = angles([(nr0, nr0 + 1, nr0 + 2, 1, 180.0, 500.0)])  # O-C#N, linear
    return a, b, ang


def instance(n_groups=3, first_resnr=5, extra_residue=True):
    """A system of unreacted cyanate groups, plus an untouched spectator residue."""
    A, B, G = [], [], []
    for k in range(n_groups):
        a, b, g = cyanate_group(1 + 3 * k, first_resnr + k)
        A.append(a); B.append(b); G.append(g)
    if extra_residue:
        nr0 = 1 + 3 * n_groups
        A.append(atoms([(nr0, 99, 'XXX', 'O1', 'oh', -0.40), (nr0 + 1, 99, 'XXX', 'C1', 'c3', 0.10)]))
        B.append(bonds([(nr0, nr0 + 1, 1, 0.140, 2.5e5)]))
    T = Topology()
    T.D['atoms'] = pd.concat(A, ignore_index=True)
    T.D['bonds'] = pd.concat(B, ignore_index=True)
    T.D['angles'] = pd.concat(G, ignore_index=True)
    return Holder(T)


def trimer_template():
    """Three cyanate groups with the triazine closed: ring bonds N(i)-C(i+1), aromatic
    types and charges, a ring angle, and the internal C-N now a ring bond."""
    A, B, G = [], [], []
    for k in range(3):
        a, b, g = cyanate_group(1 + 3 * k, 1 + k,
                                charges=(-0.30, 0.55, -0.45), types=('os', 'ca', 'nb'))
        b.loc[1, ['c0', 'c1']] = [0.134, 4.5e5]   # C#N became a ring bond
        g.loc[0, ['c0', 'c1']] = [120.0, 600.0]   # O-C-N is no longer linear
        A.append(a); B.append(b); G.append(g)
    ring, ring_ang = [], []
    for k in range(3):
        n = 3 + 3 * k              # N of group k
        c = 2 + 3 * ((k + 1) % 3)  # C of the next group
        ring.append((n, c, 1, 0.134, 4.5e5))
        ring_ang.append((2 + 3 * k, n, c, 1, 120.0, 600.0))   # C-N-C across the ring
    T = Topology()
    T.D['atoms'] = pd.concat(A, ignore_index=True)
    T.D['bonds'] = pd.concat(B + [bonds(ring)], ignore_index=True)
    T.D['angles'] = pd.concat(G + [angles(ring_ang)], ignore_index=True)
    return Mol('BCY3', T)


RESID_MAP = {1: 5, 2: 6, 3: 7}


def ring_bonds_in_instance():
    """The three instance bonds a cyclotrimerization would form, given instance()."""
    return [(3, 5), (6, 8), (9, 2)]


def form(TC, pairs):
    """Adds the new bonds with placeholder parameters, as make_bonds would."""
    rows = [(i, j, 1, 0.0, 0.0) for i, j in pairs]
    TC.Topology.D['bonds'] = pd.concat([TC.Topology.D['bonds'], bonds(rows)], ignore_index=True)


class TestAtomMap(unittest.TestCase):
    def test_maps_by_residue_and_atom_name(self):
        TC, tmpl = instance(), trimer_template()
        t2i, unmapped = atom_map(TC.Topology.D['atoms'], tmpl.TopoCoord.Topology.D['atoms'], RESID_MAP)
        self.assertEqual(len(t2i), 9)
        self.assertEqual(unmapped, [])
        # template residue 2's C1 is atom 5; instance residue 6's C1 is atom 5 as well
        self.assertEqual(t2i[5], 5)
        self.assertEqual(t2i[1], 1)

    def test_a_template_atom_with_no_counterpart_is_reported(self):
        TC, tmpl = instance(), trimer_template()
        tat = tmpl.TopoCoord.Topology.D['atoms']
        tmpl.TopoCoord.Topology.D['atoms'] = pd.concat(
            [tat, atoms([(99, 1, 'BCY', 'H9', 'ha', 0.1)])], ignore_index=True)
        t2i, unmapped = atom_map(TC.Topology.D['atoms'], tat, RESID_MAP)
        self.assertEqual(unmapped, [])
        t2i, unmapped = atom_map(TC.Topology.D['atoms'],
                                 tmpl.TopoCoord.Topology.D['atoms'], RESID_MAP)
        self.assertEqual(unmapped, [99])


class TestProductSplice(unittest.TestCase):
    def spliced(self):
        TC, tmpl = instance(), trimer_template()
        form(TC, ring_bonds_in_instance())
        stats = map_product_from_template(TC, tmpl, RESID_MAP, new_bonds=ring_bonds_in_instance())
        return TC, tmpl, stats

    def test_types_and_charges_come_from_the_template(self):
        TC, _, stats = self.spliced()
        a = TC.Topology.D['atoms'].set_index('nr')
        self.assertEqual(a.loc[2, 'type'], 'ca')
        self.assertEqual(a.loc[3, 'type'], 'nb')
        self.assertAlmostEqual(a.loc[2, 'charge'], 0.55)
        self.assertEqual(stats['types'], 6)   # C and N of three groups; O stays 'os'
        self.assertEqual(stats['charges'], 9)

    def test_the_spectator_residue_is_untouched(self):
        TC, _, stats = self.spliced()
        a = TC.Topology.D['atoms'].set_index('nr')
        self.assertEqual(a.loc[10, 'type'], 'oh')
        self.assertAlmostEqual(a.loc[10, 'charge'], -0.40)
        self.assertNotIn(10, stats['atoms'])

    def test_the_new_ring_bonds_get_template_parameters(self):
        TC, _, _ = self.spliced()
        b = TC.Topology.D['bonds']
        for i, j in ring_bonds_in_instance():
            row = b[((b.ai == i) & (b.aj == j)) | ((b.ai == j) & (b.aj == i))]
            self.assertEqual(len(row), 1)
            self.assertAlmostEqual(float(row.c0.iloc[0]), 0.134)
            self.assertAlmostEqual(float(row.c1.iloc[0]), 4.5e5)

    def test_an_intraresidue_bond_is_reparameterized_too(self):
        # the C#N triple of each group becomes a ring bond: a per-bond splice would
        # never touch it, because it contains no new bond
        TC, _, _ = self.spliced()
        b = TC.Topology.D['bonds']
        row = b[(b.ai == 2) & (b.aj == 3)]
        self.assertAlmostEqual(float(row.c0.iloc[0]), 0.134)
        self.assertAlmostEqual(float(row.c1.iloc[0]), 4.5e5)

    def test_an_intraresidue_angle_is_reparameterized(self):
        TC, _, _ = self.spliced()
        g = TC.Topology.D['angles']
        row = g[(g.ai == 1) & (g.aj == 2) & (g.ak == 3)]
        self.assertEqual(len(row), 1)
        self.assertAlmostEqual(float(row.c0.iloc[0]), 120.0)

    def test_ring_angles_are_added(self):
        TC, _, stats = self.spliced()
        g = TC.Topology.D['angles']
        self.assertEqual(len(g), 6)          # three internal, three across the ring
        self.assertEqual(stats['added']['angles'], 3)
        self.assertEqual(stats['replaced']['angles'], 3)

    def test_the_bondlist_is_rebuilt(self):
        TC, _, _ = self.spliced()
        self.assertIn(5, TC.Topology.bondlist.partners_of(3))

    def test_nothing_is_duplicated_when_spliced_twice(self):
        TC, tmpl, _ = self.spliced()
        before = (len(TC.Topology.D['bonds']), len(TC.Topology.D['angles']))
        map_product_from_template(TC, tmpl, RESID_MAP, new_bonds=ring_bonds_in_instance())
        self.assertEqual((len(TC.Topology.D['bonds']), len(TC.Topology.D['angles'])), before)


class TestProductSpliceRefusals(unittest.TestCase):
    def test_a_wrong_residue_map_is_refused(self):
        TC, tmpl = instance(), trimer_template()
        form(TC, ring_bonds_in_instance())
        with pytest.raises(ValueError, match='residue map is wrong'):
            map_product_from_template(TC, tmpl, {1: 5, 2: 6, 3: 99},
                                      new_bonds=ring_bonds_in_instance())

    def test_splicing_before_the_bonds_are_formed_is_refused(self):
        # the template has three ring bonds the instance lacks
        TC, tmpl = instance(), trimer_template()
        with pytest.raises(ValueError, match='the instance does not'):
            map_product_from_template(TC, tmpl, RESID_MAP)

    def test_strict_false_allows_it(self):
        TC, tmpl = instance(), trimer_template()
        stats = map_product_from_template(TC, tmpl, RESID_MAP, strict=False)
        self.assertEqual(stats['added']['bonds'], 3)

    def test_a_new_bond_the_template_lacks_is_refused(self):
        TC, tmpl = instance(), trimer_template()
        form(TC, [(1, 4)])      # O to O, not a template bond
        with pytest.raises(ValueError, match='not a bond of template'):
            map_product_from_template(TC, tmpl, RESID_MAP, new_bonds=[(1, 4)])
