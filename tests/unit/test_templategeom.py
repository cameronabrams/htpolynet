"""

.. module:: test_templategeom
   :synopsis: making a product template a plausible molecule before its charges are
              computed from its geometry

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import types
import unittest

import pandas as pd
import pytest

from htpolynet.core.topology import Topology, _element_of
from htpolynet.geometry.templategeom import (check_ring_closed, geometry_complaints,
                                             relax_geometry)


def cyanate_pair(ring_bond_nm=0.223, ring_order='1', cn_order='3'):
    """Two methyl cyanates with one ring-forming bond recorded between them.

    Enough of a template to exercise valence and bond-length checking: each cyanate is
    H3C-O-C#N, and the recorded bond runs from one's N to the other's C.  That leaves the
    carbon pentavalent and the nitrogen over-valent while both C#N stand, which is exactly
    what an addition does.  The two units point away from each other so the geometry is
    non-degenerate and a force field can actually run on it.
    """
    r = ring_bond_nm
    # (name, type, position in nm); unit 1 extends along -x from its N1 at the origin,
    # unit 2 along +y from its C1, which sits r away down +x
    unit1 = [('N1', 'n1', (0.0, 0.0, 0.0)), ('C1', 'c1', (-0.116, 0.0, 0.0)),
             ('O', 'os', (-0.246, 0.0, 0.0)), ('C2', 'c3', (-0.390, 0.0, 0.0)),
             ('H1', 'hc', (-0.43, 0.10, 0.0)), ('H2', 'hc', (-0.43, -0.05, 0.09)),
             ('H3', 'hc', (-0.43, -0.05, -0.09))]
    unit2 = [('N1', 'n1', (r + 0.116, 0.0, 0.0)), ('C1', 'c1', (r, 0.0, 0.0)),
             ('O', 'os', (r, 0.130, 0.0)), ('C2', 'c3', (r, 0.274, 0.0)),
             ('H1', 'hc', (r + 0.10, 0.310, 0.0)), ('H2', 'hc', (r - 0.05, 0.310, 0.09)),
             ('H3', 'hc', (r - 0.05, 0.310, -0.09))]
    rows, pos, nr = [], [], 1
    for resnr, unit in enumerate((unit1, unit2), start=1):
        for name, typ, (x, y, z) in unit:
            rows.append({'nr': nr, 'atom': name, 'type': typ, 'resnr': resnr, 'charge': 0.0})
            pos.append({'globalIdx': nr, 'posX': x, 'posY': y, 'posZ': z})
            nr += 1
    bonds = []
    def bond(ai, aj, order):
        bonds.append({'bondIdx': len(bonds), 'ai': ai, 'aj': aj, 'order': order})
    for unit in range(2):
        b = 7 * unit
        bond(b + 1, b + 2, cn_order)      # N1#C1
        bond(b + 2, b + 3, '1')           # C1-O
        bond(b + 3, b + 4, '1')           # O-CH3
        for h in (5, 6, 7):
            bond(b + 4, b + h, '1')
    bond(1, 9, ring_order)                # unit 1's N to unit 2's C
    T = Topology()
    T.D = {'atoms': pd.DataFrame(rows), 'mol2_bonds': pd.DataFrame(bonds)}
    return types.SimpleNamespace(Topology=T,
                                 Coordinates=types.SimpleNamespace(A=pd.DataFrame(pos)))


class TestElementOf(unittest.TestCase):
    def test_gaff_types_lead_with_their_element(self):
        self.assertEqual(_element_of('C1', 'ca'), 'C')
        self.assertEqual(_element_of('N1', 'nb'), 'N')
        self.assertEqual(_element_of('O', 'os'), 'O')

    def test_two_letter_halogens_are_not_carbon_or_boron(self):
        self.assertEqual(_element_of('CL1', 'cl'), 'CL')
        self.assertEqual(_element_of('BR1', 'br'), 'BR')

    def test_it_falls_back_to_the_name(self):
        self.assertEqual(_element_of('N1', ''), 'N')


class TestValenceRepair(unittest.TestCase):
    """An addition removes nothing, so an atom that already carried a multiple bond ends
    up over-valent unless that bond is reduced."""

    def test_an_addition_leaves_a_pentavalent_carbon(self):
        TC = cyanate_pair()
        complaints = geometry_complaints(TC)
        self.assertTrue(any('more than C can carry' in c for c in complaints))

    def test_repair_reduces_each_triple_bond_that_has_to_give(self):
        # the ring bond makes unit 1's N and unit 2's C both over-valent, and it is each
        # one's OWN C#N that has to drop -- two reductions, not one
        TC = cyanate_pair()
        self.assertEqual(TC.Topology.rebalance_mol2_bond_orders(), 2)
        mb = TC.Topology.D['mol2_bonds'].set_index(['ai', 'aj'])
        self.assertEqual(mb.loc[(1, 2), 'order'], '2')
        self.assertEqual(mb.loc[(8, 9), 'order'], '2')

    def test_repair_leaves_valence_satisfied(self):
        TC = cyanate_pair()
        TC.Topology.rebalance_mol2_bond_orders()
        self.assertEqual([c for c in geometry_complaints(TC) if 'carry' in c], [])

    def test_a_sound_molecule_is_not_touched(self):
        TC = cyanate_pair(cn_order='2')      # already reduced
        self.assertEqual(TC.Topology.rebalance_mol2_bond_orders(), 0)

    def test_it_can_be_restricted_to_the_atoms_just_bonded(self):
        TC = cyanate_pair()
        # atom 3 is not over-valent, so restricting to it repairs nothing
        self.assertEqual(TC.Topology.rebalance_mol2_bond_orders(atoms=[3]), 0)


class TestGeometryComplaints(unittest.TestCase):
    def test_a_bond_far_too_long_for_its_order(self):
        TC = cyanate_pair(ring_bond_nm=0.223)
        self.assertTrue(any('2.23 A long' in c for c in geometry_complaints(TC)))

    def test_a_closed_ring_draws_no_complaint(self):
        TC = cyanate_pair(ring_bond_nm=0.134, cn_order='2')
        self.assertEqual(geometry_complaints(TC), [])

    def test_nothing_to_check_without_mol2_bonds(self):
        TC = cyanate_pair()
        del TC.Topology.D['mol2_bonds']
        self.assertEqual(geometry_complaints(TC), [])

    def test_a_system_topocoord_is_simply_skipped(self):
        # no Coordinates at all, as a freshly constructed TopoCoord has
        TC = types.SimpleNamespace(Topology=Topology(), Coordinates=None)
        TC.Topology.D = {}
        self.assertEqual(geometry_complaints(TC), [])


class TestCheckRingClosed(unittest.TestCase):
    def test_it_refuses_an_open_ring_and_names_the_bond(self):
        TC = cyanate_pair(ring_bond_nm=0.223)
        with pytest.raises(RuntimeError, match='N1'):
            check_ring_closed(TC, 'MCY3')

    def test_it_says_charges_come_from_geometry(self):
        TC = cyanate_pair(ring_bond_nm=0.223)
        with pytest.raises(RuntimeError, match='geometry, not from the bond orders'):
            check_ring_closed(TC, 'MCY3')

    def test_a_closed_ring_passes(self):
        TC = cyanate_pair(ring_bond_nm=0.134, cn_order='2')
        self.assertIsNone(check_ring_closed(TC, 'MCY3'))


class TestRelaxGeometry(unittest.TestCase):
    """The repair itself: a bond recorded but never closed is pulled to a real length."""

    def test_it_closes_a_bond_that_was_only_recorded(self):
        pytest.importorskip('rdkit')
        TC = cyanate_pair(ring_bond_nm=0.223)
        TC.Topology.rebalance_mol2_bond_orders()
        self.assertTrue(relax_geometry(TC, name='pair'))
        A = TC.Coordinates.A.set_index('globalIdx')
        import numpy as np
        d = float(np.linalg.norm(A.loc[1, ['posX', 'posY', 'posZ']].to_numpy(float)
                                 - A.loc[9, ['posX', 'posY', 'posZ']].to_numpy(float)))
        self.assertLess(d, 0.18)

    def test_an_unreadable_template_is_reported_not_raised(self):
        pytest.importorskip('rdkit')
        # left pentavalent, so RDKit refuses to sanitize it; a build that cannot relax
        # its template should say so and let the guard stop it, not crash here
        TC = cyanate_pair(ring_bond_nm=0.223)
        with self.assertLogs('htpolynet.geometry.templategeom', level='WARNING'):
            self.assertFalse(relax_geometry(TC, name='pair'))

    def test_it_falls_back_to_embedding_when_minimization_will_not_start(self):
        # a bond recorded across 2.2 A can leave the line search unable to take its
        # first step; giving up there would strand the build at the guard
        pytest.importorskip('rdkit')
        TC = cyanate_pair(ring_bond_nm=0.223)
        TC.Topology.rebalance_mol2_bond_orders()
        self.assertTrue(relax_geometry(TC, name='pair'))
        self.assertEqual(geometry_complaints(TC), [])
