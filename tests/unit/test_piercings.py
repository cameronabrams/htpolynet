"""

.. module:: test_piercings
   :synopsis: finding a bond threaded through a ring in a finished network

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import types
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.piercings import (detect_rings, find_piercings, format_report,
                                          ring_composition)

BOX = np.array([2.0, 2.0, 2.0])


def hexagon(centre, radius=0.14, normal='z'):
    """Six points on a ring of `radius` about `centre`, wrapped into the box."""
    a = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    if normal == 'z':
        pts = np.stack([radius * np.cos(a), radius * np.sin(a), np.zeros(6)], axis=1)
    else:
        pts = np.stack([np.zeros(6), radius * np.cos(a), radius * np.sin(a)], axis=1)
    return [np.mod(np.asarray(centre, dtype=float) + p, BOX) for p in pts]


class System:
    """Builds a TopoCoord stub from rings and bonds placed by hand."""

    def __init__(self):
        self.pos, self.elem, self.res, self.bonds, self.rings = {}, {}, {}, [], []
        self._n = 0

    def add_ring(self, centre, elements='CCCCCC', resnr=None, normal='z', radius=0.14):
        idx = []
        for k, p in enumerate(hexagon(centre, radius=radius, normal=normal)):
            self._n += 1
            self.pos[self._n] = p
            self.elem[self._n] = elements[k]
            # a triazine spans three residues, a phenyl one: classification must not
            # depend on that
            self.res[self._n] = resnr if resnr is not None else 900 + (k % 3)
            idx.append(self._n)
        for k in range(6):
            self.bonds.append((idx[k], idx[(k + 1) % 6]))
        self.rings.append(idx)
        return idx

    def add_bond(self, p0, p1, resnr=800, element='C'):
        self._n += 1
        a = self._n
        self._n += 1
        b = self._n
        self.pos[a], self.pos[b] = np.mod(p0, BOX), np.mod(p1, BOX)
        self.elem[a] = self.elem[b] = element
        self.res[a] = self.res[b] = resnr
        self.bonds.append((a, b))
        return a, b

    def TC(self, with_ring_list=True):
        ids = sorted(self.pos)
        A = pd.DataFrame([{'globalIdx': i, 'posX': self.pos[i][0], 'posY': self.pos[i][1],
                           'posZ': self.pos[i][2]} for i in ids])
        at = pd.DataFrame([{'nr': i, 'atom': f'{self.elem[i]}{i}', 'type': self.elem[i].lower(),
                            'resnr': self.res[i]} for i in ids])
        bd = pd.DataFrame(self.bonds, columns=['ai', 'aj'])
        rings = [types.SimpleNamespace(idx=r) for r in self.rings] if with_ring_list else []
        return types.SimpleNamespace(
            Coordinates=types.SimpleNamespace(A=A, box=np.diag(BOX)),
            Topology=types.SimpleNamespace(D={'atoms': at, 'bonds': bd}, rings=rings))


class TestAnUncuredMeltHasNone(unittest.TestCase):
    """The case that caught the worst bug: a detector taking the centroid of RAW
    coordinates puts the centre of any boundary-straddling ring in the middle of the box
    and its triangle fan across the whole cell, so it reports a piercing for nearly
    everything.  One early version reported 212 on a melt where the true answer is zero
    by construction.  About 15% of rings straddle a boundary in a real system."""

    def test_rings_straddling_the_boundary_do_not_catch_distant_bonds(self):
        s = System()
        s.add_ring((0.02, 1.0, 1.0))      # straddles x = 0
        s.add_ring((1.0, 0.01, 1.0))      # straddles y = 0
        s.add_ring((1.99, 1.0, 0.03))     # straddles x = 2 and z = 0
        for c in ((1.0, 1.0, 1.0), (0.5, 1.5, 0.5), (1.5, 0.5, 1.5)):
            s.add_bond(np.array(c) - np.array([0, 0, 0.08]),
                       np.array(c) + np.array([0, 0, 0.08]))
        self.assertEqual(find_piercings(s.TC()), [])

    def test_a_melt_of_separated_rings_is_clean(self):
        s = System()
        for x in (0.3, 0.9, 1.5):
            s.add_ring((x, 1.0, 1.0))
        self.assertEqual(find_piercings(s.TC()), [])


class TestABoundaryStraddlingRingStillReportsARealPiercing(unittest.TestCase):
    """The same trap in the other direction: imaging must not make a true hit vanish."""

    def test_it_is_found(self):
        s = System()
        s.add_ring((0.02, 1.0, 1.0))
        s.add_bond(np.array([0.02, 1.0, 0.92]), np.array([0.02, 1.0, 1.08]))
        found = find_piercings(s.TC())
        self.assertEqual(len(found), 1)
        self.assertLess(found[0]['offset'], 0.01)
        self.assertLess(found[0]['angle'], 1.0)

    def test_it_is_found_when_the_ring_wraps_two_axes(self):
        s = System()
        s.add_ring((1.99, 0.01, 1.0))
        s.add_bond(np.array([1.99, 0.01, 0.92]), np.array([1.99, 0.01, 1.08]))
        self.assertEqual(len(find_piercings(s.TC())), 1)


class TestAnOverlongPiercingBond(unittest.TestCase):
    """2.4 to 3.2 A piercing bonds are the normal case in real data, not outliers;
    anything keyed to an ordinary bond length misses them."""

    def test_a_32_angstrom_bond_is_found_and_reported_as_long(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0))
        s.add_bond(np.array([1.0, 1.0, 0.84]), np.array([1.0, 1.0, 1.16]))
        found = find_piercings(s.TC())
        self.assertEqual(len(found), 1)
        self.assertAlmostEqual(found[0]['length'], 0.32, places=6)
        self.assertIn('1 of 1 piercing bond(s) are longer than', format_report(found))

    def test_a_24_angstrom_bond_is_found(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0))
        s.add_bond(np.array([1.0, 1.0, 0.88]), np.array([1.0, 1.0, 1.12]))
        self.assertEqual(len(find_piercings(s.TC())), 1)


class TestCompositionIsClassifiedOnTheWholeRing(unittest.TestCase):
    """Two ways this has been got wrong: keying on how many residues a ring spans (a
    triazine spans exactly three, a phenyl one), and reading a printout that shows only
    the first three atoms."""

    def test_a_triazine_and_a_phenyl_are_distinguished(self):
        self.assertEqual(ring_composition([1, 2, 3, 4, 5, 6],
                                          {1: 'C', 2: 'N', 3: 'C', 4: 'N', 5: 'C', 6: 'N'}),
                         'C3N3')
        self.assertEqual(ring_composition([1, 2, 3, 4, 5, 6], {i: 'C' for i in range(1, 7)}),
                         'C6')

    def test_the_first_three_atoms_do_not_decide_it(self):
        # C,C,C then N,N,N -- a detector reading ring[:3] would call this C6
        self.assertEqual(ring_composition([1, 2, 3, 4, 5, 6],
                                          {1: 'C', 2: 'C', 3: 'C', 4: 'N', 5: 'N', 6: 'N'}),
                         'C3N3')

    def test_residue_span_does_not_decide_it(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0), elements='CNCNCN', resnr=5)   # one residue, C3N3
        s.add_bond(np.array([1.0, 1.0, 0.92]), np.array([1.0, 1.0, 1.08]))
        found = find_piercings(s.TC())
        self.assertEqual(found[0]['composition'], 'C3N3')


class TestSelfVersusForeignAndInterlocking(unittest.TestCase):
    def test_a_bond_of_the_ring_s_own_residue_is_marked_self(self):
        s = System()
        ring = s.add_ring((1.0, 1.0, 1.0), resnr=42)
        s.add_bond(np.array([1.0, 1.0, 0.92]), np.array([1.0, 1.0, 1.08]), resnr=42)
        found = find_piercings(s.TC())
        self.assertTrue(found[0]['self'])

    def test_a_foreign_bond_is_not(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0), resnr=42)
        s.add_bond(np.array([1.0, 1.0, 0.92]), np.array([1.0, 1.0, 1.08]), resnr=99)
        self.assertFalse(find_piercings(s.TC())[0]['self'])

    def test_two_interlocked_rings_are_both_reported(self):
        # catenane-like: each ring threaded by a bond, as seen twice in real builds
        s = System()
        s.add_ring((1.0, 1.0, 1.0), normal='z')
        s.add_bond(np.array([1.0, 1.0, 0.92]), np.array([1.0, 1.0, 1.08]))
        s.add_ring((0.5, 0.5, 0.5), normal='x')
        s.add_bond(np.array([0.42, 0.5, 0.5]), np.array([0.58, 0.5, 0.5]))
        self.assertEqual(len(find_piercings(s.TC())), 2)


class TestRingDetectionFallback(unittest.TestCase):
    """A bare top/gro pair has no ring list; enumerating all chordless cycles does not
    finish on a percolated network, so the smallest ring through each bond is taken."""

    def test_it_finds_a_six_ring(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0))
        rings = detect_rings(s.TC().Topology.D['bonds'])
        self.assertEqual(len(rings), 1)
        self.assertEqual(len(rings[0]), 6)

    def test_a_chain_has_no_rings(self):
        bonds = pd.DataFrame([{'ai': i, 'aj': i + 1} for i in range(1, 8)])
        self.assertEqual(detect_rings(bonds), [])

    def test_a_piercing_is_found_without_a_ring_list(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0))
        s.add_bond(np.array([1.0, 1.0, 0.92]), np.array([1.0, 1.0, 1.08]))
        self.assertEqual(len(find_piercings(s.TC(with_ring_list=False))), 1)


class TestReport(unittest.TestCase):
    def test_a_clean_system_says_so(self):
        self.assertIn('No bond passes through a ring', format_report([]))

    def test_the_summary_counts_by_composition(self):
        s = System()
        s.add_ring((1.0, 1.0, 1.0), elements='CNCNCN')
        s.add_bond(np.array([1.0, 1.0, 0.92]), np.array([1.0, 1.0, 1.08]))
        self.assertIn('C3N3 1', format_report(find_piercings(s.TC())))
