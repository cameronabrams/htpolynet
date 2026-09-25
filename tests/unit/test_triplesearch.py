"""

.. module:: test_triplesearch
   :synopsis: finding the triples of cyanate groups one cyclotrimerization can join

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import numpy as np
import pandas as pd
import pytest

from htpolynet.cure.triplesearch import (bonds_by_atom, candidate_triples,
                                         reactive_sites, residue_map,
                                         ring_threading_bonds, select_disjoint,
                                         triple_bonds_dataframe)

BOX = np.array([10.0, 10.0, 10.0])


def system(groups, resname='BCY', z=1):
    """Builds an atom table of cyanate groups.

    groups: list of (resNum, molecule, donor position, acceptor position)
    """
    rows, idx = [], 1
    positions = {}
    for resnum, mol, dpos, apos in groups:
        for name, pos in (('N1', dpos), ('C1', apos)):
            rows.append({'globalIdx': idx, 'resNum': resnum, 'resName': resname,
                         'atomName': name, 'molecule': mol, 'z': z})
            positions[idx] = np.array(pos, dtype=float)
            idx += 1
    return pd.DataFrame(rows), positions


def triangle(side=0.3, resnums=(1, 2, 3), molecules=(1, 2, 3)):
    """Three groups whose donor-acceptor distances are all about `side`."""
    groups = []
    for n, (r, m) in enumerate(zip(resnums, molecules)):
        theta = 2.0 * np.pi * n / 3.0
        centre = np.array([5.0 + side * np.cos(theta), 5.0 + side * np.sin(theta), 5.0])
        # donor points at the next group's acceptor position, roughly
        groups.append((r, m, centre + np.array([0.0, 0.0, 0.02]), centre))
    return system(groups)


def dicyanate_plus_two():
    """One residue carrying two cyanate groups, plus two single-group residues, all
    close enough to ring: the only available triangle needs both arms of the dicyanate."""
    rows, positions, idx = [], {}, 1
    spec = [(1, 1, 'N1', [5.00, 5.00, 5.02]), (1, 1, 'C1', [5.00, 5.00, 5.00]),
            (1, 1, 'N2', [5.25, 5.00, 5.02]), (1, 1, 'C2', [5.25, 5.00, 5.00]),
            (2, 2, 'N1', [5.12, 5.22, 5.02]), (2, 2, 'C1', [5.12, 5.22, 5.00])]
    for resnum, mol, name, pos in spec:
        rows.append({'globalIdx': idx, 'resNum': resnum, 'resName': 'BCY',
                     'atomName': name, 'molecule': mol, 'z': 1})
        positions[idx] = np.array(pos, dtype=float)
        idx += 1
    return pd.DataFrame(rows), positions


class TestReactiveSites(unittest.TestCase):
    def test_one_row_per_group(self):
        adf, _ = triangle()
        s = reactive_sites(adf, 'BCY')
        self.assertEqual(s.shape[0], 3)
        self.assertEqual(sorted(s['resNum']), [1, 2, 3])
        self.assertEqual(sorted(s['donor']), [1, 3, 5])
        self.assertEqual(sorted(s['acceptor']), [2, 4, 6])

    def test_a_reacted_group_is_skipped(self):
        adf, _ = triangle()
        adf.loc[adf['resNum'] == 2, 'z'] = 0
        self.assertEqual(sorted(reactive_sites(adf, 'BCY')['resNum']), [1, 3])

    def test_z_can_be_ignored(self):
        adf, _ = triangle()
        adf['z'] = 0
        self.assertEqual(reactive_sites(adf, 'BCY', require_z=False).shape[0], 3)

    def test_another_residue_is_not_a_group(self):
        adf, _ = triangle()
        other = pd.DataFrame([{'globalIdx': 99, 'resNum': 9, 'resName': 'XXX',
                               'atomName': 'N1', 'molecule': 9, 'z': 1}])
        self.assertEqual(reactive_sites(pd.concat([adf, other]), 'BCY').shape[0], 3)


class TestCandidateTriples(unittest.TestCase):
    def test_a_close_triangle_is_found(self):
        adf, pos = triangle(side=0.25)
        cands = candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.6, BOX)
        self.assertEqual(len(cands), 1)
        score, ring, bonds = cands[0]
        self.assertEqual(len(bonds), 3)
        self.assertEqual(sorted(ring), [0, 1, 2])
        # each group donates once and accepts once
        self.assertEqual(sorted(a for a, b in bonds), [0, 1, 2])
        self.assertEqual(sorted(b for a, b in bonds), [0, 1, 2])

    def test_a_spread_out_triangle_is_not(self):
        adf, pos = triangle(side=2.0)
        self.assertEqual(candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.5, BOX), [])

    def test_fewer_than_three_groups_gives_nothing(self):
        adf, pos = system([(1, 1, [5, 5, 5], [5, 5, 5.02]), (2, 2, [5.2, 5, 5], [5.2, 5, 5.02])])
        self.assertEqual(candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.6, BOX), [])

    def test_two_groups_of_one_molecule_are_refused(self):
        adf, pos = triangle(side=0.25, resnums=(1, 2, 3), molecules=(1, 1, 3))
        self.assertEqual(candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.6, BOX), [])

    def test_same_molecule_allowed_when_asked(self):
        adf, pos = triangle(side=0.25, resnums=(1, 2, 3), molecules=(1, 1, 3))
        cands = candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.6, BOX, same_molecule=True)
        self.assertEqual(len(cands), 1)

    def test_two_groups_of_one_residue_are_refused(self):
        # both arms of one dicyanate: a 14-ring, and routing that pair through the
        # pairwise bondtest would hit its assert rather than be declined
        adf, pos = dicyanate_plus_two()
        sites = reactive_sites(adf, 'BCY', groups=(('N1', 'C1'), ('N2', 'C2')))
        self.assertEqual(sites.shape[0], 3)   # two arms of residue 1, one of residue 2
        self.assertEqual(candidate_triples(sites, pos, 0.9, BOX), [])

    def test_and_allowed_when_asked(self):
        # the two arms are also one molecule, so both doors have to be opened
        adf, pos = dicyanate_plus_two()
        sites = reactive_sites(adf, 'BCY', groups=(('N1', 'C1'), ('N2', 'C2')))
        self.assertTrue(candidate_triples(sites, pos, 0.9, BOX,
                                          same_residue=True, same_molecule=True))

    def test_the_score_is_the_total_bond_length(self):
        adf, pos = triangle(side=0.25)
        sites = reactive_sites(adf, 'BCY')
        score, ring, bonds = candidate_triples(sites, pos, 0.6, BOX)[0]
        expected = sum(float(np.linalg.norm(pos[int(sites.at[a, 'donor'])]
                                            - pos[int(sites.at[b, 'acceptor'])]))
                       for a, b in bonds)
        self.assertAlmostEqual(score, expected, places=6)

    def test_candidates_come_back_shortest_first(self):
        # a tight triangle and a loose one, sharing no groups
        tight, pos_t = triangle(side=0.2, resnums=(1, 2, 3), molecules=(1, 2, 3))
        loose, pos_l = triangle(side=0.45, resnums=(4, 5, 6), molecules=(4, 5, 6))
        loose = loose.copy()
        loose['globalIdx'] = loose['globalIdx'] + 100
        pos_l = {k + 100: v + np.array([2.0, 0.0, 0.0]) for k, v in pos_l.items()}
        adf = pd.concat([tight, loose], ignore_index=True)
        pos = {**pos_t, **pos_l}
        cands = candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.9, BOX)
        self.assertEqual(len(cands), 2)
        self.assertLess(cands[0][0], cands[1][0])

    def test_it_reaches_across_the_periodic_boundary(self):
        # a triangle straddling x = 0
        groups = [(1, 1, [9.95, 5, 5.02], [9.95, 5, 5]),
                  (2, 2, [0.1, 5, 5.02], [0.1, 5, 5]),
                  (3, 3, [0.0, 5.2, 5.02], [0.0, 5.2, 5])]
        adf, pos = system(groups)
        self.assertEqual(len(candidate_triples(reactive_sites(adf, 'BCY'), pos, 0.6, BOX)), 1)


class TestSelectDisjoint(unittest.TestCase):
    def candidates(self):
        return [(0.6, (0, 1, 2), [(0, 1), (1, 2), (2, 0)]),
                (0.7, (2, 3, 4), [(2, 3), (3, 4), (4, 2)]),
                (0.8, (3, 4, 5), [(3, 4), (4, 5), (5, 3)])]

    def test_no_group_is_used_twice(self):
        chosen = select_disjoint(self.candidates())
        self.assertEqual([c[1] for c in chosen], [(0, 1, 2), (3, 4, 5)])

    def test_shortest_first(self):
        chosen = select_disjoint(self.candidates())
        self.assertEqual(chosen[0][0], 0.6)

    def test_a_limit_is_honored(self):
        self.assertEqual(len(select_disjoint(self.candidates(), max_triples=1)), 1)

    def test_nothing_to_choose(self):
        self.assertEqual(select_disjoint([]), [])


class TestBondTableAndResidueMap(unittest.TestCase):
    def setUp(self):
        adf, pos = triangle(side=0.25)
        self.sites = reactive_sites(adf, 'BCY')
        self.chosen = select_disjoint(candidate_triples(self.sites, pos, 0.6, BOX))

    def test_three_bonds_per_triple_with_the_expected_columns(self):
        bdf = triple_bonds_dataframe(self.chosen, self.sites, 'BCY3')
        self.assertEqual(bdf.shape[0], 3)
        self.assertEqual(list(bdf.columns), ['ai', 'aj', 'ri', 'rj', 'order', 'reactantName', 'triple'])
        self.assertTrue((bdf['reactantName'] == 'BCY3').all())
        self.assertTrue((bdf['triple'] == 0).all())
        # a donor bonds an acceptor, never two donors
        self.assertTrue(set(bdf['ai']).isdisjoint(set(bdf['aj'])))

    def test_the_residue_map_pairs_template_and_instance(self):
        maps = residue_map(self.chosen, self.sites, [1, 2, 3])
        self.assertEqual(len(maps), 1)
        self.assertEqual(sorted(maps[0]), [1, 2, 3])
        self.assertEqual(sorted(maps[0].values()), [1, 2, 3])

    def test_a_template_of_the_wrong_size_is_refused(self):
        with pytest.raises(ValueError, match='residue'):
            residue_map(self.chosen, self.sites, [1, 2])


class TestThreadingRejection(unittest.TestCase):
    """A ring that closes around an existing bond threads that monomer through it for
    good: 15 stretched bonds over ten builds, 13 of them threading, none threading
    without being stretched.  10 of 11 were pre-existing monomer backbones, so the
    triazine closes around a monomer that was already there."""

    def hexagon(self, radius=1.0, z=0.0):
        a = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        return {i + 1: np.array([radius * np.cos(t), radius * np.sin(t), z])
                for i, t in enumerate(a)}

    def bonds(self, pairs):
        return bonds_by_atom(pd.DataFrame(pairs, columns=['ai', 'aj']))

    def test_a_bond_through_the_middle_is_found(self):
        pos = self.hexagon()
        pos[7], pos[8] = np.array([0., 0., -0.5]), np.array([0., 0., 0.5])
        found = ring_threading_bonds([1, 2, 3, 4, 5, 6], pos, self.bonds([(7, 8)]), BOX)
        self.assertEqual(found, [(7, 8)])

    def test_a_bond_beside_the_ring_is_not(self):
        pos = self.hexagon()
        pos[7], pos[8] = np.array([0.9, 0.9, -0.5]), np.array([0.9, 0.9, 0.5])
        self.assertEqual(ring_threading_bonds([1, 2, 3, 4, 5, 6], pos,
                                              self.bonds([(7, 8)]), BOX), [])

    def test_a_bond_that_stops_short_does_not_thread(self):
        # both endpoints on the same side: a segment, not an infinite ray
        pos = self.hexagon()
        pos[7], pos[8] = np.array([0., 0., 0.5]), np.array([0., 0., 1.5])
        self.assertEqual(ring_threading_bonds([1, 2, 3, 4, 5, 6], pos,
                                              self.bonds([(7, 8)]), BOX), [])

    def test_the_ring_s_own_bonds_are_not_counted(self):
        pos = self.hexagon()
        self.assertEqual(ring_threading_bonds([1, 2, 3, 4, 5, 6], pos,
                                              self.bonds([(1, 2), (3, 4)]), BOX), [])

    def test_it_works_across_a_periodic_boundary(self):
        # the same threaded arrangement, translated so the ring straddles the edge; a
        # centroid of raw wrapped coordinates is what gives false positives here
        pos = self.hexagon()
        shift = np.array([BOX[0], 0.0, 0.0])
        pos = {k: np.mod(v + shift * 0.5 + np.array([BOX[0] / 2, 0, 0]), BOX)
               for k, v in pos.items()}
        centre = np.array([0., 0., 0.]) + shift * 0.5 + np.array([BOX[0] / 2, 0, 0])
        pos[7] = np.mod(centre + np.array([0., 0., -0.5]), BOX)
        pos[8] = np.mod(centre + np.array([0., 0., 0.5]), BOX)
        found = ring_threading_bonds([1, 2, 3, 4, 5, 6], pos, self.bonds([(7, 8)]), BOX)
        self.assertEqual(found, [(7, 8)])

    def test_a_non_planar_loop_still_works(self):
        # at candidate time the loop is three groups a search radius apart, not a ring
        pos = self.hexagon(radius=2.0)
        for k in (2, 4, 6):
            pos[k] = pos[k] + np.array([0., 0., 0.8])
        pos[7], pos[8] = np.array([0., 0., -1.0]), np.array([0., 0., 1.4])
        found = ring_threading_bonds([1, 2, 3, 4, 5, 6], pos, self.bonds([(7, 8)]), BOX)
        self.assertEqual(found, [(7, 8)])
