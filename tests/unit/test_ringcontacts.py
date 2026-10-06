"""

.. module:: test_ringcontacts
   :synopsis: face-on contacts, ring-centroid RDFs and ring-flip clearance on hand-placed systems

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import math
import types
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.ringcontacts import (RingSystem, bridge_groups, centroid_rdf,
                                             face_on_contacts, flip_axis, ring_flip_clearance)

BOX = np.array([4.0, 4.0, 4.0])
MASS = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'H': 1.008}


class Builder:
    """A TopoCoord stub from atoms placed by hand."""

    def __init__(self):
        self.atoms, self.bonds = {}, []

    def atom(self, pos, element, resnr, typ=None):
        i = len(self.atoms) + 1
        self.atoms[i] = (np.mod(np.asarray(pos, dtype=float), BOX), element, resnr,
                         typ or element.lower())
        return i

    def bond(self, a, b):
        self.bonds.append((a, b))

    def ring(self, centre, elements='CCCCCC', resnr=1, radius=0.14):
        """A hexagon in the xy plane, so its normal is z.  Atom k sits at angle 60k deg,
        so atoms 0 and 3 are a para pair along x."""
        idx = []
        for k in range(6):
            a = 2 * math.pi * k / 6
            idx.append(self.atom(np.asarray(centre) + [radius * math.cos(a), radius * math.sin(a), 0.0],
                                 elements[k], resnr))
        for k in range(6):
            self.bond(idx[k], idx[(k + 1) % 6])
        return idx

    def TC(self):
        ids = sorted(self.atoms)
        A = pd.DataFrame([{'globalIdx': i, 'posX': self.atoms[i][0][0], 'posY': self.atoms[i][0][1],
                           'posZ': self.atoms[i][0][2]} for i in ids])
        at = pd.DataFrame([{'nr': i, 'atom': f'{self.atoms[i][1]}{i}', 'type': self.atoms[i][3],
                            'resnr': self.atoms[i][2], 'mass': MASS[self.atoms[i][1]]} for i in ids])
        bd = pd.DataFrame(self.bonds, columns=['ai', 'aj'])
        return types.SimpleNamespace(
            Coordinates=types.SimpleNamespace(A=A, box=np.diag(BOX)),
            Topology=types.SimpleNamespace(D={'atoms': at, 'bonds': bd}, rings=[]))


def system(b):
    return RingSystem(b.TC())


class TestRings(unittest.TestCase):
    def test_compositions_are_found_by_formula(self):
        b = Builder()
        b.ring((1.0, 1.0, 1.0), 'CNCNCN')
        b.ring((3.0, 3.0, 3.0))
        S = system(b)
        self.assertEqual(len(S.rings_of('C3N3')), 1)
        self.assertEqual(len(S.rings_of('C6')), 1)


class TestFaceOnContacts(unittest.TestCase):
    def test_a_donor_over_the_face_counts_and_one_in_the_plane_does_not(self):
        b = Builder()
        b.ring((2.0, 2.0, 2.0), 'CNCNCN', resnr=1)
        over = b.atom((2.0, 2.0, 2.33), 'O', 2)
        edge = b.atom((2.33, 2.0, 2.0), 'O', 3)
        S = system(b)
        rings = S.rings_of('C3N3')
        self.assertEqual(face_on_contacts(S, [over], rings)['observed'], 1)
        self.assertEqual(face_on_contacts(S, [edge], rings)['observed'], 0)

    def test_the_face_is_found_across_a_periodic_boundary(self):
        b = Builder()
        b.ring((2.0, 2.0, 3.95), 'CNCNCN', resnr=1)
        o = b.atom((2.0, 2.0, 0.28), 'O', 2)      # 0.33 nm above, through z = 4
        S = system(b)
        self.assertEqual(face_on_contacts(S, [o], S.rings_of('C3N3'))['observed'], 1)

    def test_expectation_is_pairs_times_the_cone_fraction_of_the_box(self):
        b = Builder()
        b.ring((2.0, 2.0, 2.0), 'CNCNCN', resnr=1)
        o = b.atom((0.5, 0.5, 0.5), 'O', 2)
        S = system(b)
        r = face_on_contacts(S, [o], S.rings_of('C3N3'), rmax=0.37, max_angle=35.0)
        vcone = (1 - math.cos(math.radians(35.0))) * 4 / 3 * math.pi * 0.37 ** 3
        self.assertEqual(r['pairs'], 1)
        self.assertAlmostEqual(r['expected'], vcone / 64.0)

    def test_a_shell_excludes_a_contact_inside_rmin_and_shrinks_the_expectation(self):
        b = Builder()
        b.ring((2.0, 2.0, 2.0), 'CNCNCN', resnr=1)
        o = b.atom((2.0, 2.0, 2.25), 'O', 2)
        S = system(b)
        rings = S.rings_of('C3N3')
        full = face_on_contacts(S, [o], rings, rmax=0.36)
        shell = face_on_contacts(S, [o], rings, rmin=0.30, rmax=0.36)
        self.assertEqual(full['observed'], 1)
        self.assertEqual(shell['observed'], 0)
        self.assertAlmostEqual(shell['expected'] / full['expected'],
                               (0.36 ** 3 - 0.30 ** 3) / 0.36 ** 3)

    def test_a_ring_the_donor_is_bonded_near_is_not_eligible(self):
        """Without the exclusion a donor scores contacts with its own substituents."""
        b = Builder()
        ring = b.ring((2.0, 2.0, 2.0), 'CNCNCN', resnr=1)
        o = b.atom((2.0, 2.0, 2.33), 'O', 1)
        b.bond(o, ring[0])
        S = system(b)
        rings = S.rings_of('C3N3')
        r = face_on_contacts(S, [o], rings)
        self.assertEqual((r['pairs'], r['observed']), (0, 0))
        # a 0-bond exclusion keeps it, so the exclusion is what removed it
        self.assertEqual(face_on_contacts(S, [o], rings, exclude_bonds=0)['observed'], 1)

    def test_residue_exclusion_drops_the_donors_own_residue_only(self):
        b = Builder()
        b.ring((2.0, 2.0, 2.0), 'CNCNCN', resnr=1)
        mine = b.atom((2.0, 2.0, 2.33), 'O', 1)
        theirs = b.atom((2.0, 2.0, 1.67), 'O', 2)
        S = system(b)
        rings = S.rings_of('C3N3')
        self.assertEqual(face_on_contacts(S, [mine], rings, exclude='residue')['pairs'], 0)
        self.assertEqual(face_on_contacts(S, [theirs], rings, exclude='residue')['observed'], 1)

    def test_an_unknown_exclusion_is_refused(self):
        b = Builder()
        b.ring((2.0, 2.0, 2.0), 'CNCNCN')
        o = b.atom((0.5, 0.5, 0.5), 'O', 2)
        S = system(b)
        with self.assertRaises(ValueError):
            face_on_contacts(S, [o], S.rings_of('C3N3'), exclude='chain')


def bridged_pair(b, centre_el='S', subst=('O', 'O'), resnr=1):
    """Two phenylenes joined at their para atoms by a bridge centre with substituents."""
    r1 = b.ring((1.0, 2.0, 2.0), resnr=resnr)
    r2 = b.ring((1.6, 2.0, 2.0), resnr=resnr)
    c = b.atom((1.3, 2.0, 2.0), centre_el, resnr)
    b.bond(c, r1[0])
    b.bond(c, r2[3])
    subs = []
    for k, el in enumerate(subst):
        s = b.atom((1.3, 2.0 + (0.14 if k == 0 else -0.14), 2.0), el, resnr)
        b.bond(c, s)
        subs.append(s)
    return r1, r2, c, subs


class TestBridgeGroups(unittest.TestCase):
    def test_a_sulfone_bridge_is_its_centre_and_both_oxygens(self):
        b = Builder()
        _, _, c, subs = bridged_pair(b)
        S = system(b)
        self.assertEqual(bridge_groups(S, 'C6'), [[c] + sorted(subs)])

    def test_an_aryl_ether_to_a_triazine_is_not_a_bridge(self):
        b = Builder()
        ph = b.ring((1.0, 2.0, 2.0))
        tz = b.ring((1.6, 2.0, 2.0), 'CNCNCN', resnr=2)
        o = b.atom((1.3, 2.0, 2.0), 'O', 1)
        b.bond(o, ph[0])
        b.bond(o, tz[3])
        S = system(b)
        self.assertEqual(bridge_groups(S, 'C6'), [])


class TestCentroidRDF(unittest.TestCase):
    def test_a_ring_lands_in_its_bin_and_the_own_residue_ring_does_not(self):
        b = Builder()
        b.ring((2.0, 2.0, 2.0), 'CNCNCN', resnr=1)     # own residue, 0.31 nm away
        b.ring((2.0, 2.0, 1.2), 'CNCNCN', resnr=5)     # foreign, 0.49 nm away
        o = b.atom((2.0, 2.0, 1.69), 'O', 1)
        S = system(b)
        x, g = centroid_rdf(S, [[o]], S.rings_of('C3N3'), rmax=1.0, dr=0.02)
        self.assertEqual(int(np.argmax(g)), int(0.49 / 0.02))
        self.assertEqual(int((g > 0).sum()), 1)


class TestRingFlip(unittest.TestCase):
    def test_the_axis_is_the_substituted_para_pair(self):
        """A hexagon has three para pairs at the same separation; only the substituted
        one is a flip axis, and picking the farthest-apart pair gets it right about a
        third of the time."""
        b = Builder()
        ring = b.ring((2.0, 2.0, 2.0))
        # substituents on atoms 1 and 4, so the first farthest-apart pair (0, 3) is wrong
        for k in (1, 4):
            p = 2.0 + 0.29 * np.array([math.cos(math.pi * k / 3), math.sin(math.pi * k / 3), 0.0])
            s = b.atom(p, 'O', 1)
            b.bond(s, ring[k])
        S = system(b)
        a, c, para = flip_axis(S, S.rings_of('C6')[0])
        self.assertTrue(para)
        self.assertEqual({a, c}, {ring[1], ring[4]})

    def test_an_ortho_methyl_widens_the_sweep_and_a_neighbour_in_it_blocks(self):
        b = Builder()
        ring = b.ring((2.0, 2.0, 2.0))
        for k, x in ((0, 1.7), (3, 2.3)):
            s = b.atom((x, 2.0, 2.0), 'O', 1)
            b.bond(s, ring[k])
        S0 = system(b)
        bare = ring_flip_clearance(S0, S0.rings_of('C6'))['R_sweep'][0]
        me = b.atom((2.07, 2.25, 2.0), 'C', 1)          # ortho methyl on ring atom 1
        b.bond(me, ring[1])
        h = b.atom((2.07, 2.36, 2.0), 'H', 1)
        b.bond(h, me)
        b.atom((2.0, 2.0, 2.30), 'C', 9)                # foreign, over the face
        S = system(b)
        res = ring_flip_clearance(S, S.rings_of('C6'))
        self.assertAlmostEqual(bare, 0.14 * math.sin(math.pi / 3), places=6)
        self.assertAlmostEqual(res['R_sweep'][0], 0.36, places=6)
        self.assertEqual(res['blocking'][0], 1)
        self.assertEqual(res['axis_fallback'], 0)


if __name__ == '__main__':
    unittest.main()
