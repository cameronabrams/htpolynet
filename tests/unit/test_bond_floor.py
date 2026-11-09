import unittest
from htpolynet.cure.curecontroller import CureController


def controller(max_bonds, **controls):
    cc = CureController({'controls': {**CureController.curedict_defaults['controls'], **controls}})
    cc.state.max_nxlinkbonds = max_bonds
    return cc


class TestBondFloor(unittest.TestCase):
    def test_the_default_is_the_fixed_count_at_any_size(self):
        # the fraction is off by default, so the enlarged examples cure as they did
        self.assertEqual(controller(400).bond_floor(10**6, 10**6), 10)
        self.assertEqual(controller(3200).bond_floor(10**6, 10**6), 10)

    def test_a_fraction_grows_the_floor_with_the_box(self):
        self.assertEqual(controller(400, min_bond_fraction_per_iteration=0.025).bond_floor(10**6, 10**6), 10)
        self.assertEqual(controller(3200, min_bond_fraction_per_iteration=0.025).bond_floor(10**6, 10**6), 80)

    def test_the_floor_is_clamped_by_what_is_left_and_what_is_accepted(self):
        cc = controller(3200, min_bond_fraction_per_iteration=0.025)
        self.assertEqual(cc.bond_floor(5, 10**6), 5)
        self.assertEqual(cc.bond_floor(10**6, 7), 7)
        self.assertEqual(cc.bond_floor(0, 0), 1)
