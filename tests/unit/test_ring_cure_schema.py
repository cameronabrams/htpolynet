"""

.. module:: test_ring_cure_schema
   :synopsis: every ring-cure knob the controller honors is one a config can actually set

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import yaml

from htpolynet.core.configuration import schema_path
from htpolynet.cure.ringcontroller import RingController


def schema_block(name):
    """The named top-level directive's attributes, as {attribute: default}."""
    with schema_path() as p:
        base = yaml.safe_load(open(p))
    for d in base['attributes']:
        if d['name'] == name:
            return {a['name']: a.get('default') for a in d.get('attributes', [])}, d
    raise AssertionError(f'no {name} directive in base.yaml')


class TestEveryControllerKnobIsReachable(unittest.TestCase):
    """A knob in RingController.defaults but not in the schema is not merely
    undocumented -- ycleptic rejects the config outright, so the knob cannot be set at
    all.  `relax` was missing for the whole life of the ring cure, and `max_accept` was
    written into CURE's `drag` block by mistake, where it was silently accepted and did
    nothing."""

    def test_the_top_level_keys_match_exactly(self):
        attrs, _ = schema_block('ring_cure')
        self.assertEqual(sorted(attrs), sorted(RingController.defaults))

    def test_the_closure_keys_match_exactly(self):
        _, block = schema_block('ring_cure')
        closure = next(a for a in block['attributes'] if a['name'] == 'closure')
        names = sorted(a['name'] for a in closure['attributes'])
        self.assertEqual(names, sorted(RingController.defaults['closure']))

    def test_scalar_defaults_agree_with_the_controller(self):
        # a schema default that disagrees documents behavior the code does not have
        attrs, _ = schema_block('ring_cure')
        for k, v in attrs.items():
            if isinstance(RingController.defaults[k], (int, float, bool, str)):
                self.assertEqual(v, RingController.defaults[k], f'default for {k}')

    def test_max_accept_is_not_still_in_cures_drag_block(self):
        _, cure = schema_block('CURE')
        drag = next(a for a in cure['attributes'] if a['name'] == 'drag')
        self.assertNotIn('max_accept', [a['name'] for a in drag['attributes']])


class TestTheCureIsNoLongerConstantVolume(unittest.TestCase):
    """233 triazines formed in a box fixed to four decimals across 22 iterations: a cure
    with no constant-pressure stage cannot densify, and cure shrinkage cannot be
    reproduced at fixed volume by construction."""

    def test_relax_ends_at_constant_pressure(self):
        self.assertEqual(RingController.defaults['relax'][-1]['ensemble'], 'npt')

    def test_it_is_above_tg_or_it_will_not_densify(self):
        npt = RingController.defaults['relax'][-1]
        self.assertGreaterEqual(npt['temperature'], 500)
        self.assertEqual(npt['pressure'], 1)

    def test_the_closure_ladder_stays_constant_volume(self):
        # NPT while restraints are actively pulling groups together is not stable
        stages = RingController.defaults['closure']['equilibration']
        self.assertNotIn('npt', [s['ensemble'] for s in stages])

    def test_it_does_not_cool_before_the_barostat(self):
        # the ladder runs at 600 K and so does the next one; relaxing at 300 K in
        # between cooled the network and handed the barostat a cold configuration
        stages = [s for s in RingController.defaults['relax'] if s['ensemble'] != 'min']
        self.assertEqual({s['temperature'] for s in stages}, {600})

    def test_it_mirrors_cures_relaxation(self):
        _, cure = schema_block('CURE')
        relax = next(a for a in cure['attributes'] if a['name'] == 'relax')
        eq = next(a for a in relax['attributes'] if a['name'] == 'equilibration')
        self.assertEqual(RingController.defaults['relax'], eq['default'])

    def test_equilibrate_matches_cures(self):
        # "do for the ring cure exactly what is already done for the binary cures" --
        # same stage, same defaults, not a new design
        _, cure = schema_block('CURE')
        theirs = next(a for a in cure['attributes'] if a['name'] == 'equilibrate')
        _, ring = schema_block('ring_cure')
        mine = next(a for a in ring['attributes'] if a['name'] == 'equilibrate')
        keyed = lambda b: {a['name']: a.get('default') for a in b['attributes']}
        self.assertEqual(keyed(mine), keyed(theirs))

    def test_the_cure_now_has_a_cold_stage_to_pull_the_box_back(self):
        # relax ends hot and constant-pressure; with nothing cold the box only grows
        eq = RingController.defaults['equilibrate']
        self.assertEqual(eq['ensemble'], 'npt')
        self.assertLess(eq['temperature'], RingController.defaults['relax'][-1]['temperature'])

    def test_the_schema_default_matches_the_controller(self):
        attrs, _ = schema_block('ring_cure')
        self.assertEqual(attrs['relax'], RingController.defaults['relax'])
