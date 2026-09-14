"""

.. module:: test_completion_bias
   :synopsis: tests the CURE.controls.completion_bias candidate ranking

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
logger=logging.getLogger(__name__)
import pandas as pd
from htpolynet.cure.curecontroller import CureController, rank_bond_candidates, residue_reaction_counts, residue_functionality

def candidates(rows):
    """Builds a minimal bond-candidate frame: (rj, r) pairs with the columns
    the ranking actually reads, plus one it must carry through untouched."""
    return pd.DataFrame({
        'ai':   [10*i+1 for i in range(len(rows))],
        'aj':   [10*i+2 for i in range(len(rows))],
        'ri':   [100+i  for i in range(len(rows))],
        'rj':   [r[0]   for r in rows],
        'r':    [r[1]   for r in rows],
        'order':[1      for _ in rows],
    })

class TestRankBondCandidates(unittest.TestCase):
    def test_no_counts_is_distance_only(self):
        bdf=candidates([(1,0.44),(2,0.31),(3,0.52),(4,0.29)])
        ranked=rank_bond_candidates(bdf)
        self.assertEqual(ranked['r'].to_list(),[0.29,0.31,0.44,0.52])
        self.assertEqual(ranked['rj'].to_list(),[4,2,1,3])

    def test_bias_promotes_the_most_complete_crosslinker(self):
        # residue 3 already carries 2 of its 3 bonds but is the farthest
        # candidate; residue 4 is unreacted and the nearest.  Distance alone
        # picks 4 -- which is exactly the mechanism that spreads bonds over
        # every crosslinker instead of completing any of them.
        bdf=candidates([(1,0.44),(2,0.31),(3,0.52),(4,0.29)])
        counts=pd.Series({1:1,2:0,3:2,4:0})
        ranked=rank_bond_candidates(bdf,counts)
        self.assertEqual(ranked['rj'].to_list(),[3,1,4,2])

    def test_bias_breaks_ties_by_distance_within_a_group(self):
        bdf=candidates([(1,0.44),(2,0.31),(3,0.52),(4,0.29)])
        counts=pd.Series({1:2,2:2,3:1,4:1})
        ranked=rank_bond_candidates(bdf,counts)
        # both 2-bond residues first, nearer one leading; then both 1-bond
        self.assertEqual(ranked['rj'].to_list(),[2,1,4,3])

    def test_residue_absent_from_counts_ranks_as_unreacted(self):
        bdf=candidates([(1,0.44),(9,0.10)])
        counts=pd.Series({1:2})
        ranked=rank_bond_candidates(bdf,counts)
        self.assertEqual(ranked['rj'].to_list(),[1,9])

    def test_ranking_leaves_the_schema_alone(self):
        bdf=candidates([(1,0.44),(2,0.31)])
        counts=pd.Series({1:2,2:0})
        ranked=rank_bond_candidates(bdf,counts)
        self.assertEqual(list(ranked.columns),list(bdf.columns))
        self.assertEqual(list(ranked.index),[0,1])

class TestResidueReactionCounts(unittest.TestCase):
    def test_counts_sum_over_a_residues_atoms(self):
        # a triazine's three ring carbons each carry at most one bond, so the
        # residue sum is the number of bridges attached to it
        adf=pd.DataFrame({
            'resNum':    [1,1,1,1,1,1, 2,2,2,2,2,2],
            'nreactions':[1,1,0,0,0,0, 1,0,0,0,0,0],
        })
        counts=residue_reaction_counts(adf)
        self.assertEqual(counts[1],2)
        self.assertEqual(counts[2],1)

class TestCompletionBiasControl(unittest.TestCase):
    def test_default_is_off(self):
        cc=CureController({})
        self.assertFalse(cc.dicts['controls']['completion_bias'])
        adf=pd.DataFrame({'resNum':[1],'nreactions':[0]})
        self.assertIsNone(cc._completion_bias_counts(adf))

    def test_enabled_yields_counts(self):
        cc=CureController({'controls':{'completion_bias':True}})
        adf=pd.DataFrame({'resNum':[1,1,2],'nreactions':[1,1,0]})
        counts=cc._completion_bias_counts(adf)
        self.assertIsNotNone(counts)
        self.assertEqual(counts[1],2)

    def test_missing_nreactions_falls_back_to_distance(self):
        # a .grx written before this attribute existed must not crash a cure
        cc=CureController({'controls':{'completion_bias':True}})
        adf=pd.DataFrame({'resNum':[1,2]})
        with self.assertLogs('htpolynet.cure.curecontroller',level='WARNING'):
            self.assertIsNone(cc._completion_bias_counts(adf))

def atoms(rows):
    """Builds an atom frame from (resNum, z, nreactions) triples."""
    return pd.DataFrame({
        'resNum':    [r[0] for r in rows],
        'z':         [r[1] for r in rows],
        'nreactions':[r[2] for r in rows],
    })

class TestResidueFunctionality(unittest.TestCase):
    def test_sum_is_conserved_as_the_cure_proceeds(self):
        # a triazine starts with three z=1 ring carbons; two have since
        # reacted.  Its functionality must still read 3.
        fresh=atoms([(1,1,0),(1,1,0),(1,1,0)])
        partly=atoms([(1,0,1),(1,0,1),(1,1,0)])
        self.assertEqual(residue_functionality(fresh)[1],3)
        self.assertEqual(residue_functionality(partly)[1],3)

    def test_bridge_and_crosslinker_are_distinguished(self):
        adf=atoms([(1,1,0),(1,1,0),          # a difunctional bridge
                   (2,1,0),(2,1,0),(2,1,0)]) # a trifunctional crosslinker
        func=residue_functionality(adf)
        self.assertEqual(func[1],2)
        self.assertEqual(func[2],3)

class TestBiasSideCheck(unittest.TestCase):
    """The failure mode is a crosslinker declared as reactant A: the bias then
    completes bridges instead.  It cannot be caught by looking for an all-zero
    bias key, because a difunctional bridge accumulates reactions too."""

    def setUp(self):
        self.cc=CureController({'controls':{'completion_bias':True}})
        # residue 1 is a difunctional bridge, residue 2 a trifunctional crosslinker
        self.adf=atoms([(1,1,0),(1,1,0), (2,1,0),(2,1,0),(2,1,0)])

    def test_silent_when_the_crosslinker_is_on_the_b_side(self):
        bdf=pd.DataFrame({'ri':[1],'rj':[2],'r':[0.3]})
        with self.assertNoLogs('htpolynet.cure.curecontroller',level='WARNING'):
            self.cc._check_bias_side(self.adf,bdf)

    def test_warns_when_the_crosslinker_is_on_the_a_side(self):
        bdf=pd.DataFrame({'ri':[2],'rj':[1],'r':[0.3]})
        with self.assertLogs('htpolynet.cure.curecontroller',level='WARNING') as cm:
            self.cc._check_bias_side(self.adf,bdf)
        self.assertIn('3 reactive sites vs 2',''.join(cm.output))

    def test_warns_only_once_per_run(self):
        bdf=pd.DataFrame({'ri':[2],'rj':[1],'r':[0.3]})
        with self.assertLogs('htpolynet.cure.curecontroller',level='WARNING'):
            self.cc._check_bias_side(self.adf,bdf)
        with self.assertNoLogs('htpolynet.cure.curecontroller',level='WARNING'):
            self.cc._check_bias_side(self.adf,bdf)

    def test_symmetric_functionality_is_not_a_complaint(self):
        adf=atoms([(1,1,0),(1,1,0), (2,1,0),(2,1,0)])
        bdf=pd.DataFrame({'ri':[2],'rj':[1],'r':[0.3]})
        with self.assertNoLogs('htpolynet.cure.curecontroller',level='WARNING'):
            self.cc._check_bias_side(adf,bdf)

    def test_no_crash_when_the_grx_predates_these_attributes(self):
        adf=pd.DataFrame({'resNum':[1,2]})
        bdf=pd.DataFrame({'ri':[2],'rj':[1],'r':[0.3]})
        self.cc._check_bias_side(adf,bdf)
        self.assertFalse(self.cc.bias_side_checked)

class TestIterationsVsFunctionality(unittest.TestCase):
    """A system whose crosslinkers did not complete has no junctions, and
    nothing else says so.  Two routes there: too few iterations to be possible
    at all, and enough iterations but a last one that formed almost nothing.
    The second is not bounded away by any iteration count, so the completed
    count is checked directly."""

    class FakeTC:
        def __init__(self, adf): self._adf = adf
        def gro_DataFrame(self, name): return self._adf if name == 'atoms' else None

    def controller(self, iterations):
        cc = CureController({})
        cc.state.iter = iterations
        return cc

    def trifunctional(self):
        return pd.DataFrame({
            'resNum':    [1, 1, 2, 2, 2],
            'resName':   ['BPA', 'BPA', 'TAZ', 'TAZ', 'TAZ'],
            'z':         [1, 1, 1, 1, 1],
            'nreactions':[0, 0, 0, 0, 0],
        })

    def complete_trifunctional(self):
        adf = self.trifunctional()
        adf.loc[2:4, ['z', 'nreactions']] = [0, 1]   # TAZ: all three sites spent
        return adf

    def taz_population(self, n_taz, n_complete):
        """n_taz trifunctional TAZ, n_complete of them fully reacted."""
        rows = []
        for i in range(n_taz):
            done = i < n_complete
            for _ in range(3):
                rows.append({'resNum': i + 1, 'resName': 'TAZ',
                             'z': 0 if done else 1, 'nreactions': 1 if done else 0})
        return pd.DataFrame(rows)

    def styrene_chains(self, n_sty, n_interior):
        """n_sty difunctional STY, n_interior of them reacted at both sites."""
        rows = []
        for i in range(n_sty):
            done = i < n_interior
            for _ in range(2):
                rows.append({'resNum': i + 1, 'resName': 'STY',
                             'z': 0 if done else 1, 'nreactions': 1 if done else 0})
        return pd.DataFrame(rows)

    def test_silent_for_a_linear_system_that_is_not_fully_reacted(self):
        # the example-1 false positive: 902 of 1000 STY mid-chain.  The
        # gel-point threshold (1/(f-1))**f is 100% at f=2, and a difunctional
        # unit is a chain interior, not a junction, so neither check applies
        cc = self.controller(20)
        with self.assertNoLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(self.styrene_chains(1000, 902)))

    def test_silent_for_a_linear_system_with_too_few_iterations(self):
        cc = self.controller(1)
        with self.assertNoLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(self.styrene_chains(10, 0)))

    def test_a_trifunctional_residue_still_triggers_the_check(self):
        # a difunctional monomer alongside a trifunctional crosslinker is a
        # network system; the f>=3 residue decides, and it is still checked
        adf = pd.concat([self.styrene_chains(4, 4), self.taz_population(40, 1)
                         .assign(resNum=lambda d: d.resNum + 100)], ignore_index=True)
        cc = self.controller(20)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING') as cm:
            cc.check_iterations_vs_functionality(self.FakeTC(adf))
        self.assertIn('do not assume this system percolates', ''.join(cm.output))

    def test_warns_when_iterations_are_below_functionality(self):
        cc = self.controller(2)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING') as cm:
            cc.check_iterations_vs_functionality(self.FakeTC(self.trifunctional()))
        text = ''.join(cm.output)
        self.assertIn('TAZ', text)
        self.assertIn('no crosslink junctions', text)
        self.assertIn('completion_bias does not change this', text)

    def test_warns_at_the_functionality_when_nothing_completed(self):
        # the route the iteration-count test passes in silence: n == f, so
        # completion is possible, but the last iteration formed almost nothing
        cc = self.controller(3)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING') as cm:
            cc.check_iterations_vs_functionality(self.FakeTC(self.trifunctional()))
        text = ''.join(cm.output)
        self.assertIn('0 of 1 TAZ', text)
        self.assertIn('do not assume this system percolates', text)

    def test_warns_with_plenty_of_iterations_when_nothing_completed(self):
        cc = self.controller(9)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(self.trifunctional()))

    def test_silent_when_the_crosslinkers_completed(self):
        adf = self.complete_trifunctional()
        cc = self.controller(3)
        with self.assertNoLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(adf))

    def test_one_in_240_warns(self):
        # the observed build: n == f, target conversion reached, and a single
        # complete crosslinker out of 240
        cc = self.controller(3)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING') as cm:
            cc.check_iterations_vs_functionality(self.FakeTC(self.taz_population(240, 1)))
        self.assertIn('1 of 240 TAZ (0.4%)', ''.join(cm.output))

    def test_silent_above_the_gel_point(self):
        # 60 of 240 is 25%, above the 12.5% an ideal f=3 network shows at its
        # gel point, so the completed count is not remarked on
        cc = self.controller(3)
        with self.assertNoLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(self.taz_population(240, 60)))

    def test_counting_message_still_wins_below_the_functionality(self):
        # below f the cause is known exactly, so say that rather than quoting
        # a fraction that could not have been anything else
        cc = self.controller(2)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING') as cm:
            cc.check_iterations_vs_functionality(self.FakeTC(self.taz_population(240, 0)))
        self.assertIn('at most one bond per residue per iteration', ''.join(cm.output))

    def test_counts_sites_already_spent(self):
        # functionality is z + nreactions, so a partly-reacted crosslinker
        # still counts as trifunctional and still triggers the check
        adf = self.trifunctional()
        adf.loc[2:3, ['z', 'nreactions']] = [0, 1]
        cc = self.controller(2)
        with self.assertLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(adf))

    def test_silent_when_nothing_is_multifunctional(self):
        adf = pd.DataFrame({'resNum': [1, 2], 'resName': ['A', 'B'],
                            'z': [1, 1], 'nreactions': [0, 0]})
        cc = self.controller(1)
        with self.assertNoLogs('htpolynet.cure.curecontroller', level='WARNING'):
            cc.check_iterations_vs_functionality(self.FakeTC(adf))

    def test_no_crash_without_the_attributes(self):
        cc = self.controller(2)
        cc.check_iterations_vs_functionality(self.FakeTC(pd.DataFrame({'resNum': [1]})))
