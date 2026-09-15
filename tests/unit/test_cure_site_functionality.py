"""

.. module:: test_cure_site_functionality
   :synopsis: functionality for the percolation check is counted on cure-reactive sites only

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
import pandas as pd

from htpolynet.cure.curecontroller import CureController, cure_site_mask
from htpolynet.cure.reaction import Reaction

LOG = 'htpolynet.cure.curecontroller'


class FakeMol:
    def __init__(self, *resnames): self.seq = list(resnames)
    def get_resname(self, resid): return self.seq[resid - 1]


class FakeTC:
    def __init__(self, adf): self._adf = adf
    def gro_DataFrame(self, name): return self._adf if name == 'atoms' else None


def rxn(name, stage, reactants, atoms):
    """atoms: {label: (reactant_idx, resid, atom, z)}; bonds A-B."""
    return Reaction({
        'name': name, 'stage': stage, 'reactants': reactants, 'product': name,
        'atoms': {k: {'reactant': r, 'resid': i, 'atom': a, 'z': z} for k, (r, i, a, z) in atoms.items()},
        'bonds': [{'atoms': ['A', 'B'], 'order': 1}],
    })


def residues(resname, n, atoms_by_residue):
    """n copies of a residue; atoms_by_residue(i) -> list of (atomName, z, nreactions)."""
    rows = []
    for i in range(n):
        for a, z, nr in atoms_by_residue(i):
            rows.append({'resNum': 0, 'resName': resname, 'atomName': a, 'z': z, 'nreactions': nr, '_i': i})
    return rows


def frame(*groups):
    rows = []
    offset = 0
    for g in groups:
        n = max(r['_i'] for r in g) + 1
        for r in g:
            r = dict(r); r['resNum'] = offset + r.pop('_i') + 1
            rows.append(r)
        offset += n
    return pd.DataFrame(rows)


class TestCureSiteMask(unittest.TestCase):
    def test_non_cure_reactions_contribute_no_sites(self):
        RL = [rxn('b1', 'param', {1: 'BPA', 2: 'HIE'}, {'A': (1, 1, 'O1', 1), 'B': (2, 1, 'C4', 2)}),
              rxn('c1', 'cure', {1: 'HIE', 2: 'HIE'}, {'A': (1, 1, 'C1', 1), 'B': (2, 1, 'C2', 1)})]
        MD = {'BPA': FakeMol('BPA'), 'HIE': FakeMol('HIE')}
        adf = pd.DataFrame({'resName': ['HIE', 'HIE', 'HIE', 'BPA'], 'atomName': ['C1', 'C2', 'C4', 'O1']})
        self.assertEqual(cure_site_mask(adf, RL, MD).tolist(), [True, True, False, False])

    def test_oligomer_reactant_resolves_the_residue_through_its_sequence(self):
        # examples 3 and 4 react oligomers like PAC~N1-C1~DGE; resid 2 is DGE
        RL = [rxn('c1', 'cure', {1: 'PAC~N1-C1~DGE', 2: 'PAC'}, {'A': (1, 2, 'C1', 1), 'B': (2, 1, 'N1', 2)})]
        MD = {'PAC~N1-C1~DGE': FakeMol('PAC', 'DGE'), 'PAC': FakeMol('PAC')}
        adf = pd.DataFrame({'resName': ['DGE', 'PAC', 'PAC'], 'atomName': ['C1', 'N1', 'C1']})
        self.assertEqual(cure_site_mask(adf, RL, MD).tolist(), [True, True, False])

    def test_none_without_cure_reactions_or_inputs(self):
        adf = pd.DataFrame({'resName': ['STY'], 'atomName': ['C1']})
        self.assertIsNone(cure_site_mask(adf, None, None))
        self.assertIsNone(cure_site_mask(adf, [rxn('p', 'param', {1: 'STY', 2: 'STY'},
                                                   {'A': (1, 1, 'C1', 1), 'B': (2, 1, 'C2', 1)})],
                                         {'STY': FakeMol('STY')}))

    def test_an_unresolvable_reactant_is_skipped_not_fatal(self):
        RL = [rxn('c1', 'cure', {1: 'GHOST', 2: 'STY'}, {'A': (1, 1, 'C1', 1), 'B': (2, 1, 'C2', 1)})]
        adf = pd.DataFrame({'resName': ['STY', 'STY'], 'atomName': ['C1', 'C2']})
        self.assertEqual(cure_site_mask(adf, RL, {'STY': FakeMol('STY')}).tolist(), [False, True])


class TestPercolationOnCureSites(unittest.TestCase):
    def controller(self, iterations):
        cc = CureController({}); cc.state.iter = iterations
        return cc

    def hie_system(self, n_both, n_one):
        """HIE as built in example 2: cure sites C1/C2, plus C4 carrying a
        param-stage bond and a leftover valence no cure reaction uses."""
        def atoms(i):
            both = i < n_both
            return [('C1', 0, 1), ('C2', 0 if both else 1, 1 if both else 0), ('C4', 1, 1)]
        adf = frame(residues('HIE', n_both + n_one, atoms))
        RL = [rxn('yy', 'cure', {1: 'HIE', 2: 'HIE'}, {'A': (1, 1, 'C1', 1), 'B': (2, 1, 'C2', 1)})]
        return adf, RL, {'HIE': FakeMol('HIE')}

    def test_example2_hie_false_positive_is_gone(self):
        adf, RL, MD = self.hie_system(133, 17)
        cc = self.controller(12)
        with self.assertLogs(LOG, level='INFO') as cm:
            cc.check_iterations_vs_functionality(FakeTC(adf), RL, MD)
        text = '\n'.join(cm.output)
        self.assertNotIn('do not assume this system percolates', text)
        self.assertIn('not assessed', text)

    def test_without_the_reaction_list_the_old_count_still_applies(self):
        # guards the premise of the test above: the same system DID warn
        adf, _, _ = self.hie_system(133, 17)
        with self.assertLogs(LOG, level='WARNING'):
            self.controller(12).check_iterations_vs_functionality(FakeTC(adf))

    def taz_system(self, n_taz, n_complete, symmetric=True):
        def atoms(i):
            done = i < n_complete
            return [(c, 0 if done else 1, 1 if done else 0) for c in ('C1', 'C2', 'C3')]
        adf = frame(residues('TAZ', n_taz, atoms))
        names = ('C1', 'C2', 'C3') if symmetric else ('C1',)
        # the runtime list is symmetry-expanded: one cure reaction per ring carbon
        RL = [rxn(f'o{c}', 'cure', {1: 'BPA', 2: 'TAZ'}, {'A': (1, 1, 'O1', 1), 'B': (2, 1, c, 1)}) for c in names]
        return adf, RL, {'TAZ': FakeMol('TAZ'), 'BPA': FakeMol('BPA')}

    def test_symmetry_expanded_trifunctional_is_still_checked_and_quiet_when_complete(self):
        # example 6: this must not be silenced
        adf, RL, MD = self.taz_system(182, 182)
        with self.assertLogs(LOG, level='INFO') as cm:
            self.controller(9).check_iterations_vs_functionality(FakeTC(adf), RL, MD)
        text = '\n'.join(cm.output)
        self.assertIn('182 of 182 TAZ are fully reacted', text)
        self.assertNotIn('percolates', text)

    def test_symmetry_expanded_trifunctional_still_warns_when_incomplete(self):
        adf, RL, MD = self.taz_system(240, 1)
        with self.assertLogs(LOG, level='WARNING') as cm:
            self.controller(9).check_iterations_vs_functionality(FakeTC(adf), RL, MD)
        self.assertIn('do not assume this system percolates', ''.join(cm.output))

    def test_crosslinker_is_not_confused_with_a_residue_that_only_looked_tied(self):
        # examples 3 and 4: counting every site tied the diepoxide with the
        # diamine at four, so the old denominator mixed both residue types
        pac = residues('PAC', 100, lambda i: [('N1', 0 if i < 83 else 1, 2 if i < 83 else 1),
                                              ('N2', 0 if i < 83 else 1, 2 if i < 83 else 1)])
        dge = residues('DGE', 200, lambda i: [('C1', 0, 1), ('C2', 0, 1), ('O9', 2, 0)])
        adf = frame(pac, dge)
        RL = [rxn('pn1', 'cure', {1: 'PAC', 2: 'DGE'}, {'A': (1, 1, 'N1', 2), 'B': (2, 1, 'C1', 1)}),
              rxn('pn2', 'cure', {1: 'PAC', 2: 'DGE'}, {'A': (1, 1, 'N2', 2), 'B': (2, 1, 'C2', 1)})]
        MD = {'PAC': FakeMol('PAC'), 'DGE': FakeMol('DGE')}
        with self.assertLogs(LOG, level='INFO') as cm:
            self.controller(15).check_iterations_vs_functionality(FakeTC(adf), RL, MD)
        self.assertIn('83 of 100 PAC are fully reacted', '\n'.join(cm.output))
