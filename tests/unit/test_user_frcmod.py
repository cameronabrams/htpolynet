"""

.. module:: test_user_frcmod
   :synopsis: a constituent's frcmod reaches its own parameterization, every template containing it,
              and the parameter-cache record, and nothing else

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import shutil
from types import SimpleNamespace

import pandas as pd
import pytest

import htpolynet.core.projectfilesystem as pfs
from htpolynet.core import paramcache
from htpolynet.core.molecule import Molecule
from htpolynet.core.runtime import Runtime
from htpolynet.core.topology import Topology, conflicting_types

FRC = ('BAF', 'MASS\n\nDIHE\nca-ca-c3-ca   1     0.3190    180.0    -2.0\n\n')


class TestCacheRecord:
    def test_no_frcmod_leaves_the_record_as_it_was(self):
        assert paramcache.FRCMOD_FIELD not in paramcache.build_key({'charge_method': 'bcc'})

    def test_frcmod_is_recorded_by_content(self):
        a = paramcache.build_key({'frcmod': [FRC]})
        b = paramcache.build_key({'frcmod': [(FRC[0], FRC[1] + 'x')]})
        assert a[paramcache.FRCMOD_FIELD] != b[paramcache.FRCMOD_FIELD]
        assert a == paramcache.build_key({'frcmod': [FRC]})

    def test_plain_cache_entry_is_a_miss_when_a_frcmod_is_requested(self):
        diffs = paramcache.describe_mismatch(paramcache.build_key({}), paramcache.build_key({'frcmod': [FRC]}))
        assert any('frcmod' in d for d in diffs)

    def test_overridden_cache_entry_is_a_miss_when_none_is_requested(self):
        assert paramcache.describe_mismatch(paramcache.build_key({'frcmod': [FRC]}), paramcache.build_key({}))

    def test_edited_frcmod_is_a_miss(self):
        assert paramcache.describe_mismatch(paramcache.build_key({'frcmod': [FRC]}),
                                            paramcache.build_key({'frcmod': [(FRC[0], FRC[1] + '\n')]}))

    def test_record_from_before_frcmods_matches_a_run_without_one(self):
        stored = {'charge_method': 'bcc', 'net_charge': 0, 'atom_type': 'gaff'}
        assert paramcache.describe_mismatch(stored, paramcache.build_key({})) == []


class TestWhichMoleculesGetIt:
    def runtime(self):
        r = SimpleNamespace(cfg=SimpleNamespace(ambertools={'atom_type': 'gaff2'}), frcmods={'BAF': FRC[1]})
        return r

    def test_the_monomer(self):
        amb = Runtime._ambertools_for(self.runtime(), SimpleNamespace(name='BAF', sequence=['BAF']))
        assert amb['frcmod'] == [FRC]
        assert amb['atom_type'] == 'gaff2'

    def test_every_template_containing_it(self):
        M = SimpleNamespace(name='BAF~O1-C1~BAF~O1-C2~TAZ', sequence=['BAF', 'BAF', 'TAZ'])
        assert Runtime._ambertools_for(self.runtime(), M)['frcmod'] == [FRC]

    def test_not_a_molecule_without_it(self):
        assert 'frcmod' not in Runtime._ambertools_for(self.runtime(), SimpleNamespace(name='TAZ', sequence=['TAZ']))

    def test_configuration_directives_are_not_mutated(self):
        r = self.runtime()
        Runtime._ambertools_for(r, SimpleNamespace(name='BAF', sequence=['BAF']))
        assert 'frcmod' not in r.cfg.ambertools


class TestReadFrcmods:
    def runtime(self, constituents):
        return SimpleNamespace(cfg=SimpleNamespace(constituents=constituents))

    def test_reads_relative_to_the_configuration(self, tmp_path):
        (tmp_path / 'baf.frcmod').write_text(FRC[1])
        pfs.pfs_setup(root=str(tmp_path), mock=True)
        got = Runtime._read_frcmods(self.runtime({'BAF': {'frcmod': 'baf.frcmod'}, 'TAZ': {'count': 1}}))
        assert got == {'BAF': FRC[1]}

    def test_a_missing_file_stops_setup(self, tmp_path):
        pfs.pfs_setup(root=str(tmp_path), mock=True)
        with pytest.raises(FileNotFoundError, match='BAF'):
            Runtime._read_frcmods(self.runtime({'BAF': {'frcmod': 'nope.frcmod'}}))


def dihedraltypes(rows):
    T = Topology()
    T.D['dihedraltypes'] = pd.DataFrame(rows, columns=['i', 'j', 'k', 'l', 'func', 'phase', 'kd', 'pn'])
    return T


class TestConflictingTypes:
    plain = [('ca', 'c3', 'ca', 'ca', 9, 0.0, 0.0, 2)]
    fitted = [('ca', 'c3', 'ca', 'ca', 9, 180.0, 1.335, 2), ('ca', 'c3', 'ca', 'ca', 9, 180.0, 0.959, 4)]

    def test_overridden_type_also_present_without_the_override(self):
        c = conflicting_types({'BAF': (dihedraltypes(self.fitted), 'd1'), 'BPA': (dihedraltypes(self.plain), None)})
        assert c == [('dihedraltypes', ('ca', 'c3', 'ca', 'ca', 9), 'BAF', 'BPA')]

    def test_an_improper_with_the_same_atom_types_is_a_different_type(self):
        improper = [('ca', 'c3', 'ca', 'ca', 4, 180.0, 0.263, 2)]
        both = self.fitted + improper
        assert conflicting_types({'BAF': (dihedraltypes(both), 'd1'), 'X': (dihedraltypes(improper), None)}) == []

    def test_reversed_type_names_are_the_same_type(self):
        rev = [('ca', 'ca', 'c3', 'ca', 9, 0.0, 0.0, 2)]
        assert conflicting_types({'BAF': (dihedraltypes(self.fitted), 'd1'), 'BPA': (dihedraltypes(rev), None)})

    def test_same_group_differences_are_not_this_check(self):
        # discrepancies among plain GAFF molecules are resolve_type_discrepancies' business
        assert conflicting_types({'A': (dihedraltypes(self.fitted), None), 'B': (dihedraltypes(self.plain), None)}) == []

    def test_agreeing_parameters_are_fine(self):
        assert conflicting_types({'BAF': (dihedraltypes(self.plain), 'd1'), 'TAZ': (dihedraltypes(self.plain), None)}) == []


_MISSING = [t for t in ('antechamber', 'parmchk2', 'tleap', 'gmx') if shutil.which(t) is None]


@pytest.mark.skipif(bool(_MISSING), reason=f'external tools not on PATH: {", ".join(_MISSING)}')
class TestTleapLoadsIt:
    def parameterize(self, root, ambertools):
        orig = os.getcwd()
        pfs.pfs_setup(root=str(root), projdir='proj', mock=False)
        try:
            M = Molecule.New('STY', None)
            M.set_sequence_from_moldict({'STY': M})
            M.generate(ambertools=ambertools, gaff={'minimize_molecules': False})
            return Topology.read_top('STY.top').D['dihedraltypes'], paramcache.read_key('STY')
        finally:
            os.chdir(orig)

    def test_the_override_replaces_the_gaff_torsion(self, tmp_path_factory):
        plain, plain_key = self.parameterize(tmp_path_factory.mktemp('plain'), {'charge_method': 'gas'})
        t = plain.iloc[0]
        names = '-'.join(str(x).ljust(2) for x in (t.i, t.j, t.k, t.l))
        frcmod = f'override\nDIHE\n{names}   1    7.7700    180.0     3.0\n\n'
        over, over_key = self.parameterize(tmp_path_factory.mktemp('over'),
                                           {'charge_method': 'gas', 'frcmod': [('STY', frcmod)]})
        rows = over[(over.i == t.i) & (over.j == t.j) & (over.k == t.k) & (over.l == t.l)]
        assert list(rows.kd.round(3)) == [round(7.77 * 4.184, 3)]
        assert list(rows.pn) == [3]
        assert paramcache.FRCMOD_FIELD not in plain_key
        assert over_key[paramcache.FRCMOD_FIELD] == paramcache.frcmod_digest([('STY', frcmod)])
