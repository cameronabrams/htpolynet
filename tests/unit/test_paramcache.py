"""Tests the provenance record that guards the parameterization cache.

The cache under ``molecules/parameterized`` is keyed on molecule name alone.
Nothing in that key reflects the charge method, so a configuration asking for
``charge_method: bcc`` used to silently reuse a ``gas`` parameterization
checked in months earlier under the same residue name -- putting two different
charge methods into one network with nothing in the output saying so.  These
tests pin the record that makes that a cache miss instead.

None of this needs AmberTools: the record is compared, not regenerated.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import json
import os
import pathlib
import shutil

from types import SimpleNamespace

import pytest

import htpolynet.core.projectfilesystem as pfs

from htpolynet.core import paramcache
from htpolynet.core.runtime import Runtime
from htpolynet.external.ambertools import AMBERTOOLS_DEFAULTS


@pytest.fixture
def library(tmp_path, monkeypatch):
    """A stand-in library that pfs.exists and pfs.checkout resolve against.

    Returns the directory that plays the part of ``molecules/parameterized``
    in the library, with the working directory set to a separate project
    directory, which is the arrangement the real code sees.
    """
    lib = tmp_path / 'lib'
    (lib / pfs.Dirs.molecules_parameterized).mkdir(parents=True)
    work = tmp_path / 'work'
    work.mkdir()
    monkeypatch.chdir(work)

    def _exists(filename):
        return (lib / filename).exists()

    def _checkout(filename, altpath=[]):
        src = lib / filename
        if not src.exists():
            return False
        shutil.copy(src, os.path.basename(filename))
        return True

    monkeypatch.setattr(pfs, 'exists', _exists)
    monkeypatch.setattr(pfs, 'checkout', _checkout)
    return lib / pfs.Dirs.molecules_parameterized


def _library_entry(library, name, **ambertools):
    """Writes a provenance record into the library for the named molecule."""
    paramcache.write_key(str(library / name), paramcache.build_key(ambertools))


def _runtime(**ambertools):
    """A stand-in for Runtime carrying only what the check reads off it."""
    r = SimpleNamespace(cfg=SimpleNamespace(ambertools=ambertools),
                        unverified_parameterizations=[], frcmods={})
    r._ambertools_for = lambda M: Runtime._ambertools_for(r, M)
    return r


def _mismatch(name='TAZ', **ambertools):
    return Runtime._cached_parameterization_mismatch(_runtime(**ambertools),
                                                     SimpleNamespace(name=name))


def _mismatch_and_runtime(name='TAZ', **ambertools):
    """As _mismatch, but also returns the stand-in so the summary list can be
    inspected."""
    r = _runtime(**ambertools)
    return Runtime._cached_parameterization_mismatch(r, SimpleNamespace(name=name)), r


class TestBuildKey:
    def test_defaults_come_from_ambertools(self):
        assert paramcache.build_key({}, version=None) == AMBERTOOLS_DEFAULTS

    def test_none_is_the_same_as_empty(self):
        assert paramcache.build_key(None) == paramcache.build_key({})

    def test_records_the_requested_directives(self):
        key = paramcache.build_key({'charge_method': 'gas', 'net_charge': -1})
        assert key['charge_method'] == 'gas'
        assert key['net_charge'] == -1
        assert key['atom_type'] == AMBERTOOLS_DEFAULTS['atom_type']

    def test_ignores_directives_that_do_not_affect_parameters(self):
        assert paramcache.build_key({'nonsense': 1}, version=None) == AMBERTOOLS_DEFAULTS

    def test_net_charge_default_preserves_the_old_hardcoded_nc(self):
        # antechamber was invoked with a literal -nc 0 before net_charge became
        # a directive; the default has to keep producing that command.
        assert AMBERTOOLS_DEFAULTS['net_charge'] == 0


class TestReadWrite:
    def test_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        key = paramcache.build_key({'charge_method': 'bcc'})
        paramcache.write_key('TAZ', key)
        assert paramcache.read_key('TAZ') == key

    def test_written_beside_the_files_it_describes(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        paramcache.write_key('TAZ', paramcache.build_key({}))
        assert (tmp_path / f'TAZ.{paramcache.CACHE_KEY_EXT}').exists()
        assert paramcache.key_filename('TAZ') == f'TAZ.{paramcache.CACHE_KEY_EXT}'

    def test_absent_record_reads_as_none(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert paramcache.read_key('TAZ') is None

    def test_malformed_record_reads_as_none(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / 'TAZ.parm').write_text('{not json')
        assert paramcache.read_key('TAZ') is None

    def test_non_object_record_reads_as_none(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / 'TAZ.parm').write_text(json.dumps(['gas']))
        assert paramcache.read_key('TAZ') is None


class TestDescribeMismatch:
    def test_agreement_is_no_mismatch(self):
        key = paramcache.build_key({'charge_method': 'bcc'})
        assert paramcache.describe_mismatch(key, key) == []

    def test_charge_method_difference_names_both_methods(self):
        diffs = paramcache.describe_mismatch(paramcache.build_key({'charge_method': 'gas'}),
                                             paramcache.build_key({'charge_method': 'bcc'}))
        assert len(diffs) == 1
        assert "'gas'" in diffs[0] and "'bcc'" in diffs[0]

    def test_every_differing_directive_is_reported(self):
        diffs = paramcache.describe_mismatch(
            paramcache.build_key({'charge_method': 'gas', 'net_charge': 0}),
            paramcache.build_key({'charge_method': 'bcc', 'net_charge': -1}))
        assert len(diffs) == 2

    def test_absent_record_is_not_a_mismatch(self):
        assert paramcache.describe_mismatch(None, paramcache.build_key({})) == []

    def test_field_absent_from_an_older_record_is_skipped(self):
        # Adding a field must not invalidate every record already written.
        stored = {'charge_method': 'bcc'}
        assert paramcache.describe_mismatch(stored, paramcache.build_key({'charge_method': 'bcc'})) == []


class TestCachedParameterizationMismatch:
    """The check as Runtime calls it, against a stand-in library."""

    def test_gas_cache_is_rejected_when_bcc_is_requested(self, library):
        # The reported bug: TAZ cached under 'gas', config asks for 'bcc'.
        _library_entry(library, 'TAZ', charge_method='gas')
        diffs = _mismatch('TAZ', charge_method='bcc')
        assert diffs, 'a gas parameterization must not satisfy a bcc request'
        assert any("'gas'" in d and "'bcc'" in d for d in diffs)

    def test_matching_cache_is_accepted(self, library):
        _library_entry(library, 'TAZ', charge_method='bcc')
        assert _mismatch('TAZ', charge_method='bcc') == []

    def test_net_charge_difference_is_rejected(self, library):
        _library_entry(library, 'ION', charge_method='bcc', net_charge=0)
        assert _mismatch('ION', charge_method='bcc', net_charge=-1)

    def test_record_absent_is_accepted_with_a_warning(self, library, caplog):
        # A library built before records existed must keep working rather than
        # re-parameterizing wholesale, but the run has to say so.
        with caplog.at_level('WARNING'):
            diffs, runtime = _mismatch_and_runtime('TAZ', charge_method='bcc')
        assert diffs == []
        assert 'TAZ' in caplog.text
        assert 'force-parameterization' in caplog.text
        # The message must say the cached values are being USED, not merely
        # that they could not be checked: correct behavior here is otherwise
        # indistinguishable from the bug it replaced.
        assert 'Reusing cached parameterization' in caplog.text
        assert 'may not be' in caplog.text
        assert runtime.unverified_parameterizations == ['TAZ'], \
            'an unverifiable entry must also reach the stage-end summary'

    def test_a_verified_entry_is_not_reported_as_unverified(self, library):
        _library_entry(library, 'TAZ', charge_method='bcc')
        _, runtime = _mismatch_and_runtime('TAZ', charge_method='bcc')
        assert runtime.unverified_parameterizations == []

    def test_a_mismatch_is_not_reported_as_unverified(self, library):
        # A mismatch re-parameterizes, so it is handled, not unverified.
        _library_entry(library, 'TAZ', charge_method='gas')
        diffs, runtime = _mismatch_and_runtime('TAZ', charge_method='bcc')
        assert diffs
        assert runtime.unverified_parameterizations == []

    def test_default_charge_method_matches_a_default_record(self, library):
        _library_entry(library, 'TAZ')
        assert _mismatch('TAZ') == []

    def test_leftover_project_record_cannot_mask_the_library(self, library):
        # A record left in the project directory by an earlier run describes
        # that run, not the library entry previously_parameterized() found.
        _library_entry(library, 'TAZ', charge_method='gas')
        paramcache.write_key('TAZ', paramcache.build_key({'charge_method': 'bcc'}))
        assert _mismatch('TAZ', charge_method='bcc'), \
            'the library record, not a stale local one, decides the question'

    def test_leftover_project_record_is_not_trusted_when_library_has_none(self, library, caplog):
        paramcache.write_key('TAZ', paramcache.build_key({'charge_method': 'bcc'}))
        with caplog.at_level('WARNING'):
            assert _mismatch('TAZ', charge_method='bcc') == []
        assert 'carries no record of how it was built' in caplog.text


class TestCheckinParameterization:
    """A record must never be checked in without the data it describes."""

    @pytest.fixture
    def cache(self, tmp_path, monkeypatch):
        """A stand-in user library that pfs.exists and pfs.checkin write to."""
        lib = tmp_path / 'lib' / pfs.Dirs.molecules_parameterized
        lib.mkdir(parents=True)
        work = tmp_path / 'work'
        work.mkdir()
        monkeypatch.chdir(work)

        def _exists(filename):
            return (tmp_path / 'lib' / filename).exists()

        def _checkin(filename, overwrite=False):
            src = pathlib.Path(os.path.basename(filename))
            if not src.exists():
                return
            dest = tmp_path / 'lib' / filename
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists() or overwrite:
                shutil.copy(src, dest)

        monkeypatch.setattr(pfs, 'exists', _exists)
        monkeypatch.setattr(pfs, 'checkin', _checkin)
        return lib

    @staticmethod
    def _produce(name, contents, **ambertools):
        """Writes the files a parameterization run leaves in the project dir."""
        for ex in ('mol2', 'top', 'tpx', 'itp', 'gro', 'grx'):
            pathlib.Path(f'{name}.{ex}').write_text(contents)
        paramcache.write_key(name, paramcache.build_key(ambertools))

    def test_new_entry_gets_its_record(self, cache):
        self._produce('TAZ', 'bcc-data', charge_method='bcc')
        Runtime._checkin_parameterization('TAZ', False)
        assert (cache / 'TAZ.gro').read_text() == 'bcc-data'
        assert json.loads((cache / 'TAZ.parm').read_text())['charge_method'] == 'bcc'

    def test_record_is_not_attached_to_data_that_was_kept(self, cache):
        # The library holds a gas entry with no record.  A bcc run that is not
        # allowed to overwrite must not leave a bcc record on the gas files.
        (cache / 'TAZ.gro').write_text('gas-data')
        self._produce('TAZ', 'bcc-data', charge_method='bcc')
        Runtime._checkin_parameterization('TAZ', False)
        assert (cache / 'TAZ.gro').read_text() == 'gas-data', 'data must be kept'
        assert not (cache / 'TAZ.parm').exists(), \
            'a record must not describe files it did not accompany'

    def test_force_checkin_replaces_data_and_record_together(self, cache):
        (cache / 'TAZ.gro').write_text('gas-data')
        paramcache.write_key(str(cache / 'TAZ'), paramcache.build_key({'charge_method': 'gas'}))
        self._produce('TAZ', 'bcc-data', charge_method='bcc')
        Runtime._checkin_parameterization('TAZ', True)
        assert (cache / 'TAZ.gro').read_text() == 'bcc-data'
        assert json.loads((cache / 'TAZ.parm').read_text())['charge_method'] == 'bcc'

    def test_record_and_data_in_the_library_always_agree(self, cache):
        # The property the two tests above are instances of: whatever the
        # library ends up holding, its record either describes its data or is
        # absent -- it is never a record for some other parameterization.
        for overwrite in (False, True):
            (cache / 'TAZ.gro').write_text('gas-data')
            self._produce('TAZ', 'bcc-data', charge_method='bcc')
            Runtime._checkin_parameterization('TAZ', overwrite)
            if (cache / 'TAZ.parm').exists():
                stored = json.loads((cache / 'TAZ.parm').read_text())
                assert stored['charge_method'] == 'bcc'
                assert (cache / 'TAZ.gro').read_text() == 'bcc-data'


class TestUnverifiedSummary:
    """The stage-end block that a person scanning a log will actually hit."""

    @staticmethod
    def _runtime(names, **ambertools):
        r = SimpleNamespace(cfg=SimpleNamespace(ambertools=ambertools),
                            unverified_parameterizations=list(names))
        return r

    def _report(self, names, caplog, **ambertools):
        with caplog.at_level('WARNING'):
            Runtime._report_unverified_parameterizations(self._runtime(names, **ambertools))
        return caplog.text

    def test_silent_when_every_entry_was_verified(self, caplog):
        assert self._report([], caplog, charge_method='bcc') == ''

    def test_names_every_unverified_molecule(self, caplog):
        text = self._report(['TAZ', 'BAF'], caplog, charge_method='bcc')
        assert 'TAZ' in text and 'BAF' in text

    def test_states_the_count_and_the_requested_charge_method(self, caplog):
        text = self._report(['TAZ', 'BAF'], caplog, charge_method='bcc')
        assert '2 parameterizations reused without provenance' in text
        assert "'bcc'" in text

    def test_singular_for_one_molecule(self, caplog):
        text = self._report(['TAZ'], caplog, charge_method='bcc')
        assert '1 parameterization reused without provenance' in text

    def test_says_the_build_carries_the_cached_charges(self, caplog):
        # Without this, a reader who upgrades with an existing library, asks
        # for bcc and gets gas numbers concludes the guard is broken.
        text = self._report(['TAZ'], caplog, charge_method='bcc')
        assert 'carries whatever charges those entries hold' in text
        assert 'may not be' in text

    def test_says_this_is_expected_rather_than_a_failure(self, caplog):
        text = self._report(['TAZ'], caplog, charge_method='bcc')
        assert 'is not a' in text and 'failure' in text

    def test_says_how_to_resolve_it(self, caplog):
        assert 'force-parameterization' in self._report(['TAZ'], caplog, charge_method='bcc')

    def test_a_molecule_reported_twice_is_listed_once(self, caplog):
        text = self._report(['TAZ', 'TAZ'], caplog, charge_method='bcc')
        assert '1 parameterization reused without provenance' in text


class TestAmberToolsVersion:
    """The same directives give different numbers under different AmberTools
    releases, and the CPU and CUDA images once shipped different ones."""

    def test_detected_from_the_running_installation(self, monkeypatch):
        from htpolynet.external import software
        monkeypatch.setitem(software.versions, 'ambertools', 'ver. 26.0 (conda)')
        assert paramcache.ambertools_version() == '26.0'
        assert paramcache.build_key({})[paramcache.VERSION_FIELD] == '26.0'

    def test_unknown_version_is_not_recorded(self, monkeypatch):
        from htpolynet.external import software
        monkeypatch.setitem(software.versions, 'ambertools', 'installed (version unknown)')
        assert paramcache.VERSION_FIELD not in paramcache.build_key({})

    def test_a_different_version_is_a_miss(self):
        diffs = paramcache.describe_mismatch(paramcache.build_key({}, version='24.8'),
                                             paramcache.build_key({}, version='26.0'))
        assert diffs == ['AmberTools 24.8 cached, 26.0 running']

    def test_the_same_version_matches(self):
        key = paramcache.build_key({}, version='26.0')
        assert paramcache.describe_mismatch(key, key) == []

    def test_a_record_without_a_version_still_matches(self):
        # records written before the field existed, like every cache in use today
        assert paramcache.describe_mismatch(paramcache.build_key({}, version=None),
                                            paramcache.build_key({}, version='26.0')) == []

    def test_a_run_that_cannot_tell_still_matches(self):
        assert paramcache.describe_mismatch(paramcache.build_key({}, version='26.0'),
                                            paramcache.build_key({}, version=None)) == []
