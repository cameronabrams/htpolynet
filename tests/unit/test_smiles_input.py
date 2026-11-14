"""Tests SMILES-to-mol2 materialization.

Every example in the depot declares its constituents as atom-mapped SMILES,
and rdkit is now a core runtime dependency, so this module sits on the path
of every normal run -- but had no coverage.

The dispatch and skip logic is pure Python and always runs.  Actual structure
generation shells out to obabel (even on the RDKit path, which round-trips
through SDF), so those tests skip when obabel is absent.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import os
import shutil

import pytest

from htpolynet.external.smiles_input import materialize_smiles_inputs, _has_atom_mapping

needs_obabel = pytest.mark.skipif(shutil.which('obabel') is None,
                                  reason='obabel not on PATH')

STYRENE = '[CH2:1]=[CH:2]c1ccccc1'


# --- atom-mapping detection ------------------------------------------------

@pytest.mark.parametrize('smiles', [STYRENE, '[CH:1]#[N:2]', 'C[C:3](C)C'])
def test_atom_mapping_detected(smiles):
    assert _has_atom_mapping(smiles)


@pytest.mark.parametrize('smiles', ['C=Cc1ccccc1', 'CCO', 'c1ccccc1', ''])
def test_atom_mapping_absent(smiles):
    assert not _has_atom_mapping(smiles)


# --- which constituents get materialized ----------------------------------

def test_creates_the_inputs_directory(tmp_path):
    d = tmp_path / 'lib' / 'molecules' / 'inputs'
    materialize_smiles_inputs({}, inputs_dir=str(d))
    assert d.is_dir()


def test_constituents_without_smiles_are_skipped(tmp_path):
    spec = {'STY': {'count': 100}, 'PAC': {'count': 50}}
    assert materialize_smiles_inputs(spec, inputs_dir=str(tmp_path)) == []
    assert list(tmp_path.iterdir()) == []


def test_non_dict_specs_are_skipped(tmp_path):
    """A bare scalar under a constituent name must not raise."""
    assert materialize_smiles_inputs({'STY': 100, 'PAC': None}, inputs_dir=str(tmp_path)) == []


def test_empty_smiles_value_is_skipped(tmp_path):
    assert materialize_smiles_inputs({'STY': {'smiles': ''}}, inputs_dir=str(tmp_path)) == []


def test_existing_mol2_is_never_regenerated(tmp_path):
    """Hand-edited mol2 files must survive a re-run; delete to regenerate."""
    existing = tmp_path / 'STY.mol2'
    existing.write_text('HAND EDITED — DO NOT CLOBBER\n')
    spec = {'STY': {'smiles': STYRENE, 'reactive_atoms': {1: 'C1', 2: 'C2'}}}
    assert materialize_smiles_inputs(spec, inputs_dir=str(tmp_path)) == []
    assert existing.read_text() == 'HAND EDITED — DO NOT CLOBBER\n'


# --- actual generation -----------------------------------------------------

@needs_obabel
def test_generates_a_mol2_from_mapped_smiles(tmp_path):
    spec = {'STY': {'smiles': STYRENE, 'reactive_atoms': {1: 'C1', 2: 'C2'}}}
    assert materialize_smiles_inputs(spec, inputs_dir=str(tmp_path)) == ['STY']
    assert (tmp_path / 'STY.mol2').is_file()


@needs_obabel
def test_generated_mol2_carries_the_requested_reactive_atom_names(tmp_path):
    from htpolynet.core.coordinates import Coordinates
    spec = {'STY': {'smiles': STYRENE, 'reactive_atoms': {1: 'C1', 2: 'C2'}}}
    materialize_smiles_inputs(spec, inputs_dir=str(tmp_path))
    names = set(Coordinates.read_mol2(str(tmp_path / 'STY.mol2')).A['atomName'])
    assert {'C1', 'C2'} <= names


@needs_obabel
def test_generated_mol2_has_hydrogens_added(tmp_path):
    """Styrene is C8H8; RDKit AddHs must run before embedding."""
    from htpolynet.core.coordinates import Coordinates
    spec = {'STY': {'smiles': STYRENE, 'reactive_atoms': {1: 'C1', 2: 'C2'}}}
    materialize_smiles_inputs(spec, inputs_dir=str(tmp_path))
    assert Coordinates.read_mol2(str(tmp_path / 'STY.mol2')).A.shape[0] == 16


@needs_obabel
def test_only_smiles_bearing_constituents_are_reported_as_generated(tmp_path):
    spec = {
        'STY': {'smiles': STYRENE, 'reactive_atoms': {1: 'C1', 2: 'C2'}},
        'PAC': {'count': 50},
    }
    assert materialize_smiles_inputs(spec, inputs_dir=str(tmp_path)) == ['STY']


@needs_obabel
def test_unparseable_smiles_raises(tmp_path):
    spec = {'BAD': {'smiles': '[CH2:1]=[[[', 'reactive_atoms': {1: 'C1'}}}
    with pytest.raises(Exception):
        materialize_smiles_inputs(spec, inputs_dir=str(tmp_path))


# --- a kept mol2 that no longer matches its config ---------------------------

@needs_obabel
def test_a_kept_mol2_from_a_different_smiles_is_warned_about(tmp_path, caplog):
    """Examples 0 and 1 both name a molecule STY: real styrene in 0, activated in 1."""
    real = {'smiles': 'c1ccccc1[CH:1]=[CH2:2]', 'reactive_atoms': {1: 'C1', 2: 'C2'}}
    activated = dict(real, activate=[[1, 2]])
    assert materialize_smiles_inputs({'STY': real}, inputs_dir=str(tmp_path)) == ['STY']
    with caplog.at_level('WARNING'):
        assert materialize_smiles_inputs({'STY': activated}, inputs_dir=str(tmp_path)) == []
    assert any('was generated from' in r.message for r in caplog.records)


@needs_obabel
def test_a_kept_mol2_from_the_same_smiles_is_not_warned_about(tmp_path, caplog):
    spec = {'smiles': 'c1ccccc1[CH:1]=[CH2:2]', 'reactive_atoms': {1: 'C1', 2: 'C2'}, 'activate': [[1, 2]]}
    materialize_smiles_inputs({'STY': spec}, inputs_dir=str(tmp_path))
    with caplog.at_level('WARNING'):
        materialize_smiles_inputs({'STY': spec}, inputs_dir=str(tmp_path))
    assert not [r for r in caplog.records if r.levelname == 'WARNING']
