"""Manages execution of AmberTools antechamber, parmchk2.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import hashlib
import logging
import os
import shutil

import parmed

from ..core import paramcache
from ..core.coordinates import Coordinates
from ..external.command import run

logger = logging.getLogger(__name__)

AMBERTOOLS_DEFAULTS = {
    'charge_method': 'bcc',
    'net_charge': 0,
    'atom_type': 'gaff',
}
"""Directives that determine the parameters GAFFParameterize produces.

These are the single source of truth for both the antechamber/parmchk2/tleap
invocations below and the provenance record written beside their output; see
:mod:`htpolynet.core.paramcache`.
"""


def leapsafe(name):
    """Returns a tleap-safe file stem for name.

    tleap reads a bare filename containing 'e' after a digit as a number, so
    names are hashed, as GAFFParameterize does for its own files.
    """
    return 'u' + hashlib.shake_128(name.encode('utf-8')).hexdigest(8).replace('e', 'x')


def GAFFParameterize(inputPrefix, outputPrefix, input_structure_format='mol2', ambertools=None):
    """Manages execution of antechamber, parmchk2, and tleap to generate GAFF parameters,
    then converts the result to Gromacs gro/top files via parmed.

    Args:
        inputPrefix (str): basename of input structure file
        outputPrefix (str): basename of output files
        input_structure_format (str): format of input structure file, defaults to 'mol2'; 'pdb' is other option
        ambertools (dict | None): ambertools configuration directives, defaults to None; an
            ``frcmod`` entry, a list of (label, text) tuples, names user frcmod files that
            tleap loads after parmchk2's, so their parameters take precedence

    Raises:
        parmed.exceptions.GromacsError: if parmed fails to convert tleap output
    """
    ambertools = ambertools or {}
    chargemethod = ambertools.get('charge_method', AMBERTOOLS_DEFAULTS['charge_method'])
    netcharge    = ambertools.get('net_charge',    AMBERTOOLS_DEFAULTS['net_charge'])
    atomtype     = ambertools.get('atom_type',     AMBERTOOLS_DEFAULTS['atom_type'])
    logger.info(f'AmberTools> generating GAFF parameters from {inputPrefix}.{input_structure_format}')

    structin  = f'{inputPrefix}.{input_structure_format}'
    mol2out   = f'{outputPrefix}.mol2'
    frcmodout = f'{outputPrefix}.frcmod'
    groOut    = f'{outputPrefix}.gro'
    topOut    = f'{outputPrefix}.top'
    itpOut    = f'{outputPrefix}.itp'

    # antechamber overwrites its input when input and output share the same path
    new_structin = structin
    if structin == mol2out:
        new_structin = f'{inputPrefix}-input.{input_structure_format}'
        logger.debug(f'AmberTools> Antechamber overwrites input {structin}; backing up to {new_structin}')
        shutil.copy(structin, new_structin)

    run(f'antechamber -j 4 -fi {input_structure_format} -fo mol2 -c {chargemethod}'
        f' -at {atomtype} -i {new_structin} -o {mol2out} -pf Y -nc {netcharge} -eq 1 -pl 10', quiet=False)
    logger.debug(f'AmberTools> Antechamber generated {mol2out}')
    run(f'parmchk2 -i {mol2out} -o {frcmodout} -f mol2 -s {atomtype}', quiet=False)
    user_frcmods = []
    for label, text in ambertools.get('frcmod') or []:
        fn = f'{leapsafe(label)}-user.frcmod'
        with open(fn, 'w') as f:
            f.write(text)
        user_frcmods.append(fn)
        logger.info(f'AmberTools> loading user frcmod for {label} into {outputPrefix}')

    # Antechamber ignores SUBSTRUCTURE records, so patch the antechamber output
    # mol2 with the original resName/resNum before passing it to tleap.
    if input_structure_format == 'mol2':
        orig    = Coordinates.read_mol2(new_structin)
        ac_out  = Coordinates.read_mol2(mol2out)
        ac_out.A[['resName', 'resNum']] = orig.A[['resName', 'resNum']]
        leap_coords = ac_out
    else:
        leap_coords = Coordinates.read_mol2(mol2out)

    # tleap rejects filenames containing 'e' (treats them as scientific notation),
    # so hash the prefix to get a safe temporary name.
    leapprefix = hashlib.shake_128(outputPrefix.encode('utf-8')).hexdigest(16).replace('e', 'x')
    logger.debug(f'Replacing "{outputPrefix}" with hash "{leapprefix}" for tleap input files.')
    leap_coords.write_mol2(f'{leapprefix}.mol2')
    shutil.copy(frcmodout, f'{leapprefix}.frcmod')

    with open(f'{inputPrefix}-tleap.in', 'w') as f:
        f.write('\n'.join([
            f'source leaprc.{atomtype}',
            # Load parmchk2's patches BEFORE checking the molecule, otherwise
            # `check` reports any GAFF-coverage gaps (e.g. h5-ce-n2 on cyanate-
            # ester dimers) as `Error!`, the run-wrapper's override needle
            # fires, and we abort even though tleap would have completed once
            # the frcmod was loaded.
            f'loadamberparams {leapprefix}.frcmod',
            # user files last, so they override both GAFF and parmchk2's guesses
            *[f'loadamberparams {fn}' for fn in user_frcmods],
            f'mymol = loadmol2 {leapprefix}.mol2',
            'check mymol',
            f'saveamberparm mymol {leapprefix}-tleap.top {leapprefix}-tleap.crd',
            'quit',
            '',
        ]))

    run(f'tleap -f {inputPrefix}-tleap.in', override=('Error!', 'Unspecified tleap error'))

    os.remove(f'{leapprefix}.frcmod')
    for fn in user_frcmods:
        os.remove(fn)
    shutil.move(f'{leapprefix}.mol2',      mol2out)
    shutil.move(f'{leapprefix}-tleap.top', f'{outputPrefix}-tleap.top')
    shutil.move(f'{leapprefix}-tleap.crd', f'{outputPrefix}-tleap.crd')

    try:
        file = parmed.load_file(f'{outputPrefix}-tleap.top', xyz=f'{outputPrefix}-tleap.crd')
        logger.debug(f'Writing {groOut}, {topOut}, and {itpOut}')
        file.save(groOut, overwrite=True)
        file.save(topOut, parameters=itpOut, overwrite=True)
    except Exception as m:
        logger.error('Unspecified parmed error')
        raise parmed.exceptions.GromacsError(m) from m

    # Record what produced these parameters, so that a later run reusing them
    # from the library can tell whether they answer the question it is asking.
    paramcache.write_key(outputPrefix, paramcache.build_key(ambertools))

