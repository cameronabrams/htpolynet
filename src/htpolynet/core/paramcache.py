"""Provenance records for cached GAFF parameterizations.

A parameterization checked into the user library is identified by molecule
name alone, but the numbers it contains depend on more than the name: the
antechamber charge method, the net charge passed as ``-nc``, and the atom-type
set.  Without a record of those, a configuration asking for
``charge_method: bcc`` silently reuses a ``gas`` parameterization checked in
months earlier under the same residue name, and the only trace is a log line
saying the cache was used.  The resulting system carries one charge method on
some residues and another on the rest, and nothing in the output says so.

This module writes a small JSON sidecar next to a parameterization recording
what produced it, and compares that record against what the current run asks
for.  A mismatch is a cache miss.  An *absent* record -- any library built
before sidecars existed -- is a cache hit with a warning, so that existing
libraries keep working rather than silently re-parameterizing wholesale.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import hashlib
import json
import logging
import os

logger = logging.getLogger(__name__)

CACHE_KEY_EXT = 'parm'
"""Filename extension of the sidecar record, alongside gro/top/itp/tpx/grx."""

_KEY_FIELDS = (
    ('charge_method', 'charge method'),
    ('net_charge',    'net charge'),
    ('atom_type',     'atom-type set'),
)
"""Fields that determine the parameters, with human-readable labels."""

FRCMOD_FIELD = 'frcmod'
"""Record field holding a digest of the user frcmod files loaded after parmchk2's.

Written only when a molecule has any, so records of molecules without one are
unchanged; its absence means "none" when compared.
"""


def frcmod_digest(frcmods):
    """Returns a digest identifying a list of user frcmod files by content.

    Args:
        frcmods (list): (label, text) tuples, in load order

    Returns:
        str or None: hex digest, or None if the list is empty
    """
    if not frcmods:
        return None
    h = hashlib.sha256()
    for label, text in frcmods:
        h.update(label.encode('utf-8') + b'\0' + text.encode('utf-8') + b'\0')
    return h.hexdigest()[:16]


def build_key(ambertools=None):
    """Returns the provenance record describing a parameterization run under
    the given AmberTools directives.

    Defaults come from :data:`htpolynet.external.ambertools.AMBERTOOLS_DEFAULTS`
    so that the record cannot drift from what GAFFParameterize actually does.

    Args:
        ambertools (dict, optional): ambertools configuration directives

    Returns:
        dict: provenance record
    """
    from ..external.ambertools import AMBERTOOLS_DEFAULTS
    ambertools = ambertools or {}
    key = {field: ambertools.get(field, AMBERTOOLS_DEFAULTS[field])
           for field, _ in _KEY_FIELDS}
    digest = frcmod_digest(ambertools.get(FRCMOD_FIELD))
    if digest:
        key[FRCMOD_FIELD] = digest
    return key


def key_filename(prefix):
    """Returns the sidecar filename for the given file basename.

    Args:
        prefix (str): file basename

    Returns:
        str: sidecar filename
    """
    return f'{prefix}.{CACHE_KEY_EXT}'


def write_key(prefix, key):
    """Writes a provenance record beside the parameterization it describes.

    Args:
        prefix (str): file basename
        key (dict): provenance record
    """
    with open(key_filename(prefix), 'w') as f:
        json.dump(key, f, indent=2, sort_keys=True)
        f.write('\n')


def read_key(prefix):
    """Reads the provenance record for the given file basename.

    A record that is missing, unreadable, or malformed is reported as absent
    rather than raising, since an unusable record must not be able to abort a
    build that would otherwise have proceeded.

    Args:
        prefix (str): file basename

    Returns:
        dict or None: provenance record, or None if absent or unusable
    """
    fn = key_filename(prefix)
    if not os.path.exists(fn):
        return None
    try:
        with open(fn) as f:
            key = json.load(f)
    except (OSError, ValueError) as m:
        logger.warning(f'Could not read parameterization record {fn} ({m}); treating it as absent')
        return None
    if not isinstance(key, dict):
        logger.warning(f'Parameterization record {fn} is not a JSON object; treating it as absent')
        return None
    return key


def describe_mismatch(stored, requested):
    """Returns human-readable descriptions of every field on which a stored
    record disagrees with what the current run requests.

    A field absent from the stored record is skipped rather than counted as a
    mismatch, so that adding a field here does not invalidate every record
    written by an earlier version.

    Args:
        stored (dict or None): record read from the cache
        requested (dict): record describing what this run asks for

    Returns:
        list: descriptions of the differing fields; empty if the record agrees
            or is absent
    """
    if not stored:
        return []
    diffs = []
    for field, label in _KEY_FIELDS:
        if field not in stored:
            continue
        if stored[field] != requested.get(field):
            diffs.append(f'{label} {stored[field]!r} cached, {requested.get(field)!r} requested')
    # unlike the fields above, an absent frcmod digest is a statement: none was loaded
    if stored.get(FRCMOD_FIELD) != requested.get(FRCMOD_FIELD):
        diffs.append(f'user frcmod {stored.get(FRCMOD_FIELD) or "none"} cached, '
                     f'{requested.get(FRCMOD_FIELD) or "none"} requested')
    return diffs
