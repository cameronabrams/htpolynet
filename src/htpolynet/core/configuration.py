"""Manages reading and parsing of YAML configuration files.

Author: Cameron F. Abrams <cfa22@drexel.edu>

Configurations are validated against ``htpolynet/schema/base.yaml`` by
ycleptic, which also fills in every documented default.  Before that, this
module read the file and ``.get()``-ed the keys it recognized, so an
unrecognized section or a misspelled key was discarded in silence -- a
``desired_converson`` typo produced a default-value build and no message.

Keys whose *presence* changes behavior are declared ``required: false`` in the
schema so ycleptic leaves them absent rather than filling them in; see the
comments on ``densification.initial_boxsize``, ``reactive_atoms`` and
``rename_atoms`` there.
"""
import json
import logging
import os

from importlib.resources import as_file, files

import yaml

from ycleptic.errors import YclepticError
from ycleptic.yclept import Yclept

logger = logging.getLogger(__name__)


def schema_path():
    """Context manager yielding a real path to the packaged base config.

    Returns:
        contextlib.AbstractContextManager: yields a pathlib.Path to base.yaml
    """
    return as_file(files('htpolynet').joinpath('schema/base.yaml'))


class Configuration:
    """Pure data container for a build configuration file.

    Reads and validates a YAML or JSON file and exposes each top-level
    section as a named attribute.  Does NOT create Molecule or Reaction
    objects; that is the responsibility of Runtime.
    """

    def __init__(self):
        self.cfgfile = ''
        self.title = ''
        self.ncpu = os.cpu_count()
        self.constituents = {}       # {name: {count: N, ...}}
        self.reaction_specs = []     # raw reaction dicts
        self.initial_composition = []  # [{'molecule': name, 'count': N}]
        self.gromacs = {}
        self.densification = {}
        self.precure = {}
        self.postcure_repair = []    # list of postcure topology-repair specs
        self.postcure = {}
        self.cure = {}
        self.ring_cure = {}
        self.gaff = {}
        self.ambertools = {}
        self.resolve_type_discrepancies = []
        self.basedict = {}

    @classmethod
    def read(cls, filename):
        """Reads a JSON or YAML configuration file and returns a populated Configuration.

        Args:
            filename (str): path to configuration file

        Raises:
            Exception: if the file extension is not .json, .yaml, or .yml

        Returns:
            Configuration: populated configuration object
        """
        _, ext = os.path.splitext(filename)
        inst = cls()
        inst.cfgfile = filename
        if ext == '.json':
            # ycleptic reads YAML; round-trip JSON through a temporary file so
            # both formats get the same validation rather than only one
            with open(filename, 'r') as f:
                raw = json.load(f)
            inst.basedict = inst._validate(raw)
        elif ext in ('.yaml', '.yml'):
            with open(filename, 'r') as f:
                raw = yaml.safe_load(f)
            inst.basedict = inst._validate(raw)
        else:
            raise Exception(f'Unknown config file extension {ext}')
        inst._parse()
        return inst

    def _validate(self, raw):
        """Validates a raw config against the packaged schema and fills defaults.

        Args:
            raw (dict): the configuration as read from disk

        Raises:
            Exception: with ycleptic's message, if the configuration is invalid

        Returns:
            dict: the validated configuration, with defaults filled in
        """
        import tempfile
        with schema_path() as base:
            with tempfile.TemporaryDirectory() as td:
                userfile = os.path.join(td, 'user.yaml')
                with open(userfile, 'w') as f:
                    yaml.dump(raw, f, default_flow_style=False, sort_keys=False)
                try:
                    y = Yclept(basefile=str(base), userfile=userfile)
                except YclepticError as e:
                    raise Exception(f'{self.cfgfile}: {e}') from None
        return y['user']

    def _parse(self):
        """Populates named attributes from self.basedict."""
        self.title = self.basedict.get('Title', 'No title provided')
        self.ncpu = self.basedict.get('ncpu', os.cpu_count())
        self.constituents = self.basedict.get('constituents', {})
        self.reaction_specs = self.basedict.get('reactions', [])
        self.gromacs = self.basedict.get('gromacs', {})
        self.densification = self.basedict.get('densification', {})
        self.precure = self.basedict.get('precure', {})
        self.postcure_repair = self.basedict.get('postcure_repair', [])
        self.postcure = self.basedict.get('postcure', {})
        self.cure = self.basedict.get('CURE', {})
        self.ring_cure = self.basedict.get('ring_cure', {})
        self.gaff = self.basedict.get('GAFF', {})
        self.ambertools = self.basedict.get('ambertools', {})
        self.resolve_type_discrepancies = self.basedict.get('resolve_type_discrepancies', [])
        self.initial_composition = [
            {'molecule': mol, 'count': rec.get('count', 0)}
            for mol, rec in self.constituents.items()
        ]
