"""

.. module:: test_config_docs_coverage
   :synopsis: every configuration section a user can write is documented for them

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import unittest

import yaml

from htpolynet.core.configuration import schema_path

# conftest chdirs each test into its own directory, so this has to be resolved from the
# test file rather than from the working directory
_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DOC = os.path.join(_REPO, 'docs', 'source', 'user-guide', 'configs', 'configs-for-run.rst')

UNDOCUMENTED = {
    # Known debt, named so it cannot grow quietly.  Each of these is a real section a
    # user can put in a config and get no guidance about.
    'ncpu',
    'ring_cure',         # the cyclotrimerization cure; ~15 settings, several added in
                         # 2.11.x and still moving
    'postcure_repair',
}
"""Sections the config guide does not cover yet.

This list exists so that a *new* undocumented section fails the suite while the existing
debt does not.  `docs/source/htpolynetpackage.rst` drifted for months in exactly this way
-- it autodoc'd a module that had ceased to exist -- and nothing noticed, because nothing
compared the document against the code.
"""


class TestEverySectionIsDocumented(unittest.TestCase):
    def sections(self):
        with schema_path() as p:
            return [a['name'] for a in yaml.safe_load(open(p))['attributes']]

    def test_no_new_section_is_undocumented(self):
        doc = open(DOC).read()
        missing = {n for n in self.sections() if f'``{n}``' not in doc}
        self.assertEqual(missing - UNDOCUMENTED, set(),
                         'a config section was added without documenting it in ' + DOC)

    def test_the_debt_list_does_not_outlive_the_debt(self):
        # when one of these is finally written up, it must come off the list
        doc = open(DOC).read()
        stale = {n for n in UNDOCUMENTED if f'``{n}``' in doc}
        self.assertEqual(stale, set(),
                         'these are documented now and should be removed from UNDOCUMENTED')

    def test_the_debt_list_names_only_real_sections(self):
        self.assertEqual(UNDOCUMENTED - set(self.sections()), set())
