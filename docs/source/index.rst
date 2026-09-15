.. htpolynet documentation master file, created by
   sphinx-quickstart on Mon May 16 14:53:57 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#########
htpolynet
#########

|tests| |pypi| |conda-forge| |python| |license| |docs| |downloads| |doi|

.. |tests| image:: https://github.com/cameronabrams/htpolynet/actions/workflows/test.yml/badge.svg
   :target: https://github.com/cameronabrams/htpolynet/actions/workflows/test.yml
   :alt: Tests

.. |pypi| image:: https://img.shields.io/pypi/v/htpolynet.svg
   :target: https://pypi.org/project/htpolynet/
   :alt: PyPI

.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/htpolynet
   :target: https://anaconda.org/conda-forge/htpolynet
   :alt: conda-forge

.. |python| image:: https://img.shields.io/pypi/pyversions/htpolynet
   :target: https://pypi.org/project/htpolynet/
   :alt: Python versions

.. |license| image:: https://img.shields.io/pypi/l/htpolynet
   :target: https://github.com/cameronabrams/htpolynet/blob/main/LICENSE
   :alt: License: MIT

.. |docs| image:: https://readthedocs.org/projects/htpolynet/badge/?version=latest
   :target: https://htpolynet.readthedocs.io/en/latest/
   :alt: Documentation status

.. |downloads| image:: https://static.pepy.tech/badge/htpolynet
   :target: https://pepy.tech/projects/htpolynet
   :alt: PyPI downloads

.. |doi| image:: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22070252-blue
   :target: https://doi.org/10.5281/zenodo.22070252
   :alt: DOI

.. image:: _static/htpolynet-landing.jpg
   :alt: A rusty robot in a straw hat sits on a driftwood log on a tropical beach, humming as it knots a net of ball-and-stick molecules strung between posts. A hand-painted sign beside it reads "HTPOLYNET: Stranded today. Better polymers tomorrow."
   :align: center
   :width: 100%

htpolynet is a command-line tool for building atomic configurations of amorphous network polymers suitable for molecular dynamics (MD) simulations.
It uses the General Amber Force Field and produces output that can be simulated 
using Gromacs.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   user-guide/index
   example-tutorials/index
   references/index
   htpolynetpackage
   changelog

Citation
========

When using htpolynet in published work, please cite :cite:t:`Huang2023`, along with the main GAFF paper :cite:t:`Wang2004DevelopmentField` and the two main Gromacs papers :cite:t:`Berendsen1995` and :cite:t:`Abraham2015`.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
