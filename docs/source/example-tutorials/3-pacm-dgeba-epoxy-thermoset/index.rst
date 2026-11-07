.. _pde_tutorial:

DGEBA-PACM Epoxy Thermoset
==========================

This tutorial covers a step-growth cure between a diepoxide (DGEBA) and a
diamine (PACM): primary amines open epoxide rings to form secondary
amines, which can then react with a second epoxide to form tertiary
amines.  Compared to :ref:`example 2 <bgs_tutorial>`, the new wrinkles are:

* **Two sequential cure reactions**.  The product of the first reaction
  (a secondary amine) is itself a reactant in the second reaction
  (forming a tertiary amine).  ``htpolynet`` chain-expands these into a
  full set of templates that cover both bond environments.
* **Symmetry-equivalent atoms**.  DGEBA has two chemically equivalent
  reactive carbons (and two equivalent oxiranes); PACM has two
  equivalent reactive nitrogens.  The ``symmetry_equivalent_atoms``
  block tells ``htpolynet`` to expand the user-supplied reactions into
  the full set of symmetry-equivalent products.
* **A cap reaction that re-forms a ring**.  Any leftover unreacted
  oxirane gets its hydroxyl-to-carbon bond re-formed as an epoxide
  three-ring at the end of the build.

.. note::

   The configuration was enlarged in November 2026 for showcasing the analysis: it
   doubles the box side, 8 times the monomers of the earlier one.  Recorded output, timings and
   measured values on these pages are from the earlier, smaller builds until they are
   re-recorded from builds at the new size.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   monomers
   reactions
   configuration
   run
   results
   postsim
