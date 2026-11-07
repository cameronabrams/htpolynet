.. _dfe_tutorial:

DFA-FDE Furan-Based Epoxy Thermoset
===================================

This tutorial covers an epoxy-amine cure using **furan-derived**
monomers: FDE, a furfuryl diepoxide, and DFA, a difurfuryl diamine.
The cure chemistry is the same step-growth epoxy-amine reaction
covered in :ref:`example 3 <pde_tutorial>` — the new content here is:

* **Furan heterocycles in the monomers.**  Both species carry one or
  two aromatic furan rings, which are an attractive "green" feedstock
  family for thermoset polymers (furans can be derived from sugars).
* **A bisphenol-A-free, fused-heterocycle alternative.**  FDE plays
  the structural role DGEBA played in tutorial 3; DFA plays the
  structural role of PACM.  The reactive atoms and atom-mapping
  scheme on FDE are deliberately identical to those on DGE, so the
  reactions block reads almost the same.

Historically this example was painful because openbabel's mol2 writer
produced inconsistent atom numbering for fused heterocycles, forcing a
PDB-based workflow with manual atom renaming.  The SMILES atom-mapping
path now in ``htpolynet`` (``[CH:1]``, ``[NH2:1]`` etc.) sidesteps that
entirely — the same self-contained-YAML approach used by the earlier
examples works cleanly here too.

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
