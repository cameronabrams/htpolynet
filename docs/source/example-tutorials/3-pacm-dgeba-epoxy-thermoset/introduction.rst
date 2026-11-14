.. _pde_introduction:

Introduction
------------

Set up a clean working directory and pull the example YAML:

.. code-block:: console

   $ mkdir my_dgeba_pacm
   $ cd my_dgeba_pacm
   $ htpolynet fetch-example 3
   Fetched 3-pacm-dgeba-epoxy-thermoset.yaml  (run with: htpolynet run 3-pacm-dgeba-epoxy-thermoset.yaml)
   $ ls
   3-pacm-dgeba-epoxy-thermoset.yaml

Self-contained YAML as in the earlier examples — both monomers are
generated from SMILES with atom-mapping tokens at the start of the run.

What's chemically new here is the epoxy + amine step-growth mechanism:

* A primary amine (``-NH₂``) opens an epoxide ring (``C–O–C``) to form a
  secondary amine and a pendant hydroxyl.  The amine donates one H to
  the oxirane oxygen.
* The resulting secondary amine can react again with another epoxide,
  forming a tertiary amine.  This second reaction is typically slower,
  which we encode with a per-reaction ``probability`` in the YAML.
* Each PACM molecule can therefore tether up to four DGEBAs through its
  two nitrogens; each DGEBA can tether up to two amines through its two
  reactive carbons.

DGEBA is written as it really is, and its ``activate`` key has
``htpolynet`` work with its active form, in which both oxiranes are
already "opened" into a secondary alcohol + sacrificial methyl
configuration, so each reactive carbon carries a sacrificial H that gets
removed when a new C–N bond forms.  The :ref:`monomer page <pde_monomers>`
walks through this in detail.
