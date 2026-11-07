.. _htpb_introduction:

Introduction
------------

Set up a clean working directory and pull the example YAML:

.. code-block:: console

   $ mkdir my_htpb_ipdi
   $ cd my_htpb_ipdi
   $ htpolynet fetch-example 5
   Fetched 5-htpb-ipdi.yaml  (run with: htpolynet run 5-htpb-ipdi.yaml)
   $ ls
   5-htpb-ipdi.yaml

Self-contained YAML as in the earlier examples — all four small
constituents are generated from SMILES, and the two long-chain HTPB
monomers (``DHT`` and ``THT``) are *assembled* from those small units
by the param- and build-stage reactions described below.

What's chemically new here is the two-tier construction:

* **Stage 1 — chain assembly.**  At param-stage, small reactions
  combine 1-butene (``OB``), trans-2-butene (``TB``), and
  1-hydroxy-trans-2-butene (``TBO``) into a handful of dimer / trimer
  / quad sub-units (``A2``, ``AO``, ``OBT``, ``OB3``, ``A3``,
  ``A4``).  Each of these is a real GAFF-parameterized molecule.  At
  build-stage, those sub-units are stitched into the two HTPB chains:
  ``DHT`` is linear (three 18-residue ``A18`` segments capped at both
  ends with ``TBO``); ``THT`` is branched (four segments meeting at a
  central ``A4`` node, capped at all four arms).  Build-stage
  reactions inherit bonded parameters from the param-stage templates
  — no new ``antechamber`` runs at build-stage.
* **Stage 2 — urethane cure.**  Each chain's terminal hydroxyl
  (``TBO.O1``) reacts with one of IPDI's two non-equivalent
  isocyanate carbons (``C1`` or ``C2`` of the active-form ``IPD``
  monomer; see the :ref:`monomers page <htpb_monomers>` for what
  "active form" means here).  The result is a urethane bond
  ``-O-C(=O)-N(H)-`` connecting the HTPB end to an IPDI residue.  IPD
  is bifunctional, so each IPD bridges two HTPB hydroxyls (potentially
  on different chains) and the network grows.

What you'll see in the build
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The build is long — both because the system is large (~111,000 atoms
final) and because the cure step-growth runs out of nearby reactive
pairs much faster than in the bulk-cure examples.  On a typical
single-CPU + single-GPU run, expect:

* ~5-10 minutes of setup (most of which is the 30+ template
  parameterizations during chain assembly);
* ~1h50m of densification + precure anneal (the long HTPB chains need
  time to pack from the low initial density, and densification runs
  until the density settles rather than for a fixed count);
* ~1h45m of cure (9 iterations to 95 % conversion, with the last two
  taking ~25 minutes each as the remaining hydroxyl / isocyanate pairs
  grow scarce);
* ~7 minutes of postcure anneal + equilibration.

The :ref:`run page <htpb_run>` gives a per-iteration breakdown and
shows where ``min_bonds_per_iteration`` and ``late_threshold`` would
matter for tuning.

The remaining pages walk through the monomer set, the
param/build/cure reaction pipeline, the YAML in full, and what to
look for in the diagnostic log and plots.
