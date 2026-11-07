.. _htpb_tutorial:

HTPB / IPDI Polyurethane Thermoset
==================================

This tutorial covers a **polyurethane thermoset** built from
hydroxyl-terminated polybutadiene (HTPB) and isophorone diisocyanate
(IPDI).  Two features are unique to this example among the depot
examples:

* **Long-chain monomers assembled at runtime.**  HTPB is a long
  polybutadiene chain (typical molecular weight 2-3 kDa) terminated at
  both ends by hydroxyl groups.  Rather than write out a single
  monstrous SMILES for it, ``htpolynet``'s param-stage and build-stage
  reactions chain together small parameterized sub-units (a 2-butene
  dimer, an 18-residue chain, a hydroxyl end-cap) into the assembled
  monomer.  Two variants are built: ``DHT`` (linear, head + tail end-
  caps) and ``THT`` (branched, three end-caps).  This is the only
  depot example that exercises the param/build-stage reaction
  pipeline at scale.

* **Urethane cure chemistry.**  Each HTPB hydroxyl bonds to one of
  IPDI's two non-equivalent isocyanate carbons to form a urethane
  linkage (``HTPB-O-C(=O)-NH-IPDI``).  The system is end-functional —
  cure only happens at the HTPB chain ends — so the network grows by
  IPDI joining HTPB chain ends rather than by attacks along the
  backbone.  The 2:1:1 stoichiometry (125 IPD : 50 DHT : 50 THT)
  ensures every IPD's two isocyanates can find HTPB hydroxyl partners.

This is the largest depot example by atom count (~111,000 atoms after
build) and by wall-clock time (cure alone takes roughly 10 hours on a
single CPU + GPU).  A full run is meaningful as an overnight job
rather than an interactive demo.

.. note::

   The configuration was enlarged in November 2026 for showcasing the analysis: it
   has 5 times the monomers of the earlier one.  Recorded output, timings and
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
