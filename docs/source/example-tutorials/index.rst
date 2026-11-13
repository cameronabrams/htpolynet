.. _example_tutorials:

Example Tutorials
=================

We present here a series of tutorials to help illustrate usage of ``htpolynet``.  Each tutorial walks through one polymerization recipe end-to-end: setting up monomer structures and reactions, running the build with ``htpolynet run``, and then performing post-build MD simulations and analyses with ``htpolynet postsim``, ``plots``, and ``analyze``.  Example ``0`` (liquid styrene) does no polymerization and is the simplest possible build; the remaining examples add cure and cap reactions of increasing complexity, with tutorial ``6`` (BADCy cyanate-ester thermoset) curing by cyclotrimerization, a three-body reaction that closes a ring.  Each tutorial number matches the example number returned by ``htpolynet fetch-example``.

.. note::

   **IMPORTANT DISCLAIMER**

   These are *not* production-level builds.  The system sizes are *way too small* and the equilibration and post-build simulation times are *way too short*.  You as a user are responsible for conducting the appropriate finite-size-effect tests and equilibration tests needed to guarantee robustness of your simulations.

.. note::

   **Every example gates its densification on measured density.**  Each one's
   final densification NPT stage carries a ``converge`` block, so instead of
   running for a fixed number of steps it repeats until the density has
   settled -- or reports that it has not.  Wall times below therefore vary
   with how quickly a given box compacts, and a densification stage may run
   two or three times the duration its ``ps`` value suggests.  See
   :ref:`the converge subdirective <configuration_run>` for what the criterion
   is and how to loosen or tighten it.

.. note::

   **Every example uses Gasteiger partial charges, which is a speed choice, not
   a modeling one.**  Each example's ``ambertools`` block sets
   ``charge_method: gas``, and none sets ``atom_type``, so all of them run
   GAFF (``gaff``, i.e. GAFF1) with Gasteiger charges.  That is also what
   ``htpolynet`` does for any configuration that does not say otherwise:
   ``gas`` is the default.

   GAFF's Lennard-Jones and bonded parameters were fitted alongside RESP
   charges, and AM1-BCC (``charge_method: bcc``) exists to reproduce RESP
   cheaply; that is the pairing GAFF is meant to be used with.  Gasteiger
   charges are not designed to reproduce RESP, so GAFF with Gasteiger is not a
   validated combination.  The examples use it because AM1-BCC runs a
   semi-empirical (``sqm``) calculation for every monomer and every reaction
   template, which makes a tutorial build much slower.

   So the densities, glass transitions and moduli these tutorials report are
   the ones their ``gas`` builds give, and are not meant to be compared with
   experiment.  For numbers you intend to publish or compare with measured
   data, set ``charge_method: bcc`` (see :ref:`the ambertools directive
   <configuration_run>`).  Note that a parameterization already in your
   library is reused only if it was made with the same ``charge_method``; see
   the same section for what happens to older library entries.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   0-liquid-styrene/index
   1-polystyrene/index
   2-bisgma-styrene-thermoset/index
   3-pacm-dgeba-epoxy-thermoset/index
   4-dfda-fde-epoxy-thermoset/index
   5-htpb-ipdi/index
   6-cyanate-ester/index
