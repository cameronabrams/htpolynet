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

   **Every example uses AM1-BCC partial charges** (``charge_method: bcc``)
   with GAFF (``gaff``, i.e. GAFF1).  GAFF's Lennard-Jones and bonded
   parameters were fitted alongside RESP charges, and AM1-BCC exists to
   reproduce RESP cheaply, so this is the pairing GAFF is meant to be used
   with.  It costs a semi-empirical (``sqm``) calculation for every monomer
   and every reaction template before the build starts: 18 minutes on 16 CPU
   cores for example 6's two molecules, and more for examples with more or
   larger templates.  ``sqm`` does not use a GPU.

   Through htpolynet 2.17 the examples used Gasteiger charges
   (``charge_method: gas``), which are fast but are not designed to reproduce
   RESP.  On example 6 the switch raised the density at 300 K by about 1.1 %
   and Young's modulus by about 20 %, and did not change *T*:sub:`g`.  **The
   results recorded on these pages, and their timings, are still from the
   Gasteiger builds** and will be re-recorded; expect your own densities and
   moduli to come out higher.  ``gas`` remains ``htpolynet``'s default for a
   configuration that does not set ``charge_method``.  A parameterization
   already in your library is reused only if it was made with the same
   ``charge_method``; see :ref:`the ambertools directive <configuration_run>`.

.. note::

   **The** *T*:sub:`g` **to quote is the two-line crossing**, the ``Tg =``
   line of the ``plots post`` log, with the per-build spread on the line
   after it when there is more than one build.  ``plots post`` also logs a
   whole-curve hyperbola fit, as a diagnostic only: its *T*:sub:`g` depends
   on whether density, its logarithm or the specific volume is fitted, and
   drifts as the ladder's top changes.  Some *T*:sub:`g` values on these
   pages come from that fit and will be replaced when the pages are
   re-recorded.

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
