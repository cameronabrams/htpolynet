.. _badcy_configuration:

The Configuration File
----------------------

The complete ``6-cyanate-ester.yaml`` from the depot:

.. literalinclude:: ../../../../src/htpolynet/resources/example_depot/6-cyanate-ester.yaml
   :language: yaml

Most blocks (``Title``, ``gromacs``, ``ambertools``, ``densification``, ``precure``,
``postcure``) follow the same conventions as the earlier tutorials.  The tutorial-
specific content is in ``constituents`` and ``reactions``, covered on the
:ref:`monomers <badcy_monomers>` and :ref:`reactions <badcy_reactions>` pages, and in
the ``ring_cure`` block, which takes the place of ``CURE``.

``ring_cure``
^^^^^^^^^^^^^

.. code-block:: yaml

   ring_cure:
     search_radius: 0.5
     radial_increment: 0.05
     desired_conversion: 0.97
     max_iterations: 150
     closure:
       nstages: 6

* ``search_radius`` (nm) -- the longest new bond a candidate ring may need at the
  start of an iteration.  ``radial_increment`` widens it when an iteration finds too
  few rings.
* ``desired_conversion`` -- the fraction of cyanate **groups** consumed, three per
  ring.  0.97 is in the range a fully post-cured BADCy reaches.
* ``closure.nstages`` -- how many restraint stages pull each candidate triple shut.

Everything else is left at its default.  Two defaults are worth knowing about:

* ``reject_threaded_rings`` and ``reject_new_threading`` are both on.  The first
  declines a ring that would close around an existing bond; the second declines a
  batch member whose dragged monomer ends up with a bond threaded through an existing
  ring.  Threading is permanent -- escaping it would need a bond to break -- so it is
  prevented rather than repaired.  ``htpolynet piercings`` counts any that remain in a
  finished build.
* ``relax``, ``equilibrate`` and ``settle`` give the ring cure the same breathing cycle
  ``CURE`` has: a 600 K relaxation and a 300 K, 1 bar re-equilibration after each
  iteration, and a minimization plus a short 300 K NVT run once the cure ends.

Every ``ring_cure`` setting, with its default and a description, is declared in the ``ring_cure`` section of ``src/htpolynet/schema/base.yaml``.

``CURE``, ``postcure_repair``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Neither appears.  The one cure reaction is ring-closing, so ``CURE`` has nothing to do
and says so in the log.  And because a cyanate group that never closes a ring simply
stays -O-C#N, there is nothing for a :ref:`postcure repair stage <postcure_repair>` to
fix.

What expansion produces
^^^^^^^^^^^^^^^^^^^^^^^

Two templates: the monomer ``BDC`` and the trimer ``BDC3``.  At setup the log reports::

   INFO> 2 molecules detected in 6-cyanate-ester.yaml
   INFO>                       explicit: 2
   INFO>     implied by stereochemistry: 0
   INFO>            implied by symmetry: 0
   ...
   INFO> Generated 2 molecule templates
   INFO> Initial composition is BDC 360

Next is :ref:`actually running the build <badcy_run>`.
