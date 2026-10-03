.. _postcure_repair:

Postcure topology repair
========================

The cure and cap stages in ``htpolynet`` are *monotonic*: every event
they perform forms a new bond and (optionally) deletes a sacrificial
hydrogen.  They cannot break bonds, delete or insert heavy atoms,
re-tag atoms into a new residue, or otherwise re-topologize the system
in any non-additive way.

For some chemistries that's a problem.  A topological model of a
network thermoset — where you describe the *cured-network connectivity*
rather than the *bond-forming mechanism* — will in general leave
artefacts at finite cure conversion that don't exist in the real
material.  The motivating case was the BADCy cyanate-ester example
as it shipped before 2.15.0: a topological A2+B3 cure model (bisphenol A
plus a pre-formed triazine) that avoided a 3-way ring closure during
cure, but at less than full conversion left free
phenolic ``BPA-OH`` groups and bare triazine C-H positions — species
that real undercured BADCy doesn't carry.  Real undercured BADCy has
``-O-C#N`` end-groups instead, because the cyanate that didn't
cyclotrimerize just *stays* as an intact cyanate.

.. note::

   Since 2.15.0 :ref:`tutorial 6 <badcy_tutorial>` builds BADCy by real
   cyclotrimerization with a ``ring_cure`` block, which leaves unreacted
   cyanate groups intact by construction and needs no repair.  The repair
   stage remains for builds made with the stand-in model, and the
   ``triazine_to_cyanate_cap`` configuration shown below is the one that
   example used.

The postcure repair stage gives ``htpolynet`` an escape valve for
exactly this kind of problem.  After cure (and capping, if any)
finishes, but before postcure MD, ``htpolynet`` calls into a *repair
driver* that may perform arbitrary topology surgery — sever bonds,
delete atoms, relocate atoms between residues, re-template affected
linkages — to convert the cured topology into a chemically realistic
final state.

Architecture
------------

The repair stage lives in :mod:`htpolynet.repair`, organized as:

* **A dispatcher** (``repair/__init__.py:run_repair``).  Iterates the
  ``postcure_repair`` config block and routes each spec to the
  driver named by its ``type`` field.

* **Surgery primitives** (``repair/topology_surgery.py``).  Operations
  on a ``TopoCoord`` that the cure machinery does not expose:
  ``delete_bonds`` with cascading angle / dihedral / 1-4-pair
  cleanup, ``set_atom_attributes``, ``reassign_residue``,
  ``add_bonds_with_template`` (a wrapper around ``make_bonds`` +
  ``map_from_templates`` plus an int-dtype rescue for the atom-index
  columns), and a few smaller helpers.

* **Concrete drivers** — one Python module per ``type:`` value.  Each
  driver receives the current ``TopoCoord``, the molecule template
  dictionary, its spec dict, and the full reaction list, and is free
  to do whatever the chemistry calls for.  The first shipped driver
  is ``triazine_to_cyanate_cap`` (in ``repair/cyanate_cap.py``),
  detailed below.

* **A new** ``reaction_stage.repair`` **enum value**.  Lets
  repair-stage reactions in the YAML ride the same setup-time
  parameterization path that produces cure-stage linked-product
  templates.  At runtime these reactions are never *executed* like
  cure / cap reactions would be; their sole purpose is to define a
  parameterized linked-product template that a repair driver can
  splice into the system for every new bond it forms.

* **Runtime integration** in :class:`htpolynet.core.runtime.Runtime`.
  ``do_repair()`` is a new stage hooked into ``do_workflow`` between
  ``do_cure()`` and ``do_postcure()``.  When the YAML carries a
  non-empty ``postcure_repair`` block, the runtime creates a
  ``systems/repair/`` working directory, invokes the dispatcher,
  writes ``repaired.{gro,top,tpx,grx}``, and runs a steepest-descent
  minimization plus a short NVT settle on the modified topology
  before the postcure MD ensemble takes over.  The minimize +
  short-NVT pair absorbs any LJ clashes from physically relocating
  atoms during the surgery.

Configuration
-------------

Postcure repair is configured via a new top-level ``postcure_repair``
list:

.. code-block:: yaml

   postcure_repair:
     - type: <driver_name>
       <driver-specific-fields>
     - type: <other_driver_name>
       <other-driver-specific-fields>

The list-of-dicts shape lets multiple drivers run in sequence; the
runtime dispatches them in order.  Each entry must carry a
``type:`` key naming a registered driver; the remaining fields are
driver-specific.

Reaction templates for repair drivers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most repair drivers need parameterized templates for the bonds they
form.  Those templates are declared by adding ``stage: repair``
entries to the existing ``reactions:`` block:

.. code-block:: yaml

   reactions:
     - name: my_repair_template
       stage: repair
       reactants: {1: SomeBridgeResidue, 2: SomeCapResidue}
       product: SomeBridge~Atom-Atom~SomeCap
       atoms:
         A: {reactant: 1, resid: 1, atom: ..., z: 1}
         B: {reactant: 2, resid: 1, atom: ..., z: 1}
       bonds:
         - atoms: [A, B]
           order: 1

The ``stage: repair`` value tells the runtime to:

1. Include this reaction in setup-time parameterization (the
   ``product`` name shows up as a parameterized linked-product
   molecule).
2. Run the standard symmetry-expansion machinery (so a single
   reaction can generate as many template variants as the
   ``symmetry_equivalent_atoms`` declarations call for).
3. Skip this reaction in the actual cure / cap iteration loops.

The driver then looks up the parameterized template by its
``product`` name and uses it to splice atom types, charges, and
bonded interactions into the system at surgery time.

The ``triazine_to_cyanate_cap`` driver
--------------------------------------

The first concrete driver, written for the BADCy example, dismantles
incomplete triazine crosslinks into independent -C#N caps.  It is
specific to triazine + bisphenol A2+B3 chemistry but illustrates the
architecture cleanly.

Spec fields
^^^^^^^^^^^

.. code-block:: yaml

   postcure_repair:
     - type: triazine_to_cyanate_cap
       crosslinker:
         residue: TAZ
         ring_carbon_atoms: [C1, C2, C3]
         ring_nitrogen_atoms: [N1, N2, N3]
         full_bond_count: 3
       bridge:
         residue: BPA
         reactive_oxygen_atoms: [O1, O2]
       cap_residue: CYN
       cap_template: BPA~O1-C1~CYN
       cap_search_radius: 0.6   # nm
       cap_min_clearance: 0.15  # nm

* ``crosslinker`` — names the trifunctional residue to be inspected.
  ``ring_carbon_atoms`` / ``ring_nitrogen_atoms`` give the
  atom-names of the six ring atoms in traversal order (C1-N1-C2-N2-
  C3-N3-C1).  ``full_bond_count`` is the count of bonded bridge-Os
  above which the ring is considered "complete" and left alone;
  rings with fewer than this many bonded bridges are dismantled.
* ``bridge`` — names the difunctional residue and which atoms are
  the reactive sites.  Free bridge sites (atoms still carrying their
  sacrificial H after cure) are the recipients of donated free-cap
  fragments.
* ``cap_residue`` / ``cap_template`` — the residue name to assign to
  each newly-formed cap, and the parameterization template whose atom
  types / charges / bonded interactions get spliced in for each
  bridge-O-to-cap-C bond.  ``cap_template`` is the ``product:`` field
  of a ``stage: repair`` reaction declared elsewhere in the YAML.
* ``cap_min_clearance`` — how much room (nm) a transferred cap tries
  to leave itself.  The cap is first placed along the bridge-O's old
  O-H vector, as it always was; if that direction is already occupied
  the driver searches for the clearest direction, holding both bond
  lengths fixed so only the orientation moves, and stops as soon as it
  reaches this clearance.  Default 0.15 nm, which is about what a
  steepest-descent minimization absorbs.  It is a demanding default,
  and deliberately so: at the heavy-atom density of a cured thermoset
  it sits well above the room a typical site has along the O-H vector,
  so the direction search runs for most caps rather than
  rubber-stamping them.  Measured on 14 independent cured
  bisphenol-A-cyanate-ester boxes at about 53 heavy atoms per
  nm\ :sup:`3` (1955 caps in total), the median clearance along the
  O-H vector was 0.120 nm — below the target — and 37 % of caps would
  have been placed within 0.10 nm of a neighbour had that direction
  been taken blind.  That is the intent, and it is not a symptom: a
  run reporting that most caps needed a search is the default working.

  Demanding is not the same as unreachable.  In those same 14 boxes
  the search found the target for all but 1 cap of 1955, so a run
  reporting more than an occasional ``n_below_target`` is saying
  something about the box, not about the threshold.

  The clearance is measured only against atoms the cap is *not* bonded
  through — the bridge oxygen itself and its aryl carbon are excluded,
  because their distance to the cap is set by bond geometry rather
  than by how crowded the site is.  Leaving either in makes the number
  stop varying: with the oxygen in it is exactly the O-C bond length
  for every candidate direction, and with the aryl carbon in it is
  bounded by the C-O-C geometry, so a comfortably placed cap reports a
  constant instead of a measurement.  The search is restricted to
  C-O-C angles between 90° and 150° instead, which is what actually
  keeps the cap off the ring it hangs from — and off the linear
  geometry that pure clearance would otherwise drift toward, since
  antiparallel maximizes distance from the rest of the molecule.

  When no direction reaches the target, the driver says so and names
  the worst offenders; it also reports separately how many preferred
  directions were abandoned because of the angle window rather than
  because of crowding, so the two do not get read as each other.  See
  :ref:`what repair reports <postcure_repair_reporting>`.
* ``cap_search_radius`` — radius (nm) within which a free-cap
  fragment looks for an unreacted bridge-O.  Expanded up to 10× on
  miss; falls back to globally nearest as a last resort.  Atom
  conservation guarantees a match exists; the search radius just
  controls how long a fragment travels.

The dismantle-and-donate algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For every crosslinker residue with fewer than ``full_bond_count``
bonded bridges:

1. **Pick a within-ring C-N matching.**  The 1,3,5-triazine ring has
   exactly two valid pairings — each ring C with one of its two
   adjacent ring N atoms.  Either matching works; the driver
   deterministically picks "each C with the N preceding it in
   traversal order".  Three ``-C#N`` fragments fall out.

2. **Sever the three ring bonds not in the matching.**
   ``delete_bonds`` removes them and cascade-deletes any
   ``[ angles ]`` / ``[ dihedrals ]`` / ``[ pairs ]`` entries that
   referenced them.

3. **Re-tag each (C, N) pair as a new cap-residue.**  Each pair gets
   a fresh ``resNum`` and the atom names get rewritten to ``C1`` /
   ``N1`` so they match the cap template.

4. **Classify each cap fragment as in-place or free.**

   * If the cap's ring-C atom was bonded to a bridge-O during cure,
     it's an **in-place cap**: the bridge-O-to-C bond stays, only
     its parameters need re-templating against the new
     ``BPA-O-C#N`` template.
   * Otherwise it's a **free cap**: it'll be physically relocated
     next to an unreacted bridge-O elsewhere in the system.

5. **Match free caps to unreacted bridge-Os.**  Greedy nearest within
   ``cap_search_radius``; radius expansion up to 10× on miss; global
   fallback if needed.  By atom conservation (one dangling ring C per
   missing bond, one unreacted bridge-O per missing bond, summed
   across the whole system), exactly one match exists per free cap.

6. **For every cap (in-place + free)**: form the bridge-O-to-C bond
   (or re-form it, for in-place caps where it already existed) via
   ``add_bonds_with_template`` against the named ``cap_template``.
   This is the same machinery the cure stage uses to splice
   template-derived parameters into the system around a new bond.

7. **Refresh the C-N bond parameters.**  ``map_from_templates``
   updates the cap atoms' GAFF types (aromatic ``ca``/``nb`` to sp
   ``c1``/``n1``) but only re-resolves bonds it's actively mapping.
   The C-N bond is not one of those, so the driver explicitly resets
   its override parameters to pick up the new atom types.

8. **Batched H deletion and charge rebalancing.**  All sacrificial
   H atoms (dangling ring C-H atoms and bridge-O-H atoms whose
   bridge-O received a free cap) are collected in a single set and
   deleted at the end of surgery, *after* all template splicing is
   done, to keep atom-index references valid throughout.  The lost
   H charges are then redistributed across the heavy-atom neighbours
   of the deleted positions via ``adjust_charges``, restoring net
   charge to 0 for Ewald.

Atom conservation
^^^^^^^^^^^^^^^^^

The matching exists because of an exact equality: across the whole
system at the end of cure,

* the **count of dangling crosslinker-C atoms** (one per
  ``(full_bond_count - k)`` of every incomplete crosslinker, summed
  over all incomplete crosslinkers) equals
* the **count of unreacted bridge sites** (every bridge atom whose
  reactive site didn't bond to a crosslinker, summed over all
  bridges).

This is just bond-conservation: every formed bond consumes exactly
one of each.

Practically: if you run with 240 triazines × 3 = 720 reactive sites
on the crosslinker side and 360 bridges × 2 = 720 reactive sites on
the bridge side, and the cure forms ``B`` bonds, you end up with
``720 - B`` of each kind unreacted.  The repair stage matches them
1:1, and no atoms are wasted.

What it does *not* do
^^^^^^^^^^^^^^^^^^^^^

* Does not modify monomers that aren't part of an incomplete
  crosslinker or an unreacted bridge.  Anything not touched by the
  surgery keeps its existing types, charges, bonds.

* Does not redistribute bonds across the network.  An incomplete
  crosslinker is dismantled in place; the bridges it was bonded to
  retain those connections (with new ``-C#N`` parameters), and the
  network's chain-extension graph is unchanged.

* Does not run during cure.  All decisions are made once, after the
  cure loop has converged and any cap reactions have fired.

* Does not (yet) ship a generic "incomplete trifunctional crosslinker
  → caps" abstraction.  The driver is hardcoded for triazine ring
  topology (6 atoms, alternating C/N, 3-way crosslinker on the C
  positions).  Other ring shapes or crosslinker sizes would need
  their own driver.

.. _postcure_repair_reporting:

What repair reports
-------------------

Repair is also where the conversion a run can honestly claim gets
decided, so the stage reports it.

The cure iterates on **bond conversion**: bonds formed over bonds
possible.  That is the number in every ``Iteration N current
conversion`` line and the number ``desired_conversion`` sets.  But
repair dismantles every crosslinker that did not fill all of its
sites, so the structure leaving this stage contains only *complete*
crosslinkers plus unreacted bridge sites.  The fraction of
crosslinkers that survive intact is a different, lower number, and it
is the one an experiment measures — for a cyanate ester it is the
FTIR ``-OCN`` conversion, because each complete triazine consumes
exactly three ``-OCN`` groups.

The gap is not small and it is not noise.  Under random placement a
trifunctional crosslinker survives only if all three of its sites are
filled, so the crosslinker conversion is roughly the **cube** of the
bond conversion: a run at a bond conversion of 0.90 lands near 0.73.
A run reported as "90 % cured" on the strength of the cure log is, as
a cured structure, closer to three-quarters converted.

**The cube law holds only once the cure has taken enough iterations.**
Below that it overstates the result, and there are two regimes with a
hard boundary between them.

* **Fewer than** ``f`` **iterations: the crosslinker conversion is
  exactly zero.**  The bond downselection admits at most one bond per
  residue per iteration, so a crosslinker with ``f`` sites cannot have
  filled all of them in fewer than ``f`` passes.  This is a counting
  constraint rather than a tendency: it holds whatever the bond
  conversion reached.  A cure that hits its target in two iterations
  builds a monomer melt with no junctions at all, however healthy its
  reported bond conversion looks.
* **At or above** ``f`` **iterations: the iteration count does not
  determine the shortfall.**  Four trifunctional runs that each took
  exactly three iterations came out at 0.06, 0.63, 0.80 and 0.94 of
  the cube-law figure — a fifteenfold spread at an identical iteration
  count, ordered by their bond conversions of 0.40, 0.50, 0.55 and
  0.58.  Measured ratios to the cube law across a series of runs: 0.00
  at two iterations, 0.06–0.94 at three, 0.86 at four, 1.02 ± 0.03
  over 14 runs at five to eight, 1.01–1.04 across the four-run
  nine-iteration cohort — where the spread at three iterations is the
  bond-conversion spread just described, not an iteration effect.
  What sets the *size* of the shortfall looks like the distribution of
  the bonds across those iterations, and specifically how many the
  last iteration formed — it is the one that has to supply a
  crosslinker's final bond, and the per-iteration bond count decays
  steeply, so the average is not what matters.  The run at 0.06 above
  spent its last iteration on 9 bonds against an average of 96, so
  almost no crosslinker could complete; a run at the same bond
  conversion built five months earlier under an older htpolynet spent
  its own last iteration on 10 and landed at the same ratio to three
  digits.  That is the best account available and it rests on few
  runs.  The one-bond-per-residue rule alone accounts for about a
  fifth of the shortfall.

htpolynet warns about both of these at the end of the cure, because
nothing else in the output would tell you.  Below ``f`` iterations it
says so by name, since the cause is then known exactly.  At or above
``f`` it counts how many crosslinkers actually ended up reacted at all
of their sites — an exact number by that point, not an estimate — and
warns when that is a small enough fraction that a percolating network
should not be assumed.  The count is logged either way.  The second
case is the one that matters in practice: a cure can run exactly ``f``
iterations, reach its target bond conversion, and still finish with
almost no complete crosslinkers if its last iteration formed few
bonds.

Proximity works the other way, and it is what closes the gap.  The
search is distance-ranked and a partly-bonded crosslinker already sits
in a bridge-rich neighbourhood, so it keeps getting re-found: a weak
version of ``completion_bias`` for free.  That is a bias of a couple
of percent, not a bound, and runs land on both sides of the cube.

The recovery is complete well before the many-iteration limit, but it
arrives as a two-sided estimate rather than as a floor.  Across 14
independent trifunctional runs at 1.7–2.7 ``f`` iterations and bond
conversions of 0.74–0.90, the crosslinker conversion sat **+2.0 %**
from the cube on average, with a standard deviation of 2.9 % and a
range of −3.7 % to +6.3 %; 3 of the 14 landed *below* the cube,
against a replicate scatter of 0.8–1.8 %, so those excursions are real
and not measurement noise.

That band is bounded from below, and the bound is sharp.  A separate
array of 12 runs at bond conversions of 0.40–0.73 does not behave like
a two-sided estimate at all: 11 of the 12 landed *below* the cube, by
−19.9 % on average, worsening monotonically as the bond conversion
falls and reaching −93.6 % at a bond conversion of 0.40 — an
overstatement by a factor of 15.  So ± 3 % is not the cube law's
general accuracy.  Carried down to a bond conversion of 0.5 it
understates the error by an order of magnitude, and in a predictable
direction: the cube is always the optimistic side down there.

The accuracy is stated against bond conversion because that is what it
tracks — over 26 runs the correlation with bond conversion is +0.81
against +0.57 for the iteration count, and the fifteenfold spread at a
fixed iteration count above rules the iteration count out as the sole
determinant.  It does not follow that the bond conversion is the
governing quantity.  Within a fixed iteration count the *average*
number of bonds formed per iteration is exactly proportional to the
bond conversion, so no set of runs sharing an iteration count can tell
those two apart, and every run here shares one with several others.
Treat the bond conversion as the variable to *read* the accuracy
against; the mechanism, as far as it is understood, is the one
described above, and the bond conversion stands in for how the bonds
were distributed rather than replacing it.

Further out than that it is not settled.  The four runs of the
nine-iteration cohort all landed above the cube, by +0.6 % to +3.5 %,
which reads like a bound — but four points cannot separate a bound
from the upper tail of a two-sided distribution, and the scatter
measured at 1.7–2.7 ``f`` cannot simply be carried to 3 ``f``, since
whether the deviation depends on the regime is the very thing at issue
here.  So the correction below applies where it was measured: at
1.7–2.7 ``f`` the cube is demonstrably not a floor.  At 3 ``f`` it is
untested, not disproved.

So the cube is a useful estimate over a specific band — good to about
3 % at bond conversions of 0.74–0.90 and 1.7–2.7 ``f`` iterations, and
across that band not a floor.  Outside it, in either direction, it
overstates the result: by an amount the iteration count will not tell
you below ``f``, and by up to a factor of 15 at low bond conversion.
The only trustworthy number anywhere is the measured one that repair
reports.

``completion_bias`` does not lift the iteration floor — it changes
which residues react, not the one-per-residue-per-iteration rule — so
it is not the remedy for a low crosslinker conversion caused by too
few iterations.  What it does buy has to be judged against what
proximity already achieves in the many-iteration limit, not against
the random-placement figure.

So at the end of the stage htpolynet logs both, together --
illustrated here with the cube-law figures for a 240-triazine box
taken to a bond conversion of 0.90::

    Crosslinker conversion after repair: 0.733 (176/240 TAZ complete);
    the cure reported a bond conversion of 0.900

and writes the same figures to ``repair-summary.yaml`` in the repair
directory, for reading across a series of runs:

.. code-block:: yaml

   bond_conversion: 0.9
   repairs:
   - residue: TAZ
     n_crosslinkers: 240
     n_complete: 176
     n_dismantled: 64
     crosslinker_conversion: 0.7333333333333333

The stage also reports how many **cap fragments it transferred**, and
how tightly they landed::

    Repair transferred 158 cap fragments; tightest placement 0.128 nm

with the full placement record in the same YAML file:

.. code-block:: yaml

     n_transferred: 158
     n_direction_searched: 31
     n_preferred_out_of_angle: 0
     n_below_target: 17
     min_clearance_nm: 0.128
     median_clearance_nm: 0.196
     blind_min_clearance_nm: 0.008
     blind_median_clearance_nm: 0.161
     n_blind_would_overlap: 9

These do not all mean the same kind of thing, and reading them as if
they did is the mistake to avoid:

* ``n_below_target`` is the one **achieved** number that describes the
  tail, which is where a placement failure actually lives.
* ``min_clearance_nm`` and ``median_clearance_nm`` are **placement
  outcomes, not measurements of the site**.  The search stops at the
  first direction that reaches the target, so both are pulled toward
  ``cap_min_clearance`` by construction and partly report the
  threshold you set rather than how crowded the box is.  The median is
  pulled; the minimum is pinned.  Whenever the search succeeds for
  every cap, the worst-placed cap is by construction one that only
  just cleared the target, so ``min_clearance_nm`` lands just above it
  — measured at 0.1503 ± 0.0006 nm against a 0.150 nm target across 14
  independent boxes, and the single box that came in below (0.1488)
  was the single box with a non-zero ``n_below_target``.  So
  ``min_clearance_nm`` carries nothing that ``n_below_target`` has not
  already said, and in particular it is not something to calibrate
  ``cap_min_clearance`` against: it is the threshold seen from above.
* ``blind_min_clearance_nm`` and ``blind_median_clearance_nm`` are the
  crowding statistics, and the blind minimum is the real tail
  statistic — it ranged 0.006–0.048 nm across those 14 boxes while the
  achieved minimum sat on the target.  Both are measured along the
  single fixed O-H direction with no search and no early exit, so they
  describe the site rather than the algorithm, and they are the ones
  to correlate across a series of runs and to calibrate the target
  against.  ``n_blind_would_overlap`` counts the caps that direction
  alone would have put inside 0.10 nm of a neighbour — which is what
  placement did before the direction search existed, and which came to
  37 % of caps on those boxes.
* ``n_preferred_out_of_angle`` separates the two reasons a preferred
  direction gets abandoned.  A large value means the C-O-C angle
  window is rejecting real O-H vectors, which would look exactly like
  a crowded box in every other field.  On the 14 boxes above it was
  6 % of caps, so the 90–150° window is not what is sending caps to
  the search; a value several times that is the one to act on.

That count is worth watching, because it is the quantity that predicts
whether this stage survives.  It is an identity — every bond the cure
did not form leaves a bridge site unreacted, and every unreacted
bridge site receives one transferred fragment::

    transferred fragments = total reactive sites - bonds formed

so it *rises as conversion falls*, and a low-conversion build can be
placing several hundred.  Each one is a two-atom group being dropped
into an already-dense box.  If the following minimization dies with a
Lennard-Jones term of order 1e15, this is where to look, and the
warning about caps that could not reach ``cap_min_clearance`` names
the specific oxygens involved.

Note that this count depends only on the bond conversion, not on how
the bonds are distributed, so ``completion_bias`` does not reduce it.
At a *matched crosslinker conversion* the bias in fact increases it,
because it reaches that conversion by forming fewer bonds.

To *raise* the crosslinker conversion rather than merely report it,
see the ``CURE.controls.completion_bias`` directive in
:ref:`the run-configuration reference <cure.controls>`: it makes the
cure finish partly-reacted crosslinkers before starting untouched
ones, which is both the more physical rule for a cyclotrimerizing
chemistry and the thing that empties this stage's work queue.

Writing a new repair driver
---------------------------

A repair driver is a Python function with the signature

.. code-block:: python

   def my_driver(TC, moldict, spec, reactions):
       # ... arbitrary topology surgery on TC ...
       return {'residue': ..., 'n_crosslinkers': ..., 'n_complete': ...,
               'n_dismantled': ..., 'crosslinker_conversion': ...}

* ``TC`` — the :class:`htpolynet.core.topocoord.TopoCoord` for the
  cured system.  Modified in place.
* ``moldict`` — the :class:`htpolynet.core.molecule.MoleculeDict` of
  all parameterized templates, including any repair-stage linked
  products.
* ``spec`` — the dict from this entry of the YAML's
  ``postcure_repair`` block.
* ``reactions`` — the full ``ReactionList`` (rarely needed; useful
  for cross-checks against the configured reaction set).

The driver should return a statistics dict.  ``n_dismantled`` is the
count of operations performed, used in the runtime log message and
summed across drivers; the remaining keys feed the conversion report
described in :ref:`what repair reports <postcure_repair_reporting>`
and may be omitted by a driver for which crosslinker completion is not
a meaningful notion.  It is free to call into
:mod:`htpolynet.repair.topology_surgery` for the heavy lifting and
into :class:`htpolynet.core.topocoord.TopoCoord` and
:class:`htpolynet.core.topology.Topology` directly for anything
finer-grained.

To register the driver, edit
``src/htpolynet/repair/__init__.py:run_repair`` to dispatch the new
``type:`` value at the top of the function.

If the driver needs a parameterized linked-product template, the
recommended pattern is to declare it via a ``stage: repair`` reaction
in the YAML (so the existing setup-time parameterization machinery
builds it), and to look it up by ``product`` name at surgery time.
