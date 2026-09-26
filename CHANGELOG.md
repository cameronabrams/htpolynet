# Changelog

All notable changes to htpolynet will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.11.4] - 2026-09-26

### Added

- **A file of your own that replaces a packaged one now says so.**  `pfs.checkout`
  searches the user library, then the user cache, then what htpolynet ships, and took
  the first hit silently.  Overriding a packaged file is supported and often the point;
  the trap is that a copy taken from an older version keeps that version's settings for
  as long as it stays on the search path, and nothing said which file a run actually
  used.

  Found the hard way.  A study's `npt.mdp`, copied from 2.6.2 and edited in one place,
  stayed on its path across four upgrades.  From 2.7.0 the packaged file sets
  `lincs_order = 8`, because at `dt = 0.002` GROMACS' default of 4 is marginal --- it
  was fatal on one chemistry in 7 of 7 builds, with exit codes 1 and 139.  The stale
  copy reverted every run to order 4, and the instability that produced was written up
  for weeks as a property of the systems being built.

  **This bears on the 2.11.3 notes below.**  "Eleven of twenty builds died in the
  anneal" is the rate from builds carrying that stale mdp, so it overstates what the
  unrelaxed final batch alone causes.  The defect and the remedy are unaffected: the
  bond energies entering the anneal are single-point energies of the structure and owe
  nothing to integration settings, and the paired seed test showed the settle helping
  even at LINCS order 4.  Only the headline failure rate is in doubt.


- **A ring cure now refuses a ring that would close around an existing bond.**  CURE
  tests every candidate bond for piercing an existing ring; the ring cure has the
  converse problem and tested nothing at all --- `ringcontroller.py` had no reference to
  `bondtest`, `pierces_ring` or `BTRC`.  A triazine could therefore form *around* a
  monomer that was already there and thread it permanently.

  That is what the stretched monomers were.  Over ten finished builds, six carried a
  threaded triazine, and the threading bonds are the same bonds found stretched past
  2 A: 15 stretched, 13 threading, none threading without being stretched.  10 of 11
  were pre-existing monomer backbone bonds.  It explains why those bonds sit on monomers
  that took no part in the ring, why a 160 ps anneal left them unchanged to two
  decimals, and why they hold roughly 700 kT without being thermal --- escaping a
  topological trap needs a bond to break, not a barrier to be crossed.

  The test runs after the closure ladder and before the bonds are formed, which is
  where the question is unambiguous: the loop is down to nearly its final size, so what
  is inside it is what the triazine will enclose.  It treats the loop as a triangle fan
  from its centroid rather than as a plane, and images every coordinate against an
  anchor atom, since a centroid of raw wrapped positions gives false positives across a
  periodic boundary.

  `ring_cure.prefilter_threaded_candidates` applies the same test at candidate time as
  well, to avoid spending a ladder on a ring that is already threaded.  It is **off by
  default**, because a candidate loop is three groups a search radius apart: a formed
  triazine's circumscribed radius is about 1.37 A and near enough rigid across 233 of
  them, while the candidate loop encloses roughly nine times that area at a 1.0 nm
  search radius and over twenty at 1.6 nm, which is where late iterations run.
  Threading would carry through the shrink if nothing moved, but the ladder runs six
  stages of NVT at 600 K, so a bond inside the candidate loop need not still be inside
  at closure.  Set `ring_cure.reject_threaded_rings: false` to reproduce the old
  behavior entirely.

  **Measured across three arms, same image, no mdp override.**  Filtered builds came out
  clean 6 of 6; unfiltered builds carried a pierced triazine in 2 of 4 --- 4, 2, 0 and 0
  piercings, a mean of 1.5 each.  There is no conversion cost: 0.971 to 0.983 in 21 to
  23 iterations against 0.971 to 0.979 in 20 to 22 without it.

  The closure-time test declined 2, 10 and 1 candidates across three builds, a mean of
  4.3 against those 1.5 piercings.  Same order of magnitude with overlapping
  distributions, and that is as much as the data supports: at three builds against four
  the smallest two-sided Mann-Whitney p obtainable is 0.057, so the observed 0.47 is no
  more evidence of agreement than of excess.  The two are not the same quantity either
  --- a declined candidate is one that *would* have threaded, and the alternative
  `select_disjoint` takes instead may itself thread and be declined in turn, so declines
  bound the piercings prevented rather than counting them.

  The candidate-time prefilter is a different matter: it discarded 1216 to 4385
  candidates per build to catch what the closure test catches in single digits, and both
  arms finished equally clean --- which is why it stays off by default.

### Fixed

- **A ring cure ran at constant volume from first ring to last, so it could not
  densify.**  Every ring-cure stage was NVT — the relaxation after each iteration, the
  closure ladder, and the settle — with no constant-pressure step anywhere.  Measured
  over 22 iterations of a 233-triazine build, the box stayed at 5.241 nm to four
  decimals and the density at 1155.9 kg/m3 from the first ring to the last, with the
  entire volume change deferred to a single postcure NPT.  CURE's density climbs
  through its cure, roughly 1085 to 1170, because its relaxation stages run NPT above
  Tg and are deliberately "the only above-Tg constant-pressure time in a cure".

  Cure shrinkage cannot be reproduced at fixed volume by construction — Snow measures
  0.0324 ml/g for BADCy — and a network that forms in a box which cannot respond
  accumulates internal stress with nowhere to put it, which is a plausible route to the
  stretched bonds recorded in `ROADMAP.md`.

  `ring_cure.relax` now mirrors CURE's `relax.equilibration` exactly --- minimization,
  NVT at 600 K, then NPT at 600 K and 1 bar --- rather than CURE's much longer
  per-iteration equilibration.  The whole sequence is above Tg, not just the
  constant-pressure stage: relaxing at 300 K in between two 600 K closure ladders cooled
  the network and would have put a cold, less mobile configuration under the barostat,
  which is a good way to measure no densification and wrongly conclude the stage does
  not work.  The closure ladder itself stays constant-volume on purpose: NPT while
  restraints are actively pulling groups together is not a stable combination.

- **A ring cure had no cold equilibration, so once it could move, its box only grew.**
  Adding the constant-pressure relax stage above let the box respond, but `relax` ends
  hot.  Measured with relax alone over a full build, the box expands thermally for six
  iterations to +27.9% in volume, then comes back down as the network forms: peak to
  minimum is -7.49% in volume and +71.9 kg/m3.  The two legs nearly cancel, so the net
  is near zero and the path length is 0.4104 nm --- a frozen box captures none of it.

  That descending leg is not a measurement of cure shrinkage: 2 ps at 600 K does not
  equilibrate, so it mixes in residual thermal relaxation, and it spans conversion 0.55
  to 0.971.  It costs iterations too --- `search_radius` is an absolute length and does
  not scale with the box, so an expanded box encloses fewer groups per sphere, and the
  hot-only arm needed 25 iterations against 22 for the same conversion.  A cure is not only shrinkage --- the hot stages enlarge
  the box and a cold stage pulls it back, and on the pairwise route the box travels
  three to nine times further in length than it ends up moving.

  `ring_cure` now has an `equilibrate` stage identical in shape and in defaults to
  `CURE.equilibrate` --- NPT, 300 K, 1 bar, 50000 steps of 2 fs --- run once per
  iteration after `relax`, through the same code path the pairwise cure uses, so it
  picks up that route's box logging and density series rather than reimplementing them.
  The per-iteration shape of the two cures is now the same.

  It is also the expensive stage, roughly 99 s against a 280 s iteration.  Whether
  100 ps is the right number for *either* route has never been measured; `ROADMAP.md`
  sketches the study.

- **Three `ring_cure` settings could not be set at all**, because they were missing from
  the schema, and ycleptic rejects an unknown key outright rather than ignoring it.
  `ring_cure.relax` and `ring_cure.closure.equilibration` were never exposed, and
  `closure.max_accept` — added in 2.11.1 — was written into `CURE`'s `drag` block by
  mistake, where a config could set it and nothing would read it.  All three now sit
  where they belong, and a test asserts that every key the controller honors appears in
  the schema with a matching default.

### Changed

- **The installation page now says not to install `htpolynet` from conda-forge.**  The
  package published there is several releases behind and cannot currently be updated:
  `htpolynet` depends on `ycleptic`, which is not yet available on conda-forge, so the
  feedstock's automatic version bumps produce a package that cannot import and none has
  been published.  What is on conda-forge does install and run, so nothing warns a user
  that they have an old `htpolynet`.

  The conda-only recipe now installs Python and the MD binaries from conda-forge and
  `htpolynet` itself from PyPI.  The warning is scoped to the `htpolynet` package alone
  --- conda-forge remains the recommended source for AmberTools and Gromacs --- and it
  gives the two commands to check the current state rather than asking the reader to
  trust a page that may itself be stale.  The container remains the recommendation, and
  the page now says outright that it requires no installation at all.

### Added

- **A build now says so when it hands on a bond that is not a bond any more.**  At the
  end of a ring cure and again in the final data, any bond longer than 0.2 nm is
  reported, the worst named by atom and residue.  Lengths use the minimum-image
  convention, so a bond across a periodic boundary is not counted.

  This exists because such a bond is otherwise silent: the topology is valid, every
  molecule is neutral, and the run completes.  htpolynet-study measured a C-C at 3.20 A
  in two of five *final* structures, unchanged to two decimals by a 160 ps anneal.  Ring
  bonds always relaxed; these did not, because they are held by network tension rather
  than by thermal motion — a monomer stretched between parts of the network that a ring
  closure pulled apart, in a ring it was not itself part of.

  This is detection and not a fix; `ROADMAP.md` has the mechanism and why no remedy has
  been chosen.  Anything computed from a structure carrying such a bond inherits it.

## [2.11.3] - 2026-09-23

### Fixed

- **A ring cure handed its last batch of rings to postcure unrelaxed.**  Every batch
  but the last is relaxed twice: once by the per-iteration relax as soon as its bonds
  exist, and again by the *next* iteration's closure ladder, which is eight
  minimizations and 16 ps at 600 K against the relax's one minimization and 4 ps at
  300 K.  The final batch only ever gets the first, and then goes straight into a
  500 K anneal.

  htpolynet-study measured it on three 2.11.1 builds, all of which had ended on a
  ring-forming iteration and had therefore been relaxed.  Bond energy entering the
  anneal tracked how many rings that last iteration formed -- 1, 1 and 4 rings for
  18849, 26113 and 35959 kJ/mol, against 13619 +/- 226 for a route that enters already
  relaxed.  An inventory of one such structure found five bonds beyond 2.0 A: one
  triazine's three ring bonds at 2.24 A, and two monomers whose own backbones the
  ladder had stretched to 2.34 and 3.20 A while pulling their arms toward different
  rings.  Eleven of twenty builds died in the anneal.

  The cure now runs a minimization and a short 300 K NVT once on exit, mirroring what
  the postcure repair stage already did for the route that never saw this.  It is
  unconditional -- the cure has four ways to end and which of them leaves strain behind
  is not worth reasoning about per exit -- and configurable as `ring_cure.settle`.

  Confirmed paired, same eight velocity seeds against the same two structures on the
  same hardware: two anneals exploded without the settle and none with it.  The two
  came from different source structures, one of which had never failed in production,
  so the effect is not particular to a pathological build, and nothing exploded during
  the settle itself -- the minimization absorbs the strain rather than the failure
  moving earlier.  Eight seeds is a small sample and that comparison alone would be
  weak; it is corroboration of the bond energies above, not the case on its own.

## [2.11.2] - 2026-09-23

### Fixed

- **A ring cure near its target kept re-proposing the triples it had just declined.**
  Declining a ring that the closure ladder could not pull shut leaves its three groups
  unreacted, which is the point -- but the next search then offered the same groups,
  ran the same eight-stage ladder, and declined them again.  htpolynet-study measured
  the cost on a 2.11.1 build: at conversion 0.97 with 24 groups left, nine iterations
  running produced 77 declines and not one ring, at about 67 s each.

  A declined combination is now passed over while any untried ring is within the search
  radius.  It is a preference and not a prohibition: that cure escaped on iteration 41
  and finished at 0.971, because the relaxation between iterations moves the groups, so
  a triple that failed once can still close later.  When the memory is the only thing
  left keeping an iteration idle, the search lifts it and offers those rings again
  rather than stand still.

## [2.11.1] - 2026-09-22

### Fixed

- **A cyclotrimerization's charges were computed on a ring that was never closed.**  The
  trimer template reached antechamber with its three ring bonds recorded but still about
  2.2 A long: the closure ladder cannot beat the reactants' own C#N, whose sp parameters
  hold each group linear, and no force constant wins that.  Charges come from the
  geometry, so the semiempirical step charged three separate cyanates rather than a
  triazine, leaving **every ring carbon about 0.21 e low** with the deficit smeared over
  the rest of the molecule.

  Nothing caught it.  The atom types were right, because antechamber perceives those from
  the connectivity, which was right; the molecules were neutral, because the error sums
  to zero; and the system's own ring bonds relaxed to 1.347 A after the splice, because
  the bonded parameters were fine.  It surfaced only when htpolynet-study ran a campaign
  with `bcc` charges and half their builds died in a 500 K anneal -- the release smoke
  test had used `gas` charges, which are computed from the topology and could not have
  seen it.

  A product whose reaction closes a ring is now relaxed against its own connectivity
  before it is parameterized, and refused if the ring is still open afterward.
  Measured on a methyl cyanate trimer, ring-carbon charge: **+0.699 before, +0.909
  after, against +0.909** for the same molecule built independently from SMILES.
  htpolynet-study measured the same 0.21 e deficit on a bisphenol A dicyanate trimer.

  Two things fell out of the investigation that are worth knowing:
  - **Correcting the bond orders does not fix it.**  An addition bond leaves the cyanate
    carbon pentavalent -- it keeps its C#N, gains the ring bond, and still carries its
    ester oxygen -- and that is now repaired too, since RDKit cannot read a molecule with
    an impossible valence.  But repairing the orders alone produced charges identical to
    the broken ones, because the semiempirical step reads coordinates, not bond orders.
  - **antechamber's output geometry is its input geometry.**  It is not evidence about
    what the charge calculation did or did not relax, and it should not be read as such.

- **A ring the closure ladder could not pull shut is now declined rather than bonded
  long.**  Late in a cure the matrix is rigid enough that a ring sometimes stops well
  short -- one batch at conversion 0.95 ended 0.327 nm apart against a target of 0.150 --
  and forming the bond anyway leaves it long for good: two ring bonds survived relaxation
  and a 500 K anneal still 2.6 A apart, strained but stable, a defect the network then
  carries permanently.  Declining costs almost nothing, since the groups stay unreacted
  and are offered again next iteration once the neighbourhood has moved.  The limit is
  `ring_cure: closure: max_accept`, default 0.30 nm; 0 restores the old behavior.  Only
  visible now because the charge fix let a build reach conversion 0.95 for the first
  time.

- `Topology.rebalance_mol2_bond_orders` reduces a multiple bond when an addition would
  otherwise leave an atom over-valent, preferring a bond that is over-valent at both ends
  because reducing it repairs two atoms at once.

## [2.11.0] - 2026-09-22

### Added

- **A three-body cure reaction: cyclotrimerization** (see `ROADMAP.md`).  A
  `reactions` entry whose bonds close a ring among three reactants, plus a
  `ring_cure` configuration block, now build a network by the real cyanate-ester
  chemistry -- three -O-C#N groups closing one triazine -- instead of the
  pre-formed-triazine surrogate that needs postcure repair.  CURE is untouched
  and remains the default; a configuration without a ring-closing reaction
  behaves exactly as before.

  Verified end to end twice.  On 90 methyl cyanates, which can only trimerize:
  19 triazines, 63% conversion, every ring six-membered with three aromatic
  carbons and three aromatic nitrogens, bonds averaging 1.351 A (a triazine's
  is 1.34) and planar to 0.065 A, the 33 unreacted groups keeping their C#N at
  1.153 A, and all 52 molecules neutral.  Then on 100 ethylene dicyanates,
  which gel: 50 triazines, 75% conversion, ring bonds averaging 1.347 A, 55
  monomers reacted through both arms and 40 through one, no chemically
  impossible bond, and the largest of the 8 remaining molecules holding 1068 of
  1200 atoms -- a percolating network.

  Then at scale, on 150 bisphenol A dicyanates -- the real monomer, 35 atoms
  rather than 12, with a rigid aromatic backbone and its two arms about 0.9 nm
  apart: 61 triazines at 61% conversion, every unit six-membered with three
  aromatic nitrogens and three ring carbons still carrying their ester oxygen,
  ring bonds averaging 1.347 A, no bond anywhere longer than 2.5 A, 54 monomers
  reacted through both arms and 75 through one, all 150 molecules neutral, and
  the largest of the 30 remaining components holding 1295 of 5250 atoms.

  A monomer with more than one reactive site needs two things beyond the
  single-site case, both of which the difunctional build exposed:
  - **Only the reacting site's share of a residue is spliced.**  A dicyanate can
    join two rings, and the trimer template shows one arm reacted, so splicing
    the whole residue reset the other ring's arm -- 36 bonds between an aromatic
    nitrogen and an unreacted cyanate carbon in the first attempt.  Each atom is
    assigned to whichever site is nearer through the residue's own bonds.
  - **The template's atom names are translated for a residue that rings through
    its other site.**  That exchange moves every atom of the arm, hydrogens
    included, so it is derived as an automorphism of the monomer's bond graph
    rather than read off `symmetry_equivalent_atoms`, which names only the
    reactive atoms.  The sets are still consulted for which sites exist.

  A ring-closing reaction is also exempt from symmetry expansion, which permutes
  each reactive atom independently: for three dicyanates that produced 64
  reactions, most of them mixing the two arms of one site, and a template for
  each.
  - A reaction's residue-offset arithmetic now works for any number of
    reactants.  It was written for exactly two, so a bond between the first and
    third reactant of a three-reactant reaction resolved to the wrong residue.
  - `reaction.is_ring_closing` recognizes a reaction whose bonds close a cycle
    among its reactants, which a cyclotrimerization does (1-2, 2-3, 3-1) even
    though its new bonds alone form no cycle -- the ring closes through bonds
    each reactant already had.
  - The pairwise bond search now skips such a reaction with a warning naming
    it, rather than treating its bonds as three unrelated pairwise candidates
    and forming them in separate iterations, which would not build a ring.
  - A reaction bond can declare `sacrificial_h: false`, for an addition where
    nothing is lost.  Every reaction htpolynet had before is a condensation
    that deletes one hydrogen from each of the bond's atoms, which stays the
    default; cyclotrimerization forms three bonds and deletes no atom.
  - `reaction.spanning_and_closing_bonds` separates the bonds that assemble a
    product from the bond that closes its ring.  A template build positions
    each incoming piece as it bonds it, which a closing bond cannot do: the
    residues it joins are already placed, so it is skipped there and its
    geometry handled separately.
  - `Runtime.do_ring_cure` runs it as a workflow stage, with restart state
    beside the cure's.  `do_cure` now returns early when every cure reaction
    closes a ring: the pairwise search skips those, so it would otherwise
    iterate to its limit forming nothing.  (The `CURE` block is present even
    when a configuration omits it, since the schema supplies its defaults.)
  - `cure.ringcontroller` runs the loop: search, pull every chosen ring shut,
    form its bonds, splice the trimer template over its three residues, settle
    the charge, and count what was consumed.  Conversion is counted in reactive **groups**, three per
    ring, which is what this chemistry reports -- the pre-formed-triazine route
    has to infer it from a bond count instead.  The search widens until an
    iteration yields `min_rings_per_iteration` rings, as the pairwise cure's
    does for bonds.  Its state round-trips to YAML for restart, carrying the
    widened radius with it.
  - `cure.triplesearch` finds the triangles of reactive groups one ring-closing
    reaction can join.  The pairwise search cannot: it asks which A atom is near
    which B atom, while a cyclotrimerization needs three groups mutually close,
    and two well-chosen groups can still leave the third far off.  Candidates
    are scored as whole triangles by their total new-bond length, both ring
    orientations are tried, each group may join only one triple, and a group is
    identified by residue *and* atom names -- a bisphenol dicyanate carries two,
    which is what makes it a crosslinker.  Two groups of one residue, or of one
    molecule, are excluded by default.  On a 720-molecule box (1440 groups) at
    a 0.5 nm radius it finds 358 candidates in half a second and packs 148
    triples from them.
  - A reaction bond that consumes no hydrogen is positioned by
    `Molecule.place_for_addition` rather than `transrot`, which aligns two
    reactants on the hydrogens they are about to lose and so cannot place an
    addition at all: in a cyclotrimerization neither the cyanate carbon nor its
    nitrogen carries one.  The new placement puts the incoming piece a bond
    length away, in the direction and at the turn that leave the most room,
    using one shared copy of that geometry in `geometry.placement`.
  - `Molecule.close_ring_geometry` pulls a ring shut before its bonds exist,
    by restraining every pair of the ring and stepping the restraints down to a
    bond length under minimization.  The merged reactants' topology is already
    valid for MD, and a type-6 restraint needs no atom types, so this needs
    nothing parameterized.
  - `Topology.set_restraint_parameters` prescribes a restraint's length and
    stiffness outright.  `attenuate_bond_parameters` cannot serve here: it
    reads each bond's reference values from the row it overwrites, which is
    stable for a real bond, whose reference comes from its atom types, but
    compounds on a restraint carrying explicit parameters -- an eight-stage
    ladder ended 400 times too weak to pull anything.
  - `core.productsplice.map_product_from_template` makes the residues of one
    reaction event match a product template exactly -- atom types, charges, and
    every bond, angle, dihedral and 1-4 pair among them, replacing the
    instance's parameters where they exist and adding them where they do not.
    This is what a ring closure needs and what the per-bond splice cannot do:
    among the interactions a triazine closure changes are ones wholly inside a
    residue (a cyanate C#N becoming a ring bond), which contain no new bond and
    so are never revisited by `map_from_templates`.  It takes the residue
    correspondence from its caller rather than inferring it from local context,
    refuses a residue map that does not explain the bonds it was given, and
    refuses by default to run before those bonds exist.

### Fixed

- **A ring cure now widens its search radius when an iteration is short of
  rings, not only when it finds none.**  The pairwise cure has always grown its
  radius within an iteration until it has `min_bonds_per_iteration` bonds; the
  ring loop instead widened only after an iteration that found nothing at all.
  A floppy monomer never noticed, but a rigid one keeps trickling one or two
  rings at the starting radius, so the radius never grew: 150 bisphenol A
  dicyanates stopped at conversion 0.40 after eight iterations, having asked for
  0.60 and never once widened.  With the same rule the pairwise cure uses --
  counted in rings, and clamped against the rings still needed to reach the
  target -- the same system reaches 0.61 in eight iterations.  The floor is
  `ring_cure: min_rings_per_iteration`, default 4.

- **A GPU build died at its first minimization when the configuration set
  `ntomp`.**  GROMACS refuses OpenMP threads together with a GPU unless it is
  also told how many ranks to run: "Please specify the number of thread-MPI
  ranks as well (option -ntmpi)."  Nothing supplied one, so every such build
  failed about a minute in, after parameterization had already been paid for.
  htpolynet now passes `-ntmpi 1` when a configuration asks for a GPU, sets
  `ntomp`, and gives no rank count of its own.  A rank count in the
  configuration (`ntmpi`, `nt` or `npme`) is left as it is, and so is an MPI
  `mdrun`, which takes its ranks from `mpirun` or `srun` and rejects
  `-ntmpi`.  Reported by htpolynet-study from a 9-task GPU build where every
  task failed.

- **`ntmpi` passed as a keyword argument to `grompp_and_mdrun` was dropped**,
  because the passthrough list spelled it `ntpmi`.  Configurations were
  unaffected: `mdrun_options` is copied whole.

## [2.10.3] - 2026-09-17

### Fixed

- **The CPU and CUDA container images shipped different AmberTools.**  The
  CUDA images, `:cuda-v2.10.0` through `:cuda-v2.10.2`, had AmberTools 24.8,
  while the CPU images had 26.0. Nobody chose that: conda-forge can't install
  AmberTools 26 alongside CUDA GROMACS, and the unpinned Dockerfile line fell
  back silently. So the same configuration was parameterized by different
  AmberTools releases depending on which image ran it. Both images now
  install AmberTools in their own conda environment, pinned to one version
  (26.0, set by the `AMBERTOOLS_VERSION` build argument). A future conflict
  now fails the build instead of quietly downgrading.

- **A cached parameterization is no longer reused across AmberTools
  versions.**  The cache record now includes the AmberTools version, and a
  run on a different version re-parameterizes instead of reusing the entry.
  Records written before this, and runs that can't tell their version, still
  match. A cache filled with a CUDA image up to v2.10.2 (AmberTools 24.8)
  has no version in its records, so it would still be reused by a newer
  image. Refill it if the charges must come from one AmberTools release.

## [2.10.2] - 2026-09-17

### Fixed

- **Templates for crowded monomers could start with atoms on top of each
  other, and AM1-BCC/ABCG2 charges then failed.**  When htpolynet built a
  linked-product template, it attached the incoming molecule along the new
  bond but never turned it about that bond. Its clash check also compared
  the moved piece with its own pre-move position. With tetramethylbisphenol
  F (TMB), whose methyls flank each reacting oxygen, every triazine template
  started with a methyl hydrogen 0.65-1.26 A from a ring atom. On one triply
  substituted template, sqm (AmberTools 24.8, as shipped in the `:cuda`
  image) then failed to converge ("Cannot properly run sqm"), which stopped
  parameterization.

  The placement now scans turns about the new bond, in 10 degree steps, and
  keeps the one farthest from the residues that stay put. Every TMB template
  now starts at 1.50 A or more. The six triply substituted TMB templates
  were run with ABCG2 on that same AmberTools 24.8: five, including the
  geometry that failed before, completed, and the sixth converged its SCF
  and was still optimizing when this was written. On these crowded molecules
  the optimization is slow, 40 minutes to 2.6 hours each, run three at a
  time. Charges and atom types of all 123 templates in the shipped examples,
  which use gas charges, are unchanged. A template that still starts with a
  contact under 1 A is now logged as a warning.

## [2.10.1] - 2026-09-16

### Fixed

- **Asking for `pme: gpu` or `update: gpu` in `mdrun_options` stopped the
  build at its first minimization.**  The same options go to every `mdrun`,
  and GROMACS refuses both tasks on a GPU without a dynamical integrator.
  For minimization stages they now become `auto`, and so does `bonded`, which
  hasn't been checked with a minimizer.  Dynamics gets the options unchanged.

- **The parameterization summary undercounted and named the wrong place.**
  It was printed before the symmetry-, chain- and second-shell-derived
  templates were built, so a run that reused 29 templates reported 11.  It
  also always named `~/.htpolynet`, even when `HTPOLYNET_CACHE` pointed
  elsewhere.  It now comes after every template exists, reports "N of M
  reused", and names the library and cache actually searched.

### Changed

- **The "percolation not assessed" message says when it means "cannot
  assess".**  It read as though the outcome were always the correct answer
  for a linear polymer, with the multi-residue crosslinker case as an aside.
  For example 2 that aside is the actual situation: bis-GMA crosslinks
  through junctions a residue-level count cannot see.  Noticed by a sweep of
  the shipped examples.


## [2.10.0] - 2026-09-15

### Added

- **A constituent can carry its own frcmod.**  `frcmod: file.frcmod` on a
  constituent loads that Amber parameter file into `tleap` after GAFF and
  `parmchk2`, so its parameters override theirs.  It is loaded for the
  molecule and for every template containing it, including cure, capping,
  repair and multiply-substituted templates.  The file's contents are part of
  the parameter cache record, so an edited frcmod re-parameterizes instead of
  reusing stale numbers.  A system has one table per bonded type, so setup
  stops with the type's name if an overridden type also appears, with other
  parameters, in a molecule built without the file.

## [2.9.0] - 2026-09-15

### Added

- **The build log reports molecules that carry net charge.**  After cure,
  after postcure repair, and when the final results are written, htpolynet
  checks every covalently bonded molecule.  It logs one line when all are
  neutral to within 0.01 e, and a warning naming the count and the largest
  charge when some are not.  The system total is zero by construction, so
  until now a template or repair step that moved charge from one molecule to
  another left no trace.  The check can't see charge misplaced *inside* a
  percolated network, because the whole network is a single molecule.

### Fixed

- **Later bonds on a triazine reset the charges of its earlier ones.**
  Each cure template shows a single BPA on the ring, and splicing it in
  rewrites the charges of atoms up to two bonds from the new bond, which
  includes the other ring carbons.  So every ring carbon that reacted before
  the last one got its *unreacted* charge back.  Example 6's fully reacted
  rings had carbons anywhere from 0.14 to 0.31 e, where the template value is
  0.28.  A cured bisphenol-A dicyanate network showed a clean 2:1 split,
  +0.70 against +0.90 e.  htpolynet now builds extra templates at setup for
  rings that already carry one or two BPAs; example 6 gets 18 of them.  Each
  cure bond picks the one that matches the ring as it stands: templates
  record which atoms two bonds from the bonding atom have already reacted,
  and matching prefers the one that records the most.  This applies to any
  reactive atom whose `symmetry_equivalent_atoms` partners sit two bonds
  away.  None of the other shipped examples has one, and example 5's
  build-stage molecules are unchanged, charge for charge.

  Replaying example 6's recorded cure bonds offline puts every bonded ring
  carbon at its template value: 0.275, 0.277 or 0.279 e for one, two or
  three BPAs.  All 13 molecules after cure are now neutral, where before
  they ranged from -0.37 to +0.29 e.  Cyanate-ester builds should be
  rebuilt.

- **A template that matched a bond only in reverse crashed the cure** with
  `mapping mismatch -- bug`.  After reversing the bond to line up with the
  template, the mapping check still compared it the original way round.  No
  shipped example had reached this path.

- **Cyanate-cap repair left molecules charged.**  `triazine_to_cyanate_cap`
  kept the system neutral, but not each molecule.  It settled charge twice,
  both times system-wide:
  - once while splicing the cap template, when the sacrificial H's still
    existed;
  - once after deleting them, onto only the atoms those H's had been bonded
    to.

  Caps formed in place and caps transferred to a free BPA-OH therefore ended
  with charges 0.24 e apart on the phenolic O and the cap C.  In a
  zero-conversion bisphenol-A dicyanate melt, every molecule came out at
  either +0.21 or -0.27 e.  In example 6, the 27 molecules left after repair
  ranged from -3.59 to +0.62 e.  Repair now brings each molecule it touched
  back to neutral, and makes the correction only on the atoms it changed.
  Example 6's molecules are all neutral, and its cap residues sit at
  0.05-0.09 e, near the template's 0.069, instead of splitting between 0.04
  and 0.20.  Builds that used `postcure_repair` should be rebuilt.

- **A false percolation warning on linear polymers.**  The end-of-cure check
  warned that too few residues had reacted at every site to assume the system
  percolates, and it fired on chain-growth systems like example 1's
  polystyrene -- "only 902 of 1000 STY (90.2%) ... below the 100.0% ... at its
  gel point".  For a difunctional monomer that threshold, `(1/(f-1))^f`, is
  exactly 100%, and a fully reacted difunctional unit is a chain interior
  rather than a crosslink junction, so neither the warning nor its premise
  applied.  The check now runs only when some residue has at least three
  reactive sites.

- **The percolation check counted sites the cure never uses.**  It summed
  every reactive site a residue was built with, including bonds made while
  assembling molecules and valences no cure reaction can consume.  Example 2
  warned that "0 of 150 HIE" had reacted at all four sites, which was
  impossible: two of those sites sit on an atom only the molecule-building
  reactions touch.  Example 5's warning about TB was the same mistake, since
  every TB site is a chain bond made before the cure.  In examples 3 and 4 the
  diepoxide and the diamine tied at four sites, so the reported fraction
  mixed both residue types ("83 of 300 DGE").  Functionality and completion are
  now counted only on atoms a cure reaction can bond, taken from the
  symmetry-expanded reaction list the bond search itself uses.  Example 3 now
  reports "83 of 100 PAC", and example 6's check is unchanged.  When no
  residue has three or more cure sites, the check says crosslink percolation
  was not assessed, rather than staying silent: systems that crosslink through
  a multi-residue molecule, such as a dimethacrylate, are outside what a
  residue-level count can see.

- **Bondless MOL2 files for symmetry-sibling molecules.**  A molecule built by
  copying a parent's topology (for example `GMAS-4`, from `GMA`) has its bonds
  in the topology but no MOL2 bond table, so writing it as MOL2 produced a file
  with no bonds and a "Cannot write any bonds to MOL2 file" warning -- twelve
  times in example 2 on a fresh parameterization.  Those writes are now filled
  from the topology's bonds, as single bonds.  In every shipped example the
  affected files were diagnostic `*-prebonding.mol2` outputs that nothing reads,
  so no parameterization input was ever affected.

- **The CPU container image failed to build for 2.8.0**, so
  `ghcr.io/cameronabrams/htpolynet:v2.8.0` was never published (the CUDA image
  was).  The base image `condaforge/miniforge3:latest` moved to mamba 2.9.0,
  whose `mamba clean -afy` fails on a symlink inside the CPU AmberTools
  package.  The install was succeeding; only the cleanup step failed.  The
  Dockerfile now cleans with `conda clean`.

### Documentation

- **The installation page now recommends the container as the default way to
  run htpolynet**, in a callout at the top with the two commands you need on a
  desktop (Docker Compose) or on a cluster (Singularity/Apptainer), plus the
  note that the default image is CPU-only and GPU users want the `:cuda` tag.
  Installing htpolynet and its tools natively is still documented below it,
  for development and for systems without containers.

- **The documentation has a logo and a landing-page illustration.**  The
  logo appears in every page's sidebar, in black-on-transparent for the light
  theme and white-on-transparent for the dark one, so it never shows as a white
  box in dark mode.  The browser tab shows it too, as a favicon on a white disc
  so it stays visible whether the tab bar is light or dark.

## [2.8.0] - 2026-09-12

### Added

- **Configuration files are now validated.**  `htpolynet run` checks every
  configuration against a schema before doing anything else, and fills in
  every documented default from that same schema.  Previously a configuration
  was read by picking out the keys htpolynet recognized and silently discarding
  the rest, so a misspelled key did nothing and said nothing: a
  `desired_converson: 0.9` typo produced a default-value build -- a half-cure --
  with no message.  Errors now name the offending location, for example
  `Attribute 'smilse' invalid ... under 'constituents[STY]'` or
  `Attribute 'stage' of 'reactions[1]' must be one of build, param, cure, cap,
  repair`.  Validation covers every section, including molecule records,
  reactions and postcure repair specs; only `gromacs.mdrun_options` is left
  open on purpose.  This adds a runtime dependency on
  [ycleptic](https://pypi.org/project/ycleptic/) 2.4.1 or later.

- **Densification can gate on measured density convergence instead of a fixed
  step count.**  An NPT record in `densification.equilibration` may carry a
  `converge` block; the stage then repeats until the density settles or a
  ceiling is hit.  It is off unless configured, so existing configurations
  reproduce exactly.  The criterion is an autocorrelation-corrected standard
  error, `sigma/sqrt(N/tau_int)` -- an NPT cell density is correlated over
  hundreds of steps, so successive frames are not independent samples and a
  naive standard error is optimistic by roughly 1.5x, which would stop the
  gate early.  It also makes the tolerance size-aware for free.  A small
  standard error alone is not accepted as convergence, because a trace that is
  still climbing looks tight in every short window; the window's two halves
  are compared as well.  **Reaching the ceiling is reported as a failure**,
  with a warning saying the number is an unsettled box rather than the
  system's density.

  Densification is the one place this belongs: the 200 to ~1100 kg/m^3
  compaction happens once and involves a large volume change.  Gating the CURE
  relaxation loop was measured and rejected -- past the gel point the
  unreacted species are bonded into the network and topologically constrained,
  the effective diffusion exponent falls to about 0.14, and restoring the
  mobility a gate would wait for takes of order 10^4 times the relaxation
  time.  Note also that the gate certifies convergence of the Berendsen
  barostat `npt.mdp` currently uses, which does not sample a correct NPT
  ensemble; treat it as a reproducibility criterion, not a physical one.

### Changed

- **A configuration that 2.7.0 accepted may now be rejected.**  This is the
  point of the validation above, but it is a change in behavior: a key
  htpolynet does not know, a value of the wrong type, a value outside its
  allowed choices, or a missing required field now stops the run with a
  message instead of being ignored.  Every shipped example validates.  If an
  existing configuration of yours is rejected, the message says where; the fix
  is almost always a misspelled or obsolete key that was already having no
  effect.

- **Every shipped example now gates its densification on density
  convergence**, using the `converge` block described under Added, so their
  densification stages run until the box settles rather than for a fixed
  duration.

- **The conda-forge package for this release will follow later.**  ycleptic,
  the new dependency, is not yet available on conda-forge, so the conda-forge
  build of this version cannot be published until it is.  In the meantime the
  release is available from PyPI (`pip install htpolynet`) and in the container
  image (`ghcr.io/cameronabrams/htpolynet`); conda-forge users will stay on
  2.7.0 until then.

- **The single moleculetype htpolynet writes is now named `whole_system`, not
  `None`.**  An entire build is written as one moleculetype on purpose -- a
  cured network is one covalently connected molecule -- but the default name
  was the literal string `None`, which reads like an unset field or a bug.  It
  sits on exactly the line users inspect after GROMACS warns about
  "inconsistent shifts", so it was sending them looking for a topology error
  that is not there.  Purely cosmetic: GROMACS accepts either, and existing
  topologies that say `None` still read back unchanged.

### Fixed

- **A config that omitted `CURE.relax.increment` crashed at the first
  relax.**  The default was `0.0`, and `_distance_attenuation` derives its
  stage count as `int(maxL/increment)` with no guard, so the build died with
  `ZeroDivisionError` rather than any message about the config.  The
  documented default -- 0.08 -- is now the actual one.  Dragging's `increment`
  stays 0.0, which is a sentinel rather than the same bug: its own 0.0 `limit`
  disables dragging before anything divides.  Every shipped example sets
  `relax.increment` explicitly, which is why this went unnoticed.

### Documentation

- **New page: analyzing trajectories of periodic networks.**  Explains why
  `gmx trjconv -pbc whole` reports "There were N inconsistent shifts" on a
  cured network (the bond graph wraps through the periodic boundaries, so no
  consistent unwrapped image exists), and that this is expected rather than a
  topology error -- on example 3 the count is zero until the network
  percolates, then 2 at 69% conversion and 140 at 95%.  It gives tested recipes
  for visualization, MSD and free volume.  The MSD one carries a real
  warning: `gmx msd`'s default `-rmpbc` tries to make the network whole every
  frame and silently inflates the MSD, by about 2.4x at 100 ps on example 3,
  while still exiting normally; use `-pbc nojump` and then `-normpbc`.  Prompted
  by a user report.

- The `htpolynet analyze` free-volume docs now say that the per-molecule lines
  in `gmx freevolume`'s output describe the whole box as one molecule, and two
  instances of `poststim` now read `postsim`.

- **The docs landing page carries the standard badge row, and the release
  history is gone from the table of contents.**  `release-history.rst`
  duplicated what this file already records from 1.0.8 forward, so it is
  deleted and `changelog` moves to the *end* of the contents tree rather than
  sitting second.  The seven pre-1.0.8 entries it held -- 1.0.7.2 back to
  0.0.1, which this file did not cover -- were migrated here first, so nothing
  is lost by the deletion.  The landing page now shows the same eight badges
  as the README: tests, PyPI, conda-forge, Python versions, license, docs,
  downloads and DOI.

- **Two more `CURE`/`densification` defaults in the config tables were
  wrong.**  `densification.initial_density` is 200.0 kg/m^3, not the 300.0 the
  table gave; and the prose under `min_allowable_bondcycle_length` still said
  "setting it to zero (the default)" after the table had been corrected to -1.
  Found while transcribing the defaults into a machine-readable schema, which
  is the point of doing so.

## [2.7.0] - 2026-09-08

### Changed

- **The per-iteration equilibration now constrains hydrogen bonds only, not
  all bonds.**  The packaged `npt.mdp` and `nvt.mdp` paired `dt = 0.002` with
  `constraints = all-bonds` and set no `lincs_order`, so GROMACS used its
  default LINCS accuracy (order 4).  That is marginal whenever heavy-atom
  bonds are constrained, and for a halogenated monomer it was fatal: a
  fluorinated bisphenol died at `4-cure_equilibrate-npt` in **7 of 7** build
  attempts, at cure iteration 5-7 of ~10, with exit codes 1 and 139.  `LINCS`
  appeared 14 times in that bridge's `diagnostics.log` and **zero** times in
  each of seven other bridges built identically -- a perfect discriminator
  across eight chemistries.  In a melt of pristine, uncrosslinked monomers at
  production density, with no cure at all, step-0 pressure was -1.34e5 bar;
  `h-bonds` took it to -495 bar and the LINCS constraint rmsd from 4.5e-4 to
  7.3e-7.  With the new settings the whole eight-bridge series rebuilt **32 of
  32**, the fluorinated one reaching a bond conversion of 0.900 in all four
  replicates.  `h-bonds` is also the conventional pairing with a 2 fs
  timestep.

  **This changes the sampled ensemble**, so structures and densities from
  earlier versions are not strictly comparable with new ones, and because
  `postsim` inherits `npt.mdp` it changes the production measurement too, not
  just the cure.  The drag and relax ladders are unaffected -- they always ran
  unconstrained at 1 fs.

- **`lincs_order = 8` in the packaged `npt.mdp` and `nvt.mdp`.**  Independent
  of the constraint change: it improves the accuracy of the constraint solve
  without changing which bonds are constrained, so it does not itself alter
  the ensemble.  On the diagnostic melt above, order 8 alone reduced the
  pressure artifact 4.8-fold.

### Added

- **The CURE relax stages now report the density they produce.**  The relax
  ladder is the only above-`Tg` constant-pressure time in a cure -- roughly
  120 ps of it, at the defaults -- and nothing looked at the density it
  produced: `_do_relax` delegates to `_distance_attenuation`, which never
  calls `TopoCoord.equilibrate()`, the only method that traced Density.  The
  per-stage relax table gains a `Density (kg/m3)` column, read from the NPT
  `.edr` each stage already wrote.  Drag is deliberately excluded, since it
  runs under restraints and its density is not comparable.

- **The CURE relax stages now report reactive-species mobility (Varshney's
  criterion).**  After each relax ladder, htpolynet reports the rmsd
  displacement of atoms that still carry an unused reactive site, and what
  fraction of them moved at least one `CURE.controls.search_radius` during
  the window.  Varshney sized the original 40 ps relaxation window on exactly
  this requirement -- that unreacted species diffuse far enough between
  reactions to find new partners -- but the requirement decays over a cure as
  those species are bonded into the growing network, and until now nothing
  reported when a window had stopped satisfying it.  A build where fewer than
  25 % of still-reactive atoms cross a search radius now says so, because its
  later bonds are being chosen from a nearly frozen neighborhood.

  Both reports are pure observation: they change no simulation input and
  gate nothing.  Every failure path is swallowed and logged at debug level,
  so instrumentation cannot fail a build.  A run resuming mid-ladder skips
  the mobility report rather than measuring only the tail of its window.

- **`repair-summary.yaml` now reports the pre-repair bond histogram.**  A new
  `prerepair_bond_counts` key gives how many crosslinkers carried 0, 1, ... up
  to `full_bond_count` bonds *before* repair dismantled any of them.  Repair
  rewrites the topology and the final structure does not record which cap came
  from which ring, so this distribution was previously unrecoverable after the
  fact -- only `n_complete`, its top bin, survived.  It is the statistic that
  makes the independence assumption behind "crosslinker conversion = bond
  conversion cubed" directly testable, and that relationship is now known to
  be wrong in a way that matters: audited against 54 builds, the effective
  exponent runs about 3.5 near a bond conversion of 0.55, 3.0 near 0.73 and
  2.6-2.7 near 0.90, so the deviation from the cube changes sign and no single
  power law fits.  The histogram is zero-filled, so the shape of the summary
  does not depend on the box, and it is cross-checked against the completion
  count on every build -- the two are computed independently, and a
  disagreement now warns.

### Documentation

- **The constraint trap is documented where someone hitting it will look.**
  The failure above surfaces at the equilibration step immediately after the
  relax ladder, so the natural diagnosis is that the relax schedule is too
  coarse.  It is not: the drag and relax ladders run unconstrained at a 1 fs
  timestep and cannot be responsible, and refining them makes matters worse --
  a three-variant array (increment 0.08 to 0.04 to 0.02, plus a double-MD arm)
  confirmed that, with quadrupling the stage count degrading the result.  The
  CURE section of the program-flow page now says which stages are constrained
  and which are not, and says to look at the constraints rather than the
  ladder when LINCS warnings appear at `cure_equilibrate`.

- **Four defaults in the `CURE.controls` table were wrong.**  The table and
  `curedict_defaults` had drifted apart: `radial_increment` is 0.05 and was
  documented as 0.25, `max_iterations` is 100 and was documented as 150,
  `desired_conversion` is 0.5 and was documented as 0.95, and
  `min_allowable_bondcycle_length` is -1 rather than 0 (any value <= 0
  disallows all cycles, so the documented *behavior* was right and only the
  literal was wrong).  `desired_conversion` is the costly one: a user who
  read "default 0.95" and omitted the key got a half-cure.  Every shipped
  example that cures sets it explicitly, which is why nothing caught this.

- **`desired_conversion` now says which conversion it means.**  The
  `CURE.controls` table described it as "target conversion", which every
  reader takes to be the crosslinker conversion an experiment measures.  It
  is the *bond* conversion -- bonds formed over bonds possible -- and the
  crosslinker conversion is lower: at a nominal 90 % cure, an audit of 30
  builds found 75.2 triazines complete per 100, not 90.  The row now says
  so and links to the fuller explanation under postcure repair.

## [2.6.2] - 2026-08-31

### Fixed

- **The too-few-iterations warning was silent on the failure it exists to
  catch.**  `check_iterations_vs_functionality` returned early whenever the
  cure ran at least `f` iterations, so it fired only for the counting
  impossibility `n < f`.  But `n >= f` only makes completion *possible*.  The
  bonds a cure forms are not spread evenly across its iterations -- the
  per-iteration count decays steeply, and it is the last iteration that has to
  supply each crosslinker its final bond -- so a cure that runs exactly `f`
  iterations and spends its last one on almost nothing arrives at the same
  place as one that ran too few.  A build did exactly that: three iterations
  against a functionality of three, target bond conversion reached, last
  iteration spent on 9 bonds, and **one complete crosslinker out of 240**,
  with nothing in the log about it.  The system it wrote is a monomer melt
  that equilibrates and reports a sensible density.

  The check now tests the outcome rather than the precondition, which needs no
  new measurement: by the time it runs, how many crosslinkers reacted at all
  of their sites is exactly known.  It logs that count on every cure and warns
  when the fraction falls below the crosslinker conversion an ideal randomly
  branching network of the same functionality would show at its gel point --
  12.5 % for a trifunctional crosslinker.  That figure is an idealization and
  is used only to decide when to raise the volume; the count itself is exact
  and is reported either way.  The `n < f` case keeps its own message, because
  there the cause is known exactly and worth naming.

  The regression test suite had the bug written into it: a test asserted that
  a run at exactly `f` iterations with no complete crosslinkers stays silent.

### Changed

- **`min_clearance_nm` is documented as a placement outcome, not a tail
  statistic.**  2.6.1's docs said `min_clearance_nm` and `n_below_target` were
  the tail and were what `cap_min_clearance` should be calibrated against.
  The first half is wrong for the same reason the median was: the direction
  search exits at the first direction reaching the target, so whenever it
  succeeds for every cap the worst-placed cap is one that only just cleared,
  and the reported minimum is pinned to the threshold by construction.  Across
  14 independent real boxes it came in at 0.1503 +/- 0.0006 nm against a
  0.150 nm target, and the single box that fell below it was the single box
  with a non-zero `n_below_target`.  The minimum therefore carries nothing
  `n_below_target` does not already say.  `blind_min_clearance_nm`, which
  ranged 0.006-0.048 nm over the same boxes, is the real tail statistic, and
  the `blind_*` fields are what to calibrate the target against.

- **The cube law is documented as an estimate over a measured band, not as a
  floor.**  The docs said crosslinker conversion sits at or above the cube of
  bond conversion in the many-iteration limit.  It does not, and the band over
  which the cube is even a good estimate is narrower than that claim implied.

  *Where it holds*: across 14 trifunctional runs at 1.7-2.7 `f` iterations and
  bond conversions of 0.74-0.90, the crosslinker conversion sits +2.0 % from
  the cube with a standard deviation of 2.9 %, and 3 of the 14 land *below*
  it, to -3.7 %, against a replicate scatter of 0.8-1.8 % -- so the excursions
  below are real and it is a two-sided estimate, good to about 3 %.

  *Below that band*: on 12 further builds at bond conversions of 0.40-0.73 it
  is neither mild nor two-sided.  11 of the 12 fall below the cube, by -19.9 %
  on average, worsening monotonically as the bond conversion falls and
  reaching -93.6 % at 0.40.  Carrying +/- 3 % down to a bond conversion of 0.5
  understates the error by an order of magnitude, in a predictable direction.

  *Above it*: the four runs of the nine-iteration cohort all landed +0.6 % to
  +3.5 % above the cube, which reads like a bound but cannot establish one --
  four points cannot separate a bound from the upper tail of a two-sided
  distribution, and the scatter measured at 1.7-2.7 `f` cannot be carried to
  3 `f` when regime dependence is exactly what is at issue.  The docs say 3 `f`
  is untested rather than either way.

  *What sets the shortfall*: not the iteration count.  Four runs that each took
  exactly three iterations span 15x in ratio-to-cube (0.06, 0.63, 0.80, 0.94).
  It is how the bonds were distributed across those iterations, and
  specifically how many the last one formed, since that is the iteration which
  has to supply each crosslinker its final bond -- per-iteration bond counts
  decay steeply, so the average is not the operative number.  The run at 0.06
  spent its last iteration on 9 bonds against an average of 96; a run at the
  same bond conversion five months earlier, under an older htpolynet, spent its
  own last iteration on 10 and reached the same ratio to three digits.

  *Which variable to read it against*: bond conversion, because that is what
  the shortfall tracks -- r = +0.81, against +0.57 for the iteration count,
  over 26 runs.  The docs are explicit that this is not the mechanism.  Within
  a fixed iteration count the *average* bonds formed per iteration is exactly
  proportional to the bond conversion, so no set of runs sharing an iteration
  count can distinguish those two, and none here does.

- **The placement documentation now quotes real-box numbers instead of
  synthetic ones.**  `cap_min_clearance`'s "demanding default" was described
  from a synthetic sweep that put the blind median at 0.143 nm and had about
  half of caps needing a search.  On real cured boxes at the same heavy-atom
  density the blind median is 0.120 nm -- ~16 % tighter, in the direction the
  sweep's own caveat predicted -- so the search runs for well over half of all
  caps and 37 % of them would have been placed inside 0.10 nm blind.  It
  nonetheless reaches the target for all but 1 cap in 1955, so the default is
  demanding and reachable at once.  `n_preferred_out_of_angle` came in at 6 %,
  which the docs now give as the scale for reading that field: the 90-150
  degree window is not what is sending caps to the search.

## [2.6.1] - 2026-08-27

### Fixed

- **Cap-placement clearance was pinned at the O-C bond length and measured
  nothing** (regression in 2.6.0).  The neighbour set a transferred `-C#N` cap
  was scored against excluded the cap atoms being moved and the hydrogens
  about to be deleted, but not the bridge oxygen the cap bonds *to*.  The cap
  carbon sits at exactly `oc_len` = 0.136 nm from that oxygen in every
  candidate direction, so the reported clearance was 0.136 nm whichever
  direction was scored.  Two independent acceptance builds on different
  topologies both reported `min_clearance_nm` 0.136, median 0.136, and 71 of
  71 caps below the 0.150 nm target -- identical to three decimals, which is
  what a constant looks like when it is being read as a distribution.

  Consequences: the reported clearance carried no information about crowding;
  the default target was unreachable by construction, so the "could not reach
  clearance" warning fired on every cap in every run and
  `n_direction_searched` was always 100 %; and the search, while not dead,
  could only reject directions worse than 0.136 nm rather than pick the
  clearest one.  Placement was degraded, not disabled, and the repair
  chemistry was never affected -- every conservation identity held through
  2.6.0.

  Each cap now drops the atoms it is bonded *through* -- its attachment oxygen
  and that oxygen's aryl carbon -- and nothing else.  Every *other* cap's
  oxygen is a real atom in the way, as is every cap already placed.  Both
  exclusions are for the same reason: a distance fixed by bond geometry is not
  a measurement of how crowded the site is.  The oxygen is the hard case, at
  exactly 0.136 nm in every direction.  The aryl carbon is the soft one, and
  keeping it would have left `median_clearance_nm` saturated at the C-O-C
  geometry -- so a healthy run would still have reported a constant, just a
  larger one, which is the same defect one bond further out.

  Keeping the cap off the ring it hangs from is now an explicit constraint
  rather than a side effect of the metric: the direction search is restricted
  to C-O-C angles between 90 and 150 degrees.  That guards an end the aryl
  carbon never did -- clearance alone is *maximized* by a linear ether, since
  antiparallel puts the cap as far as possible from the rest of the molecule,
  so a pure clearance search drifts toward a geometry no aryl ether adopts.
  Bond lengths are still held fixed and open-space caps still land on the O-H
  vector, so nothing moves for a cap that was already comfortable.

  Because a systematically rejected preferred direction would reproduce the
  exact symptom being fixed here -- everything searched, everything flagged,
  for a reason unrelated to crowding -- repair now reports
  `n_preferred_out_of_angle` separately from `n_direction_searched`.

  The 2.6.0 calibration missed the original bug because it was done on a
  jittered cubic lattice with no bonded attachment oxygen -- the fixture
  lacked the one feature that defeats the fix.  The regression test added here
  supplies it, and asserts the pinning at 0.136 nm that the released code
  produces.

- **Container releases now carry a version tag.**  `docker.yml` pushed only
  the moving tag (`:latest` / `:cuda`) and the commit sha, so
  `docker pull ...:v2.6.0` -- the obvious thing to type, and what the release
  notes imply -- failed with `manifest unknown`.  That error reads like the
  image failed to build rather than like the tag is named something else, and
  it sent a user off to resolve the release through the GitHub tag API by
  hand.  A `v*` push now also tags `:v<version>` and `:<version>`, in both
  spellings because nothing tells a user which convention a registry chose,
  with the `cuda-` prefix for the CUDA variant.  A `d*` tag is a
  Dockerfile-only rebuild rather than a release and gets neither.

  Releases before this one are unaffected and still have no version tag; the
  container docs now say so and say how to resolve one to its commit sha.

### Changed

- **The docs now say which placement fields measure the box and which measure
  the search.**  `repair-summary.yaml` carries nine numbers about cap
  placement and they are not the same kind of thing, which is worth knowing
  before treating any of them as physics.  `min_clearance_nm` and
  `n_below_target` describe the tail, and are what `cap_min_clearance` should
  be calibrated against.  `median_clearance_nm` is a placement *outcome*: the
  direction search stops at the first direction that reaches the target, so
  the number is pulled toward the threshold you set and partly reports it back
  to you.  `blind_min_clearance_nm` and `blind_median_clearance_nm` are the
  crowding statistics -- one fixed direction, no search, no early exit -- and
  are the ones to correlate across a series of runs.

  `cap_min_clearance` is also now described honestly as a demanding default.
  At the heavy-atom density of a cured thermoset it sits slightly above the
  room a typical site has along the O-H vector, so the direction search runs
  for roughly half of all caps.  That is the intent -- the search exists for
  the crowded half -- but a run reporting that most caps needed a search is
  the default working, not a symptom.

- **The docs now bound the cube law by what has actually been measured.**
  Crosslinker conversion goes as roughly the cube of bond conversion, but only
  in the many-iteration limit, where proximity lifts real runs a couple of
  percent *above* it -- the search is distance-ranked and a partly-bonded
  crosslinker sits in a bridge-rich neighbourhood, so `completion_bias` gets a
  weak version of itself for free.  Below that limit the page now says two
  things and no more than two.  Under `f` iterations the crosslinker
  conversion is *exactly zero*, which is a counting constraint and holds
  unconditionally.  At or above `f` it falls short of the cube by an amount
  set by how the bonds were distributed across the iterations, not by how many
  there were: two trifunctional runs at three iterations apiece came out at
  6 % and 50 % of the cube-law figure.  An earlier draft of this page claimed
  the iteration count decides where a run lands and that the
  one-bond-per-residue rule accounts for the shortfall; the eightfold spread
  refutes the first, and the rule is worth only about a fifth of the second.
  The remaining mechanism is unidentified and the page now says so rather than
  guessing.

### Added

- **htpolynet warns when a cure ran too few iterations for its crosslinkers to
  complete.**  The bond downselection admits at most one bond per residue per
  iteration, so a residue with `f` reactive sites cannot be fully reacted in
  fewer than `f` iterations -- whatever bond conversion was reached.  A cure
  that hits its target in two iterations therefore leaves *every*
  trifunctional crosslinker incomplete, and a system with no complete
  crosslinker has no junctions: it is a monomer melt that equilibrates,
  reports a sensible density, and looks in every other respect like a cured
  network.

  Nothing else surfaces this.  The conversion the cure reports is a bond
  conversion and is perfectly happy with it; only a `postcure_repair` stage
  would notice, and only if one is configured.  The warning also says that
  `completion_bias` is not the fix -- it changes which residues react, not the
  one-per-residue-per-iteration rule -- because it is the first thing a user
  would reach for on seeing a low crosslinker conversion.

- **Repair now records what blind cap placement *would* have given, on the
  real box.**  The O-H direction is tried first anyway, so its clearance is
  free to compute; logging it turns "blind placement was putting caps in
  overlap" from a claim measured on a synthetic lattice into a number measured
  on the system actually being built.  Reported alongside the achieved
  clearance and carried in `repair-summary.yaml` as `blind_min_clearance_nm`
  `blind_median_clearance_nm`, and `n_blind_would_overlap`, with
  `median_clearance_nm` alongside them -- a minimum says whether one build is
  in danger, but correlating placement quality against outcomes across a
  series of runs needs a statistic that is not an extreme value.  This is what tells someone whether structures
  built before v2.6.0 were shaped by overlap resolution rather than by
  placement -- a question that could not be answered from those runs' own
  output.

## [2.6.0] - 2026-08-26

### Fixed

- **Transferred `-C#N` caps are no longer placed blind, which is what has been
  killing low-conversion builds.**  The repair stage relocates a cap by
  putting it along the bridge oxygen's old O-H vector, with no regard for
  what is already there.  In a box at polymer density that direction is
  usually occupied: in a synthetic 13000-atom box at 93 atoms/nm³, blind
  placement put 91 of 158 caps within 0.10 nm of a neighbour and the worst at
  0.008 nm -- effectively superimposed, which is the step-0 Lennard-Jones
  term of order 1e15 that no minimization recovers from.

  The cap is now still tried along the O-H vector first, so a cap in open
  space lands exactly where it always did; only when that direction is
  occupied does the driver search the sphere for the clearest one, holding
  both bond lengths fixed so only the orientation moves.  Caps already placed
  are part of the neighbourhood the next one sees.  On the same synthetic box
  no cap lands closer than 0.128 nm, at either 158 or 378 transfers.

  This matters more the lower the conversion, because the number of caps to
  place is an identity -- `total reactive sites - bonds formed` -- so it rises
  as conversion falls.  `completion_bias` does not reduce it.

- **The stage now says how many fragments it transferred, and complains about
  the ones it could not place well.**  The transfer count is the quantity that
  predicts whether repair survives and it was buried in one mid-log INFO line;
  it is now reported next to the conversion, along with the tightest placement
  achieved.  A cap that could not reach `cap_min_clearance` (new, default
  0.15 nm) in any direction is named, with its oxygen, while that is still
  cheap to act on -- the alternative was working back to cap placement from a
  GROMACS internal error at step 0.

### Added

- **`completion_bias` now says so when it is ranking the wrong side of the
  reaction.**  The bias ranks on the `B` reactant, because htpolynet's A2+B3
  idiom puts the crosslinker there.  Declare the crosslinker as `A` instead
  and the bias quietly starts completing the *bridges* -- a different claim
  about which partly-reacted species is a reactive intermediate, and almost
  certainly not the one that was wanted.  The first cure iteration now
  compares the initial functionality of the two sides and warns if the `A`
  side is the more functional one.

  Worth recording why the obvious detector does not work: looking for an
  all-zero bias key catches nothing, because a difunctional bridge
  accumulates reactions perfectly well and the key is non-zero.  What
  separates the two cases is which side carries more reactive sites, and that
  is readable straight off the atom dataframe -- forming a bond decrements an
  atom's `z` and increments its `nreactions` in the same operation, so their
  sum is conserved and reads the same at iteration 0 as at the end.

## [2.5.0] - 2026-08-26

### Added

- **The conversion a cyanate-ester run reports is now the conversion the
  structure actually carries.**  The cure iterates on bond conversion --
  bonds formed over bonds possible -- and that is the only number it printed.
  But `postcure_repair` then dismantles every crosslinker that did not fill
  all of its sites, so the structure leaving the repair stage contains only
  *complete* crosslinkers, and the fraction of those is what an experiment
  measures.  For a trifunctional crosslinker under random placement the
  second number is roughly the cube of the first, so a run at a bond
  conversion of 0.90 leaves a cyanate conversion near 0.73 -- and nothing
  said so.  The repair stage now logs both figures side by side and writes
  them, with the counts behind them, to `repair-summary.yaml` in the repair
  directory.

- **`CURE.controls.completion_bias`, an opt-in change to how bond candidates
  are ranked.**  Candidates have always been ordered by pair separation
  alone.  Separation is uncorrelated with how many bonds a crosslinker
  already carries, and the downselection that follows admits at most one bond
  per residue per iteration, so bonds spread evenly across every crosslinker
  in the box instead of finishing any of them.  With `completion_bias: true`
  the number of bonds already on the candidate's `B`-side residue becomes the
  primary sort key and separation the tie-break within each group, so a
  crosslinker with two of three sites filled is completed before an untouched
  one is started.  Nothing else about the search changes: same radius growth,
  same dragging and relaxation, same probability application, same cycle
  handling.

  This is a modelling option, not a fix, and it is **off by default** so that
  every existing config and every shipped example behaves exactly as before.
  Where a partly-reacted crosslinker is a stable species the old ranking is
  the more faithful one; where it is a reactive intermediate -- a
  cyclotrimerizing cyanate ester, whose triazine ring either closes or does
  not -- the new one is.  Turning it on changes what `desired_conversion`
  means physically: a config that reaches a crosslinker conversion of 0.76 at
  `desired_conversion: 0.90` unbiased reaches about 0.90 biased, so runs
  either side of the setting must be compared at matched crosslinker
  conversion, not matched `desired_conversion`.

  It also reduces the repair stage's dismantling work: at a bond conversion
  of 0.90 the driver dismantles roughly 24 crosslinkers instead of 58, and
  places 72 caps instead of 174.

  **Corrected 2026-08-26, after this section was published:** the original
  text here claimed the bias also relieves the cap-placement blowup that has
  killed builds at low conversion.  It does not.  Only *transferred* caps are
  placed geometrically -- `_place_cyn_along` runs solely where
  `bonded_o is None` -- and the number of those is an identity,
  `total_sites - bonds`, with no term for how the bonds are distributed.  So
  the bias cannot change it at a matched bond conversion, and at a matched
  crosslinker conversion it makes it *worse*, because it reaches that
  conversion by forming fewer bonds and every bond not formed is a fragment
  that must be transferred.  The blowup is a separate defect in
  `repair/cyanate_cap.py`; see `ROADMAP.md`.

- **A CUDA image, published as `ghcr.io/cameronabrams/htpolynet:cuda`.**  The
  only image until now installed conda-forge's default linux-64 Gromacs,
  which is an OpenCL build; Gromacs no longer drives NVIDIA devices through
  OpenCL, so that image cannot use a GPU at all.  The new tag installs
  `gromacs=*=nompi_cuda*` from the same Dockerfile and can.

  The CPU image remains `:latest` and remains the default, because the CUDA
  build pulls in the CUDA toolkit and inflates the image substantially --
  `docker run` on a laptop should not download a toolkit it cannot use.
  Both tags are built by the same workflow on the same triggers, and both
  carry the per-commit tag scheme (`:<sha>` and `:cuda-<sha>`).

  `external.software.gpu_unusable_reasons()` needed no change: it already
  reconciled detected hardware against what the gmx build can drive, so the
  CUDA image simply starts passing checks the CPU image fails.

  Verified on hardware, since nothing in CI can check this: on a Picotte
  V100 node the image reports `GPU support: CUDA` where `:latest` reports
  OpenCL, `htpolynet info` detects the device instead of calling it
  `unusable`, and a complete `fetch-example 1` build ran every one of its
  `gmx mdrun` steps with `-gpu_id 0`, using cuFFT for PME.  The CUDA 12.9
  runtime in the image runs against that node's 12.4-era driver (550.127.05)
  under CUDA minor-version compatibility, so a site does not need a
  bleeding-edge driver to use the tag.

### Changed

- **Repair drivers now return a statistics dict rather than an operation
  count**, and `htpolynet.repair.run_repair` returns `(total, stats)` rather
  than a bare total.  This is what carries the crosslinker-conversion figures
  out to the runtime for reporting.  Only affects code that calls a repair
  driver directly; the `postcure_repair` config surface is unchanged.

### Fixed

- **The container documentation told users to give the image a GPU, which it
  cannot use.**  `container-usage.rst` carried a "GPU support" section
  walking through the NVIDIA Container Toolkit and a `deploy.resources`
  block reserving nvidia devices, ending "htpolynet will detect the
  available GPU(s) automatically at startup" -- and then, sixty lines later
  in the Singularity section, correctly warned that the image's conda-forge
  Gromacs is an OpenCL build and cannot drive NVIDIA devices at all.  The
  page contradicted itself, and `README.md` repeated the wrong half with a
  `docker run --gpus all` example.

  Both now say the same thing: exposing a GPU to this image starts the
  container and changes nothing about how it computes, so target CPU
  partitions and do not hold a device another job could use.  The reason and
  the detection behavior are stated once, under Docker, and the HPC warning
  points at it for the `--nv`/`--gres=gpu` specifics.

### Changed

- `htpolynet setup-claude` is now discoverable where a new user will meet
  it: the installation page and the README, rather than only the subcommand
  reference.

## [2.4.0] - 2026-08-25

### Added

- **`htpolynet setup-claude` installs the bundled Claude Code skill**, so it
  reaches users who `pip install` or `conda install` htpolynet rather than
  only those working inside a clone.  The skill previously lived at
  `.claude/skills/htpolynet/`, which the tool finds only when the working
  directory is the repository -- that is contributors, and not most users.
  It now ships as package data and is copied to
  `~/.claude/skills/htpolynet/SKILL.md` on request; `--skills-dir
  ./.claude/skills` scopes it to one project instead, and `--force`
  overwrites an existing copy after an upgrade.

  Nothing happens at install time: installing the package never writes to
  `~/.claude/`.  The repository's `.claude/skills/htpolynet/SKILL.md` is now
  a symbolic link to the packaged file, so a clone and an install get the
  same skill and the two cannot drift apart.

  The skill itself was rewritten to stand alone.  It used to be a router --
  its first instruction was to read
  `docs/source/user-guide/building-a-system.rst`, a path that does not exist
  for an installed user -- so it now carries the procedure inline and cites
  Read the Docs once as the full reference.

### Fixed

- **The container image now reports the commit it was built from.**  A
  published image had no way to say what code it contained: there is no git
  in the image and no `.git` beside the installed package, so
  `htpolynet info` fell back to the installed version -- which is right for
  pip and conda and *misleading for the container*, because the weekly
  scheduled rebuild builds from `main` HEAD and reports whatever
  `pyproject.toml` last said.  An image built four commits past a release
  claimed to be that release, and `:latest` is the default thing people
  pull.  A user could pull `:latest`, trust the version string, and record
  the wrong version in a methods section.

  The build now passes the commit as a `HTPOLYNET_COMMIT` build argument,
  the image carries it in its environment, and `htpolynet info` reports it
  in preference to the version fallback.  A real git checkout still wins
  over both, since it reflects the working tree including uncommitted
  changes.  The container-usage guide explains why `:latest` moves and how
  to pull by digest or per-commit tag when provenance has to be stateable
  later.

- The container-usage guide now distinguishes two habits that are easy to
  conflate: **pulling once** gives a campaign a constant tool chain, while
  **recording the digest** is what lets you state afterwards what that tool
  chain was.  The image pins more than the htpolynet code -- Gromacs and
  AmberTools come unpinned from conda-forge at build time, so two images
  built days apart can carry different versions of either while running
  identical htpolynet code.

## [2.3.1] - 2026-08-25

### Fixed

- **`htpolynet postsim` and `htpolynet analyze` no longer die in a directory
  without a `lib/`.**  Both took `-lib` with a default of `'lib'`, but only
  `htpolynet run` creates that tree, so the default value of a flag the user
  never typed reached `UserLibrary` and tripped its existence assertion:
  `AssertionError: lib is not a directory`, three seconds in.  The documented
  container invocation failed on it -- running `postsim` against a staged
  directory holding only finished results, which is the normal shape of a
  cluster workflow, where the cure ran somewhere else.

  `-lib` now defaults to nothing for these two subcommands, and the
  conventional `./lib` is used when it happens to be there.  An explicitly
  supplied `-lib` is still validated, so a typo fails loudly rather than
  silently degrading to no library.  The assertion message now names the flag
  and the resolved path instead of echoing a bare `lib` the user never
  supplied.

### Changed

- Documented that **the density at the end of `postcure` is not an
  equilibrated density**, in the `postcure` directive reference and at the
  point in the example 6 tutorial where a reader would take a number off the
  plot.  A crosslinked network relaxes only above its glass transition, and
  a postcure anneal peaking near *Tg* spends almost no time where the
  network can move, so the plateau reports an under-relaxed structure
  however long the plateau runs.  Measured on four independent BPA builds:
  the plateau gives 1.1712 +/- 0.0029 g/cm3 against 1.1983 +/- 0.0037 for
  the same systems melted and slowly re-cooled -- 2.31% apart, with the
  plateau 2.6-2.8% below experiment and the re-cooled value within 1%.
  Nothing is computed incorrectly; the protocol simply does not equilibrate
  what a reader would assume it does.  The example's own anneal peak is
  unchanged pending a test that raising it actually helps; see `ROADMAP.md`.

- `htpolynet info` no longer reports `HTPolyNet git commit: unknown` for an
  installed copy.  pip, conda and the published container have no `.git`
  beside the package, so the lookup could never succeed and a build had no
  way to identify its own code from inside itself -- the same question the
  parameterization records answer for molecules.  It now falls back to the
  installed distribution version, which for a released install maps to a tag.

- The warnings about a cached parameterization with no provenance record
  now say that the cached values are being **used**, not merely that they
  could not be checked.  The old wording ("cannot be checked against the
  requested 'bcc' charge method") left the decisive fact implicit: the
  build proceeds with whatever charges that entry holds.  A user who
  upgrades with an existing library, asks for `bcc` and measures `gas`
  numbers would reasonably conclude the 2.3.0 cache guard does not work --
  correct behavior is otherwise indistinguishable from the bug it replaced.
  Both the per-molecule line and the stage-end block now state that the
  build carries the cached charges, that they may not match the requested
  method, and that this is expected for a pre-2.3 entry rather than a
  failure.  The same clarification is made in the user guide and the
  `ambertools` directive reference.

## [2.3.0] - 2026-08-23

### Fixed

- **A cached parameterization is no longer reused for a run that asked for
  different AmberTools directives.**  `molecules/parameterized` is keyed on
  molecule name alone, and nothing in that key reflected the charge method.
  A configuration specifying `charge_method: bcc` therefore reused, without
  warning, a `gas` parameterization checked into the user library under the
  same residue name, logging only `Using cached parameterization for TAZ`.
  Measured on a cyanate-ester triazine, the cached `gas` entry carries ring
  C +0.1185 / N -0.2249 where a real `bcc` run gives +0.6539 / -0.7160 --
  5.5x on the charge of the crosslink node -- and the result was a network
  built with `bcc` on one monomer and Gasteiger on another, with nothing in
  the output recording it.  The failure was silent and invalidated results
  without failing the build.

  Each parameterization now writes a `.parm` record beside its
  gro/top/itp/tpx/grx files listing the `charge_method`, `net_charge` and
  `atom_type` that produced it, and that record is checked into the library
  with them.  A run whose directives disagree with the record treats the
  cache as a miss and re-parameterizes, saying which directives differed.

  A library entry written before this release carries no record.  Those are
  still used -- an existing library keeps working rather than
  re-parameterizing wholesale -- but each one logs a warning naming the
  charge method that could not be verified, and the parameterization stage
  ends with a block listing every such molecule and the count.  Rebuild them
  with `--force-parameterization` if you need certainty about what a build
  used.

  Re-parameterizing after a mismatch checks its output in only under
  `--force-checkin`, so a library entry is never silently replaced by one
  built with different directives.

### Added

- `ambertools.net_charge` and `ambertools.atom_type` configuration
  directives, defaulting to `0` and `gaff`.  The net charge was previously
  hardcoded as `-nc 0` in the antechamber invocation, so an ionic or
  zwitterionic monomer was parameterized as though it were neutral with no
  way to say otherwise.  The defaults reproduce the previous commands
  exactly.

- A procedural page in the user guide, **Building a System, Start to
  Finish**.  The guide was strong on reference and had nothing on order of
  operations: start from the nearest bundled example rather than an empty
  file, describe monomers in their active form, check what you can before
  spending compute, and a short list of failures that are known rather than
  mysterious.  It states plainly what `input-check` does not yet verify and
  that there is no seed control, so builds are not reproducible run to run.

- A Claude skill (`.claude/skills/htpolynet/`) for users working in a clone
  of the repository.  It routes to the guide rather than restating it, and
  carries the subcommand table and the handful of facts that are easy to get
  wrong early.

- Test coverage for the parameterization record: 34 tests comparing records
  without external tools, and 6 that run antechamber under both charge
  methods to confirm the record describes what actually ran, that the
  directive it guards changes the charges, and that a library holding a
  `gas` entry refuses a `bcc` request while still reusing it for a `gas`
  one.  The latter skip when the AmberTools chain is absent.

### Changed

- The `-lib` description in the user guide corrected.  It said htpolynet
  would "check-in the results of parameterized molecules ... in
  `lib/molecules/parameterized`", which is false: `-lib` governs lookup
  only, and check-in always goes to the per-user cache at `~/.htpolynet`
  (or `$HTPOLYNET_CACHE`).  A user who set `-lib` expecting to contain a
  run's output was not contained, and nothing announced the writes.
  `HTPOLYNET_CACHE` is now documented as the variable that actually governs
  where products land.

- The user guide's "Parameterization caching" section corrected to match the
  new behavior.  It stated that a stale entry is silently reused when you
  change "SMILES, atom-naming, or charge method" and that the cache is
  "keyed by molecule name only"; the charge-method half of the first and all
  of the second are no longer true.  The warning now covers exactly what the
  record does not: a constituent's structure, atom-naming or geometry
  changing without a rename.

- `--force-checkin`'s help text corrected.  It said "force check-in of
  generated parameter files to the system library", which reads as though
  check-in happens only when the flag is given.  It does not: a molecule
  whose name the library does not yet hold is checked in either way, and the
  flag governs only whether an entry already there is *overwritten*.  It
  also named the wrong library -- check-in goes to the per-user library at
  `~/.htpolynet` (or `$HTPOLYNET_CACHE`), not the system library.  Behavior
  is unchanged.

## [2.2.0] - 2026-08-23

### Changed

- Container docs' HPC section now leads with `htpolynet gen-slurm-script`
  instead of a hand-written batch script -- the subcommand is already
  Apptainer-aware (`--sif`), but was documented only in `usage.rst`,
  so container users had no reason to find it.  The old example
  defaulted to `--gres=gpu:1` and `--nv`, which is actively wrong for
  this image (see the GPU entry under Fixed); replaced with a warning
  explaining why CPU partitions are the right target, plus guidance on
  sizing cores against system size and keeping the submit directory off
  NFS.
- `rdkit` promoted from the `[smiles]` optional extra to a core
  runtime dependency.  Every depot example uses atom-mapped SMILES
  (`[CH:1]`, `[NH2:2]`, etc.), so RDKit is required for any normal
  user workflow; the obabel-only fallback that the extra was
  guarding remains supported but isn't exercised by anything we
  ship.  `pip install htpolynet` (or `uv pip install -e .` from the
  repo) now installs RDKit automatically; `'htpolynet[smiles]'` is
  no longer needed (and is gone from `install.rst`).
- Example 5 (`5-htpb-ipdi.yaml`) retuned for shorter wall-clock.
  The saving comes from a smaller system: the monomer pool drops
  from 125/50/50 to 50/20/20 (IPD/DHT/THT), keeping IPD at
  (2·DHT + 3·THT)/2 so every crosslinker still has both NCO groups
  spoken for at full conversion.  Precure anneal segments go 500 →
  200 ps and postcure postequilibration 1000 → 200 ps.
  Densification was re-balanced for the smaller box in the other
  direction — `initial_density` 50 → 10 kg/m³ to give the long HTPB
  chains room to relax without overlap, and NPT `repeat` 20 → 50 so
  the looser start still reaches target density — so densification
  itself does more work, on a much smaller system.  The build still
  converges.

### Added

- **CI now runs the unit suite** (`.github/workflows/test.yml`), on pushes
  to `main` and on every pull request, across Python 3.10 and 3.13.  Nothing
  ran the tests automatically before, which is how a broken import sat in
  `test_resources.py` aborting collection indefinitely.
- Tests that shell out to `gmx` / `antechamber` / `tleap` / `parmchk2` now
  skip when those binaries are absent instead of failing, so a runner with
  no MD toolchain still reports the ~250 tests that do not need one (in
  under 4 seconds).
- Test coverage for `external/slurm.py` (0% -> 98%),
  `external/smiles_input.py` (0% -> 84%), `utils/inputcheck.py`
  (0% -> 70%), and `analysis/plot.py` (6.7% -> 34%), none of which had
  any.  Overall coverage 34.4% -> 38.8%.  The plot smoke tests were
  checked against the pre-fix module under matplotlib 3.11 and do fail
  there, so they would have caught the `cm.get_cmap` removal.
- Docker image now carries an `org.opencontainers.image.source` label,
  linking the published GHCR package back to the repository.  Without
  it the package is orphaned: it doesn't appear on the repo page and
  doesn't inherit repository-based access permissions.
- `test` optional-dependency extra (`uv run --extra test pytest tests/unit`),
  with `dev` kept as an alias so both spellings work.
- `scripts/run_all_examples.sh`: fail-fast preflight that checks
  every required native tool (`htpolynet`, `antechamber`,
  `parmchk2`, `tleap`, `gmx`, `obabel`, `dot`) is on `PATH` before
  starting any build.  Better than hitting the first missing
  binary hours into a partial run.  Header docstring also gains a
  Prerequisites block pointing at `install.rst` for setup.
- Docker image now built on `condaforge/miniforge3:latest` (was
  `continuumio/miniconda3:latest`).  Miniforge is community-
  maintained, conda-forge only, no Anaconda Inc. terms-of-service
  exposure.  Package installs switched from `conda` to `mamba` for
  faster solves.  Verified end-to-end: built container's
  `antechamber`, `gmx 2025.4-conda_forge`, `obabel`, `parmed`,
  `rdkit`, and `htpolynet 2.1.0` all callable; `htpolynet
  fetch-example 6 && htpolynet input-check` round-trips.
- `install.rst` rewritten around the uv + Miniforge workflow: per-
  repo `uv venv` + `uv pip install -e .` for the Python side,
  separate `mamba create -n gromacs` / `mamba create -n ambertools`
  envs for the native MD binaries (both env bins appended to PATH
  in `.bashrc`).  Documents `uv tool install --editable .` as the
  way to get a global `htpolynet` command callable from any shell.
  Legacy conda-only one-stop install demoted to a subsection.
- `htpolynet.utils.profiling` (moved from `htpolynet.profiling`).
  Small utility module — fits utils/ scope; keeps the package root
  focused on actual subpackages.  Three internal call sites
  updated.

### Fixed

- **API reference documented a module that no longer exists.**
  `docs/source/htpolynetpackage.rst` autodoc'd `htpolynet.driver`,
  removed in the 2.0 refactor, so the package's top-level API section
  rendered empty; it now documents `htpolynet.cli`.  The same page was
  missing nine modules that do exist -- most notably the entire
  `repair` subpackage (`repair.cyanate_cap`, `repair.topology_surgery`,
  i.e. the postcure-repair machinery), plus `external.slurm`,
  `external.smiles_input`, `geometry.lattice`, `utils.profiling`, and
  `utils.vmd_viz`.  A duplicated `htpolynet.core` heading was merged.
- `release-history.rst` appeared in no toctree, so the pre-2.0 release
  history (1.0.7.2 back to 0.0.1, which `CHANGELOG.md` does not cover)
  was unreachable from the docs.  Linked from `index.rst`; its 2.0.0
  date corrected from 2026-04-15 to 2026-05-07 to match the tag.
- Two docstrings (`BondTemplate.matches`,
  `utils.profiling.classify_command`) opened bullet lists with no
  preceding blank line, which docutils rejects; both rendered as
  errors.  `conf.py` also pointed `html_static_path` at a
  `docs/source/_static` that did not exist.  The docs now build with
  zero warnings.
- **The unit suite could not run at all.**  `tests/unit/test_resources.py`
  imported `RuntimeLibrary` from `htpolynet.utils.projectfilesystem` and
  `Software` from `htpolynet.external.software`; neither symbol nor that
  module path survived the 2.0 refactor.  Because the imports were at
  module scope, collection aborted for the *entire* `tests/unit` tree, so
  `pytest tests/unit` had been failing outright rather than reporting
  results.  Rewritten against the current `SystemLibrary` API (15 tests).
- Two `test_chain.py` tests asserted exceptions (`'This is a bug - no
  i-chain!'` / `'no j-chain!'`) that no longer exist anywhere in the
  source: `cure/chain.py` deliberately replaced them with graceful chain
  extension, since bonding to a chain-less atom is legitimate for
  non-vinyl chemistry such as HTPB assembly.  Rewritten to assert the
  current semantics, plus a new test for the `create_if_missing=False`
  branch.  These had been invisible behind the collection failure above.
- `test_write_top` wrote its scratch file into the repository.  The
  autouse `change_test_dir` fixture chdirs each test into a directory
  inside the source tree, and cleanup only ran on the success path, so
  any failure or interrupt left `tests/unit/test_topology/write_test.top`
  behind.  Now uses `tempfile.TemporaryDirectory`.
- **Every plot call crashed on matplotlib 3.11.**  `analysis/plot.py`
  called `matplotlib.cm.get_cmap`, deprecated in 3.7 and removed in
  3.11, at five sites.  Since `pyproject.toml` floors matplotlib at
  `>=3.5` with no ceiling, any reasonably fresh install -- including
  the container image, which tracks latest conda-forge -- died with
  `AttributeError: module 'matplotlib.cm' has no attribute
  'get_cmap'` at the first density trace, i.e. *after* densification
  had already burned its compute.  Replaced with a `_get_cmap()`
  helper that prefers the `matplotlib.colormaps` registry and falls
  back to the legacy call only below 3.5.  Caught by running example 6
  on Picotte through the container.
- GPU usability is now judged on whether the gmx build can actually
  drive the detected devices, not merely on whether its GPU support is
  non-`disabled`.  conda-forge (and hence our container) ships an
  OpenCL Gromacs build; `gpu_ids` is populated from nvidia-smi and so
  only ever lists NVIDIA devices, which Gromacs no longer drives via
  OpenCL.  The previous check passed that combination through, so a
  `gpu_id` from the config reached an `mdrun` that could not honor it.
  New `software.gpu_unusable_reasons()` centralizes the predicate and
  is used by `_mdrun_cmd`, `_enforce_gpu_consistency`, the startup
  banner, and the `grompp_and_mdrun` backstop, which previously
  duplicated a weaker hardware-only version of the test.
- Container image was missing the `graphviz` system package, so the
  `dot` binary `htpolynet.analysis.plot.draw_reaction_dag` shells out
  to was absent.  `pyproject.toml` declares the `graphviz` Python
  binding but the Dockerfile only apt-installed `openbabel` and
  `gosu`.  Failure was silent-ish -- `cure/reaction.py` catches the
  exception and logs `reaction_network.png render failed` -- so
  container builds simply came out with no reaction-network figure.
  Found while porting example 6 to Picotte via Apptainer.
- `htpolynet plots diag` parser templates: the module-path token
  the matcher keyed on was `HTPolyNet.runtime.my_logger` /
  `HTPolyNet.curecontroller.do_iter` from the pre-2.0 namespace.
  After the module reorganization into `htpolynet.core.runtime` /
  `htpolynet.cure.curecontroller`, both lines silently stopped
  matching and the diag parser produced an empty dataframe →
  `IndexError` at first row access.  Templates refreshed to the
  current module paths, and the module-name token dropped from
  `pat_idx` so future renames don't break it again.
- Reaction-network plot (`plots/reaction_network.png`) replaced
  with a bipartite DAG rendered via graphviz `dot` (was a
  spring-layout networkx render that produced tangled, label-
  overlapping diagrams; example 5 was a 30+ node hairball).  Each
  molecule is a rounded box; each reaction is a diamond with edges
  from its reactants and an outgoing edge to its product; nodes are
  colored by role (constituent / intermediate / final) and reaction
  stage (param / build / cure / cap / repair).  Procession-
  expanded reactions (e.g. example 5's `polymerization` with
  `procession.count: 15`, which `parse_reaction_list` explodes
  into 16 sequential reactions + 15 `A18_I*` intermediates) are
  collapsed back into one node labeled `(×N)` so the diagram
  matches what the user wrote.  New runtime dep: `graphviz` (the
  Python wrapper; also needs the system `dot` binary, a separate
  install).

## [2.1.0] - 2026-06-01

### Added

- New `htpolynet.repair` package implementing a postcure topology-repair stage that sits between cure and postcure. Drivers can do bond-breaking, atom deletion, atom transfer between residues, and re-templating — operations the monotonic cure/cap reaction machinery cannot perform. `repair/__init__.py` dispatches each `postcure_repair` config entry by its `type:` field; `repair/topology_surgery.py` provides the generic edit primitives (`delete_bonds` with cascading angle/dihedral/14-pair cleanup, `set_atom_attributes`, `reassign_residue`, `add_bonds_with_template` wrapping `make_bonds` + `map_from_templates` + an int-dtype rescue for atom-index columns that pandas float-promotes via NaN-tainted concat); `repair/cyanate_cap.py` carries the first concrete driver. A new `reaction_stage.repair` enum value lets repair-stage reactions ride the existing symmetry-expansion and parameterization paths so the cure-template lookup at surgery time uses a properly parameterized linked-product Molecule. The runtime gains `cfg.postcure_repair`, `Dirs.systems_repair`, and a `do_repair()` hook wired into `do_workflow` between cure and postcure, including a steepest-descent + short NVT relaxation pass to absorb LJ clashes from relocated cap atoms.
- New `triazine_to_cyanate_cap` repair type: the BADCy-specific driver in `repair/cyanate_cap.py`. At finite cure conversion, the topological A2+B3 BADCy model in example 6 leaves artifacts that don't exist in a real undercured thermoset — free BPA-OH groups and bare triazine C-H sites instead of -O-C#N end-groups. Atom-conservation (the count of unreacted bridge-OH atoms equals the count of dangling crosslinker C atoms across the whole system, exactly) lets the driver dismantle every incomplete triazine (`k < full_bond_count` bonded BPAs) into three independent -C#N fragments via a within-ring C-N matching; the `k` fragments already bonded to a BPA become BPA-O-C#N caps in place, and the remaining `3 - k` are transferred to the nearest unreacted BPA-OH within `cap_search_radius` (greedy matching with radius expansion + global-nearest fallback). After the surgery the heavy-atom neighbors of each deleted sacrificial H absorb its lost charge via `adjust_charges`, keeping the system net-neutral for Ewald. Topology-level outcome on the small test: 19 incomplete TAZ → 57 CYN residues + 1 surviving TAZ + exact heavy-atom conservation, with the C-N bond resolving to GAFF c1-n1 (0.115 nm sp triple) and the BPA-O-C bond to os-c1 (0.132 nm aryl-cyanate ether).
- Example 6 (`6-cyanate-ester.yaml`) rebuilt around the topological A2+B3 + postcure-repair architecture. The BPA-O-C#N cyanate-ester core is now represented topologically: BPA (90 → 360 at 4x scale, two reactive phenolic O atoms) reacts with bare 1,3,5-triazine TAZ (60 → 240, three reactive ring C-H atoms; ring N atoms additionally atom-mapped to N1/N2/N3 so the repair driver can refer to them by name) in a simple cure-stage aryl-ether substitution — no in-cure ring closure, no `bondcycle_collective` bypass needed because the triazine ring is pre-formed in the TAZ monomer rather than constructed via 3-way cyclotrimerization during cure. A new auxiliary CYN building block (`[CH:1]#[N:2]`, hydrogen cyanide; not inserted into the box, exists only as a parameterization template) plus a `repair`-stage `cap_with_cyanate` reaction supplies the auto-generated `BPA~O1-C1~CYN` linked-product template the repair driver splices into the system for every cap. A `postcure_repair: [{type: triazine_to_cyanate_cap, ...}]` block at the end of the YAML drives the conversion. The header comment block explains the topological model, its tradeoff vs. the previous cyclotrimerization model (no cure-kinetics realism, faithful final-network structure), and how the repair stage restores BADCy residual chemistry.
- New `htpolynet/profiling.py` module: a `RunProfile` with a stage-stack context manager (`profiling.stage('name')`) and a subprocess-attribution path. Every external command routed through `external/command.run` (and the two raw `subprocess.run` sites in `external/smiles_input`) is timed and classified — `gmx-mdrun`, `gmx-grompp`, `antechamber`, `parmchk2`, `tleap`, `obabel`, `rdkit`, etc. `do_workflow` wraps each stage (`setup`, `initialization`, `densification`, `precure`, `cure` with one nested frame per `iter-K`, `capping`, `postcure`, `final`) so subprocess time is attributed to whichever stage was active when the call happened. At end-of-run a formatted table is written to the log (one line per `logger.info` call, no `my_logger` asterisk padding) and a machine-readable `proj-N/profile.json` is dumped beside `final.top`.
- New `CURE.controls.min_bonds_per_iteration` knob (default `10`). The bond-search loop now grows the radius until at least this many bonds have been found, falling through to whatever count exists at `max_search_radius`. The effective floor is clamped against `bond_target` (remaining bonds needed to reach `desired_conversion`) and `bond_limit` (the `max_conversion_per_iteration` cap), so demanding e.g. `min_bonds_per_iteration: 50` near end-of-cure never stalls the build. The post-loop "if `nbonds > 0` proceed, else `search_failed`" branch is preserved — accepting fewer bonds than the floor (when max radius is reached) still triggers relax + equilibrate as before. Empirically on the DGEBA/PACM example, `min_bonds_per_iteration=10` cuts the cure iteration count from 41 (with `=1`) to 15; raising further to `=20` saves only one more iteration. The default of 10 was picked off that diminishing-returns curve.
- The "Radius increased to N nm" log line now also reports the cumulative bond count and the iteration's min-bonds floor as `(X/Y eligible bonds so far)`. Makes it visible at a glance whether the floor or `bond_target` is the active constraint as the search radius grows.
- `scripts/run_all_examples.sh` — runs every depot example sequentially in its own subdirectory under `./examples-runs/`. Discovers the example ID list by parsing `htpolynet fetch-example --help` (with a `0..4` fallback). Idempotent: skips `fetch-example` if a YAML is already in the per-example directory. Reports per-example exit status and exits non-zero if any example failed. Pass-through after `--` is forwarded to `htpolynet run`.
- `scripts/run_all_examples.sh --force-reparameterize` — convenience flag that forwards `--force-parameterization --force-checkin` to every `htpolynet run` call. Each example re-runs antechamber/parmchk/tleap on its monomers and overwrites the user cache (`~/.htpolynet/molecules/parameterized/`). Appropriate rigor when consecutive examples share monomers but differ in reaction sets — sidesteps the cache-poisoning interaction between e.g. example 0 (no reactions) and example 1 (cure reactions on STY).
- All five bundled examples (0–4) are now self-contained YAMLs using the RDKit atom-mapping path on each constituent, so the user names reactive atoms by chemical identity (e.g. `[CH2:1][CH3:2]`) instead of by obabel's output ordering. The legacy `.sh` and `.tgz` siblings have been removed; `htpolynet fetch-example N` now delivers a single `.yaml` for any N. Example 4 (DFDA/FDE) additionally gets a `reactive_atoms` entry for `O1`/`O2` that was missing in its prior shell script (the cap reaction references them). `htpolynet fetch-example 1` delivers the YAML directly; usage collapses to `htpolynet run 1-polystyrene.yaml`. The legacy `1-polystyrene.sh` and `1-polystyrene.tgz` have been removed. `fetch-example` now prefers `.yaml` > `.sh` > `.tgz`.
- Final-stage save now emits `final.viz.psf` (real bond topology, written via parmed from `final.top` + `final.gro`) and `final.viz.tcl` (drops any bond longer than 3 Å from the display) alongside the existing `final.gro` / `final.top` / `final.tpx` / `final.grx`. Load with `vmd final.viz.psf final.gro -e final.viz.tcl` to view a crosslinked network without the "long bonds across PBC" artifact. The TCL uses `topo getbondlist both` / `topo setbondlist both $list` — the valid topotools 1.x flag values are `type`, `order`, `both`, `none`; the earlier `all` returned an empty list silently and reported "PSF appears to carry no bonds". The TCL also prints the first bond's measured length so the user can verify VMD loaded coordinates in Å.
- New CLI subcommand `htpolynet make-viz` regenerates `final.viz.psf` + `final.viz.tcl` from any `final.top` + `final.gro` pair without re-running the full workflow. Defaults assume the current directory has `final.top` and `final.gro` (i.e. you've `cd`'d into `systems/final-results/`); override with `-top` / `-gro` / `-prefix`.
- VMD viz now ships a sidecar `<prefix>.viz.macros.tcl` of constituent-keyed `atomselect` macros, sourced automatically from `<prefix>.viz.tcl`. Two layers: `<NAME>` selects every atom of every instance of constituent `<NAME>` (e.g. `GMA` picks all 75 bis-GMAs in example 2; `DHT` picks all 50 linear HTPB chains in example 5), and `<NAME>_<NNN>` picks one specific instance by global molecule index. The macros are built from `final.grx`'s `molecule` / `molecule_name` columns and compress contiguous atom ranges into VMD `index A to B` tokens (so e.g. example 5's 11 KB macros file covers 225 instances). Lets a user highlight chemical entities like bis-GMA or HTPB whose internal residue scheme reflects building blocks (`BPA`+2×`HIE`, `OB`+`TB`×n+`TBO`×2) rather than the assembled molecule. The residue-level view is untouched — `resname TBO` etc. still work — the new macros are additive. `htpolynet make-viz` gains a `-grx` flag (auto-detected from the `-gro` stem) so the macros are also generated when invoked standalone.
- Two follow-ons to the `find_template` bystander relaxation, both needed so the small-fragment cure idiom works end-to-end:
    - `Molecule.idx_mappers` previously asserted that the template and instance had the same bystander count on each side, and built atom-pair mappings by flat-concatenating the per-side bystander lists into one zip. With subset-bystander matching the lists can legitimately differ in length, and the flat concatenation misaligns side-A and side-B bystanders across the zip. Pair each region (bonded residues, side-A bystanders, side-B bystanders, oneaways) in its own zip so a length mismatch in one region doesn't shift the alignment of others, and drop the exact-count assertion.
    - `map_from_templates` copies the template's angle / dihedral / pair tables into the system after mapping template atom indices through `temp2inst`. For cyanate-ester cure, the template's CY has more atoms than the post-build system's CY (the cure-reactive C consumes one H during build, so a system CY has one fewer H on that side than the fresh-from-SMILES cure-template CY does); the template's H atom that doesn't exist in the system maps to NaN. Filter rows whose mapped atom indices contain any NaN before concatenating into the system topology — those rows are force-field parameters for atoms that don't exist in the cured system. No effect on chemistries where all template atoms have system counterparts (examples 1-5).
- `find_template` now uses subset semantics on bystanders. A parameterization-stage template `T` matches a system-instance bond `B` if every (bystander_resname, bystander_atomname) pair declared by `T` also appears in `B`; `B` is allowed to carry additional bystanders that `T` doesn't mention. Oneaway context, atom names, residue names, and the `intraresidue` flag still require exact equality. When multiple templates match, the one with the **most** bystanders declared wins, so chain-extension templates produced by `bondchain_expand_reactions` (which carry specific bystanders) still beat the bare dimer template when their additional context is exactly the in-chain instance's. This unblocks the "small-fragment cure reactant" idiom in cases like cyanate-ester cure, where `CY.C1` is intramolecularly bonded to `BPA.O1` in every BCY-embedded instance — the bare CY+CY dimer template carries no BPA bystander, but the BPA bystander is structural context, not bond chemistry, and the subset rule lets the small template match anyway. Verified: example 1's iter-2 chain-context bond still picks the trimer chain-extension template over the bare dimer template via strict-oneaway discrimination; example 6's iter-1 CY-CY-in-BCY bond now matches the CY+CY cure template via subset-bystander relaxation.
- Fixed: tleap-input ordering in `external/ambertools.py` ran `check mymol` *before* `loadamberparams <frcmod>`, so any GAFF-coverage gap that parmchk2 had already patched (e.g. `h5-ce-n2` on cyanate-ester C=N–C=N dimer templates) still showed up in tleap's output as an early `Error!`. The run-wrapper's override needle then fired and aborted the parameterization, even though tleap actually completed and the `.top`/`.crd` files were valid. Reordered to load the frcmod *before* the check so the patched parameters are in scope when the molecule is validated. Unblocks new depot example 6 (cyanate-ester thermoset) whose C=N–C=N open-chain cure dimer falls in a GAFF coverage gap that parmchk2 patches by analogy. Other examples are unaffected — for chemistry where parmchk2 emits no patches, the reorder is a no-op.
- Fixed: cached monomer `.grx` files in `~/.htpolynet/molecules/parameterized/` carry *reactivity-related* attributes (`z`, `sea_idx`, `bondchain`, `bondchain_idx`) that are YAML-dependent — they reflect the reactions defined for the *run that wrote the cache*, not anything intrinsic to the monomer. Running example 0 (liquid styrene, no reactions) wrote `STY.grx` with all-zero z; running example 1 (polystyrene) afterward then loaded that cache and produced 0 candidate atoms in the cure bond search, silently stalling at "Radius increased to N nm (0/10 eligible bonds so far)" all the way out to the max radius. Fix: in the cache-hit branch of `_generate_molecule`, for monomers (no generator), re-run `initialize_monomer_grx_attributes()` against the current run's `zrecs` so z/sea_idx/bondchain are derived from this YAML rather than inherited from a stale cache. The cache itself can still be written with run-specific z values; only the load-time interpretation is hardened.
- Fixed: corollary of the monomer-cache-poisoning fix. When STY was first poisoned with z=0 by example 0, the subsequent example-1 run generated and cached the cure-stage dimer (`STY~C1-C2~STY.grx`) and the cap (`STYCC.grx`) with *empty* `bondchain` data — the dimer's `chain_manager.injest_bond` no-ops when neither atom is in a chain, which is exactly what happens when the upstream monomer's chain_manager was empty at the time. On the next example-1 run, the cached dimer loaded with 0 chains, `bondchain_expand_reactions` found no 4-atom chains to extend, and zero chain-context templates were generated — so CURE iteration 1 worked (only the dimer template was needed) but iteration 2 raised "you have a bond for which I cannot find a template" because the C1-C2 bond now had a `oneaway` STY chain partner that no available template captured. Fix: after loading a build product from cache, compare the total chain-atom count carried by the cached `chain_manager` against the sum across the product's reactants' `chain_manager`s; if the cache carries fewer atoms (either zero chains, or a partial chain — e.g. example 2's hetero-dimer `STY~C1-C2~HIE` came out length-3 instead of length-4), treat as stale, reset the molecule's `TopoCoord`/`chain_manager`/`bond_templates`/`reaction_bonds`/`sequence` to a blank state, re-parameterize via the normal `generate()` path, and overwrite the cache. The state reset is needed because the cache-load steps populate `TopoCoord` (which `generate()` will then re-merge reactants into) and the half-loaded state ends up float-promoting the `globalIdx` column on the merged dataframe, crashing the prebonding-mol2 writer. Verified on examples 1 and 2: regenerates the affected dimers, after which `bondchain_expand_reactions` produces the expected chain-extension templates (3 for example 1, 32 for example 2 — up from 12 before).
- New depot example `6-cyanate-ester.yaml`: bisphenol-A dicyanate ester (BADCy) thermoset. The BCY constituent is assembled at param-stage from a BPA bisphenol-A core plus two single-carbon `CY` cyanate end-groups (formaldimine, `[CH2:1]=[NH:2]` — drawn in the sp2 imino-formate active form so the cure-stage triazine-forming C-N bonds have one sacrificial H pre-allocated on each side). Mirrors example 2's `GMA = BPA + 2 HIE` build pattern. The cure stage forms C-N bonds between cyanate end-groups on different BCYs via a single `cyclize` reaction; three such bonds among three monomers close into the 1,3,5-triazine ring (the characteristic crosslink of a cured cyanate ester). The `bondcycle_collective` ring-suppression check is C-C-specific via the `ChainManager`, so the heteroatom triazine ring is allowed to close unhindered. Available as `htpolynet fetch-example 6`. Note: pair with `--force-parameterization --force-checkin` when extending the YAML to cover atoms not previously named in any reaction (e.g. CY's N1) — cached build products inherit `zrecs`-derived `z` values from the prior YAML's reactions and won't pick up newly-added reactivity otherwise.
- `write_top` now casts known-int columns (atom indices, function codes, dihedral periodicities, `nrexcl`, etc.) to pandas' nullable `Int64` before serialization, so they emit as e.g. `2` rather than `2.0`. The float form had been silently accepted by `gmx grompp` but rejected by `parmed`'s gromacs top reader, which broke the new `.viz.psf` generation.
- `gmx --version` output is now parsed for `GPU support:` (CUDA, OpenCL, SYCL, disabled) and shown alongside the version line in the startup banner.
- Consistency check: if the YAML config sets `mdrun_options.gpu_id` but the installed gmx was built without GPU support, or no GPU devices are visible on the host, the option is dropped and a warning is logged. This prevents the runtime crash that `gmx mdrun -gpu_id 0` produces when zero devices are detected.
- Cache hits during parameterization are now logged at INFO ("Using cached parameterization for `<name>`") instead of DEBUG, plus a post-loop summary line tallying reused vs freshly-parameterized molecules and a reminder of the `--force-parameterization --force-checkin` flags to invalidate stale entries. Pairs with the new "Parameterization caching" section in the user-guide.
- The `-restart` flag now emits a prominent runtime warning that resumption is experimental and known to fail at the first cure-stage topology update; the argparse help string is annotated likewise, and `docs/source/user-guide/usage.rst` carries an expanded warning explaining the root cause (in-memory cure state is not fully reconstructible from `cure_state.yaml` + on-disk topology files). Pre-cure stages still resume correctly; this section is parked pending a redesign of cure-state persistence.
- Fixed: `htpolynet run -restart` failed in `CureState.from_yaml` with a `ConstructorError` for the `!!python/object:` tag, because curecontroller's loader was `yaml.FullLoader` (which recent PyYAML tightened to reject Python object tags) while the matching `yaml.dump(self)` writes those tags. Switched to `yaml.Loader` to match what `checkpoint.py` already uses.
- Fixed: rebuild `self.chain_manager` from the reloaded coordinates on restart. `do_initialization` is correctly skipped by the checkpoint decorator on resume, but it's also where `chain_manager` was being constructed, so subsequent stages (`do_cure`) hit `AttributeError`. `do_workflow` now reconstructs it from the loaded TopoCoord whenever a checkpoint payload is present.
- Fixed: `htpolynet run -restart` could die with `shutil.SameFileError` when the userlibrary search fell through to `projPath` and the cwd already lived inside it (so the source file IS the destination). `projectfilesystem.py` now uses a `_safe_copyfile` helper that no-ops when src and dst resolve to the same path.
- Fixed: example 2 HIE constituent's SMILES used `[C:1]` (zero implicit H by SMILES bracket-atom rules) instead of `[CH:1]`, so the α-carbon was emitted at valence 3, antechamber typed it `c2`, and tleap failed with "no angle parameter for o - c2 - os". Documented the bracket-atom H-count gotcha in `docs/source/user-guide/molecular-structure-inputs.rst`.
- Fixed: the RDKit SMILES path now goes through an SDF (molfile) intermediate to obabel rather than PDB. PDB does not carry bond orders, so obabel had to re-infer them and frequently mis-assigned a carbonyl carbon as the alkene sp2 type (`C.2` from a `C-O` single bond rather than `C.2`+`O.2` double-bond pair), which propagated to GAFF as `c2` instead of `c` and broke tleap with "no angle parameter for o - c2 - os" on monomers with ester groups (e.g. HIE in example 2). SDF preserves bond orders, so obabel emits the right sybyl types and antechamber assigns the correct GAFF types.
- Fixed: SMILES-generated mol2 files were being written to `projPath/lib/molecules/inputs/<NAME>.mol2` because `Runtime.__init__` runs after `pfs._setup_project_dir` has `chdir`'d into the project directory. `pfs.checkout()` looks in the *user library* (`rootPath/lib/...`), so the files were unfindable and molecule generation fell through to a `BPA.pdb`/`STY.pdb` assertion. `materialize_smiles_inputs` is now invoked with an absolute `inputs_dir` rooted at the user library, and `htpolynet run` pre-creates `lib/molecules/{inputs,parameterized}/` at startup so the library is wired up even when the user's working directory had no `lib/`.
- Constituents in the YAML config may now carry a `smiles:` key.  When present, htpolynet generates `lib/molecules/inputs/<NAME>.mol2` itself before parameterization, eliminating the obabel/sed boilerplate that example shell scripts have historically duplicated.  Reactive atom names are set via either `rename_atoms: {<1-based-index>: <name>}` (obabel path, always available) or `reactive_atoms: {<smiles-map-num>: <name>}` (RDKit path, used when the SMILES contains `[*:N]` atom-mapping tokens and RDKit is importable).  RDKit is an optional extra: `pip install 'htpolynet[smiles]'`; the container ships it by default.

- New `docker-entrypoint.sh` that auto-detects the host owner of the `/work` bind mount and drops privileges via `gosu` before invoking `htpolynet`. Users no longer need to set `--user`, `HOST_UID`/`HOST_GID`, or any other env vars — output files are written with host ownership automatically. The script also writes an `/etc/passwd` entry for the runtime uid so `gosu` resolves `HOME` to `/home/htpolynet` rather than falling back to `/`.
- Entrypoint dispatches by inspecting the first argument: if it resolves to an executable on `PATH` (`bash`, `python`, `obabel`, ...) it is exec'd directly; otherwise it is treated as an `htpolynet` subcommand. This makes it possible to run `docker compose run --rm htpolynet bash 1-polystyrene.sh --run` (i.e. drive the example shell scripts that themselves call `obabel`/`htpolynet`).

### Changed

- Docker-related files moved from the repo root into a new `docker/` subdirectory: `docker/Dockerfile`, `docker/compose.yml`, `docker/docker-entrypoint.sh`. The GH Actions build workflow and the docs `curl -O` URL are updated to match.
- Dockerfile installs `gosu` and uses the new entrypoint script.
- `compose.yml` simplified: no `user:` field; no `HOST_UID`/`HOST_GID` substitution. The entrypoint handles uid mapping at runtime.
- `compose.yml` bind mount uses `${PWD}:/work:Z` so SELinux-enforcing hosts (Fedora, RHEL, openSUSE Tumbleweed, ...) relabel the host directory to `container_file_t`; without this the container is denied writes regardless of POSIX permissions. Harmless on non-SELinux systems.
- `compose.yml` sets `MPLCONFIGDIR=/tmp/matplotlib` in the environment to silence matplotlib's "not a writable directory" warning.

### Fixed

- `_do_pap` previously applied the `_nonempty_directives` guard only to the anneal branch; preequil and postequil ran whenever the section was truthy. That broke for any YAML that defined `precure: preequilibration: ...` without an explicit `precure: postequilibration` — `_apply_runtime_defaults` injects the runtime-default postequilibration block (`ps: 0`, i.e. "no postequilibration"), `TC.equilibrate` returns `None` for `ps: 0`, and the trailing `trace('Density', edr_list, ...)` crashes with `'NoneType' object is not iterable`. Same guard now wraps all three branches, mirroring the existing anneal pattern, with a comment explaining why the postcure-default ps=0 default behaves as it does.
- Example 6 (postcure anneal) first cycle segment changed from `ps: 20` to `ps: 0`, matching example 3's pattern. The earlier value pushed `annealing_time[0]` to 20 ps while `init_t` was 0, tripping `gmx grompp`'s "First time point for annealing > init_t" fatal error at the postcure stage; the leading ps=0 segment exists purely to anchor the annealing protocol at simulation time 0.
- `Topology.rep_ex` was shifting the per-copy `resnr` by `c` (the copy index) instead of `c * residues_per_copy`, so when a multi-residue molecule (e.g. the assembled HTPB chains DHT/THT in example 5, ~41 residues each) was replicated `N` times, the resid ranges of successive copies overlapped almost completely. Single-residue monomers like IPD happened to work because `c * 1 == c`, which is why no prior example hit this. The collision surfaced downstream in `Molecule.idx_mappers` (CURE iter-1 topology update): the resid-keyed filter pulled atoms from several distinct instances at once, the atomName merge fanned out, and the same template atom got mapped to multiple instance atoms, tripping the "temp_idx N already claimed" sanity check. Fix: capture `resnr_per_copy = max(resnr)` before the concat and use it as the shift unit. Existing init.gro/init.top files with corrupted resids must be regenerated (delete `proj-*` and re-run `do_initialization`).
- `TopoCoord.bondcycle_collective` crashed with `AttributeError: 'NoneType' object has no attribute 'added_bonds'` when a cure bond's endpoints weren't part of any vinyl C-C bondchain (e.g. the urethane O-C linkages between TBO oxygen and IPDI formyl carbon in example 5). The chain-manager's `injest_bond` is a no-op when both atoms are outside the chain graph (and `create_if_missing=False`), so the subsequent `chain_of(r.ai)` returns `None`. Such bonds can't form a C-C bondcycle by chemistry — guard the loop and skip them.
- `Molecule.generate_conformers` now reuses existing conformer `.gro` files instead of regenerating them on every invocation. The check is on file presence in the cwd (`molecules/parameterized/`), so it kicks in equivalently on `-restart` and on re-runs in an existing proj dir. Changing `count` upward still triggers regeneration because the expected files won't all be present. For example 5 (HTPB chains use 6 gromacs-generated conformers per stereoisomer per long-chain monomer), this skips a substantial amount of work on every restart.
- Dockerfile: pre-create `/home/htpolynet` with mode `0777` so named docker volumes mounted there inherit a world-writable initial state. Without this, a fresh `htpolynet-home` volume came up owned by root and the non-root container user could not create `~/.htpolynet`, `~/.config/matplotlib`, etc.
- ParmEd `GromacsWarning: The [ pairs ] section contains N exceptions that aren't 1-4 pairs; make sure you know what you're doing!` at end of run. CURE can shorten the topological distance between two atoms previously templated as a 1-4 pair down to 1-3 (a new cure bond shortcuts the original 3-bond path), but nothing was pruning those now-invalid pair entries. New `Topology.prune_stale_14_pairs()` walks the bondlist and drops any `[ pairs ]` entry whose endpoints aren't actually 1-4 in the post-cure graph (uses the existing `bondlist.partners_of`, so it's O(degree³) per pair — fast). Wired into `Runtime.save_data` just before `write_top`; an INFO log line reports the count when any entries get pruned. On the DFA/FDE example 4 build (~23k pairs total) it dropped 40 stale 1-3 pairs.
- Silent no-op config keys in the depot example YAMLs: `initial_search_radius` (schema is `search_radius`) and `late_threshhold` (schema is `late_threshold`, no extra `h`) were being silently swallowed by `dict.get(key, default)` because the values happened to match defaults. Renamed in all four CURE-bearing examples (1, 2, 3, 4). Behavior unchanged but the docs can now quote the keys honestly.
- Silent no-op `nconformers:` keys on FDE and DFA in example 4. The runtime reads conformer settings from a `conformers:` sub-block (with `count`, `generator`, `minimize` keys); the flat `nconformers:` line was never read. Dropped.
- Default `CURE.controls.min_bonds_per_iteration` is now `10` (was effectively `1`, gated by `while nbonds == 0`). This is silent for users with custom YAMLs that don't pin the key; rerun behavior shifts toward fewer, larger CURE iterations.

### Documentation

- Container-usage page rewritten around the entrypoint-driven uid mapping; covers SELinux + `:Z` and the role of the `htpolynet-home` named volume.
- Full rewrite of every short-build tutorial against the current self-contained-YAML workflow. Tutorial directory names now match the depot stem: `2-DGEBA-PACM/` → `3-pacm-dgeba-epoxy-thermoset/`, `3-VE-STY/` → `2-bisgma-styrene-thermoset/`, and a brand-new tutorial 4 (`4-dfda-fde-epoxy-thermoset/`) was written. All cross-references (`ve_*` → `bgs_*`, `dgeba_*` → `pde_*`, new `dfe_*`) renamed accordingly. Each tutorial walks the YAML block-by-block via `literalinclude`, drops the obsolete `run.sh`/`obabel`/`sed` monomer-prep machinery, and points readers at the new `profile.json` and `min_bonds_per_iteration` knob where relevant. Tutorial 0's results page now ships the `final-box.png` VMD render.
- `scripts/run_all_examples.sh` documented in-script via its header block.
- Removed the legacy `src/htpolynet/resources/cfg/` directory (12 orphaned config snippets with no Python references; pre-example-depot artifacts). Updated the one `:download:` reference in `usage.rst` that pointed to a file in there.

## [2.0.1] - 2026-05-12

### Changed

- `compose.yml`: bind-mount source switched from `.` to `${PWD}` so a single shared `compose.yml` referenced via `docker compose -f` mounts the caller's working directory rather than the file's directory.
- `compose.yml`: `user:` field now reads `${HOST_UID}` / `${HOST_GID}` instead of `${UID}` / `${GID}` — bash's `UID` is read-only and `GID` is not exported, so the original form silently fell back to `0` (root) and broke writes into the bind mount under rootless Docker.
- `compose.yml`: container now gets a persistent `HOME` via a named `htpolynet-home` docker volume — without this, `~/.htpolynet` resolved to `/.htpolynet` (root of the container fs) for the non-root user and the user cache could not be created. The named volume also keeps parameterized monomers/oligomers around across `docker compose run --rm` invocations.

### Fixed

- Dockerfile header comments: removed duplicated `htpolynet` token in the example `docker run` invocations (the `ENTRYPOINT` already provides it) and added `--user $(id -u):$(id -g)` so output files are not owned by root.

### Documentation

- Container-usage page now notes that the image is published only to GHCR (a bare `docker run htpolynet` resolves against Docker Hub and fails) and shows a `docker tag` shortcut for a local alias.
- Added a `curl -O` one-liner for fetching `compose.yml` directly from the repo.

## [2.0.0] - 2026-05-07

### Changed

- Package renamed from `HTPolyNet` to `htpolynet` (fully lowercase) for PEP 8 compliance and PyPI consistency.
- Runtime now logs the HTPolyNet git commit hash at startup, with a warning when uncommitted changes are present.

### Added

- Apptainer/Singularity container support: distributed as a `.sif` image for reproducible execution on HPC clusters.
- New `gen-slurm-script` subcommand generates a ready-to-submit SLURM batch script from an htpolynet YAML config file.

### Fixed

- Chain-expansion bug: bond-chain `ChainManager` was not rebuilt for monomers on the fetch path, causing `bondchain_expand_reactions` to produce no chain-extended oligomers in runs that reused cached parameterizations.

## [1.0.9] - 2025-01-01

### Added

- `minimum_bondcycle_length` parameter to allow for cyclic polymerization above a certain threshold length.

### Fixed

- Rings not transferred from monomer templates if they are pre-parameterized.
- Atom indexes in bondchain structure not remapped after atom deletion.

## [1.0.8] - 2024-01-04

### Changed

- Uses `chordless_cycles` to find rings; `ringidx` is no longer a unique atom attribute; improved ring-pierce detection.

## [1.0.7.2] - untagged

### Changed

- Moved the Library package to the `resources` subpackage of `htpolynet`.

## [1.0.6] - 2023-06-21

### Added

- `gmx`-style `analyze` subcommand.

## [1.0.5] - 2022-09-29

### Added

- Post-build MD simulations and plotting functionality.

## [1.0.2] - 2022-09-16

### Changed

- Enhanced molecule-network graph drawing in the `plot` subcommand.

## [1.0.1] - 2022-09-07

### Fixed

- Atom index assignment for systems with more than 100,000 atoms.

## [1.0.0] - 2022-09-03

- First release.

## [0.0.1] - 2022-08-29

- Initial beta version.
