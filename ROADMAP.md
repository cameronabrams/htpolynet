# Roadmap

Ideas worth doing that we haven't done yet. This is a living list, not a
commitment or a schedule. When something here gets done, it moves to
`CHANGELOG.md` and comes off this page.

Rough ordering within each section is by value, not by effort.

## Container and deployment

- **The container's Gromacs is generic `AVX2_256`, on both tags.** The
  `:cuda` image fixes the GPU half of this problem and not the SIMD half:
  conda-forge builds for a portable baseline, so on an AVX-512 host both
  images leave single-core throughput on the table relative to a natively
  built or module-provided Gromacs. This is inherent to installing Gromacs
  from conda-forge and cannot be fixed by choosing a different build string
  -- it would take building Gromacs in the image, which trades the weekly
  rebuild's freshness for a long build and a host-specific artifact.

  Confirmed on hardware 2026-08-25: the `:cuda` image on Picotte's gpu001
  reports `SIMD instructions: AVX2_256` on a node whose CPUs support
  AVX-512, and Gromacs itself prints the hint that "AVX_512 ...
  instructions will perform best on this hardware". Still unquantified, and
  that is what would decide whether it matters: nobody has run the same
  system against Picotte's own `abramsGrp-gromacs/2021.2/cpu-gpu` module on
  the same node. Note Gromacs also observes that AVX2 is often the better
  choice for runs that offload to a GPU anyway, so the penalty may be small
  for exactly the case the `:cuda` image serves.

- **Releases before v2.6.1 have no version tag in the registry.** From
  v2.6.1 on, `docker.yml` pushes `:v<version>` and `:<version>` (and the
  `cuda-` variants) alongside the commit sha, but the earlier releases were
  pushed as `:latest` and a sha only. So `pull ...:v2.5.0` returns `manifest
  unknown`, which reads like a failed build rather than a differently-named
  tag. The docs now say how to resolve an old release to its sha, which is
  the cheap half. Retagging the existing images would be better and is a
  registry operation, not a code change: `docker buildx imagetools create -t
  ghcr.io/cameronabrams/htpolynet:v2.5.0 ghcr.io/cameronabrams/htpolynet:<sha>`
  for each past release, which needs a token with `write:packages`. Worth
  doing if anyone reports hitting it a second time.

- **The image's Gromacs and AmberTools are unpinned, so two images with
  identical htpolynet code can compute different numbers.** `docker/Dockerfile`
  installs `ambertools`, `gromacs`, `parmed` and `rdkit` from conda-forge with
  no version constraints, and the weekly rebuild exists precisely to pick up
  whatever is newest. So the commit stamp now tells you which htpolynet
  produced a build, and still does not tell you which force-field tool chain
  did. Observed live: the image built 2026-08-25 carries AmberTools 26.0 and
  Gromacs 2026.3, while panacea's native environment is on Gromacs
  2025.4 -- a different major version against the same htpolynet.

  Unpinning is deliberate and mostly right: pinning would freeze the image on
  old Gromacs and defeat the point of a weekly rebuild. The gap is that
  nothing *records* what a given image resolved to, so the versions are
  discoverable only by running `htpolynet info` inside it and writing the
  answer down by hand. The better of two candidates is an image
  **label** carrying a `conda list` export, written at build time: `docker
  inspect` then answers the question without running the image, so the record
  survives reaching someone who cannot run the container at all -- a reviewer,
  an archive, a future reader holding only the digest. The alternative, making
  `htpolynet info` machine-readable so a build can capture it, requires the
  image to still be runnable, which is the weaker guarantee. This is the same
  requirement as the build manifest below, one layer further down -- what
  produced this build, all the way to the compilers.

  Raised by the calibration study, which found the 2025.4/2026.3 split only
  because it went looking after an unrelated prompt.

- **Retire `ghcr.io/abramsgroup/htpolynet`.** Superseded by the
  `cameronabrams` package; still public and still serving a June image to
  anyone with an old link.

## Testing and CI

Coverage as of the last measurement: **38.8%** overall.

- **`repair/`'s driver has no test.** `test_cap_placement.py` and
  `test_repair_conversion.py` now cover the placement search and the
  reported statistics, but `triazine_to_cyanate_cap` itself and all of
  `topology_surgery.py` (125 statements) are still untouched. This is the
  highest-value gap: the postcure repair stage makes the strongest
  correctness claim in the project ("atom conservation is exact"), and
  right now the only thing checking it in-tree is reading a residue census
  at the end of a multi-hour build. (An external audit of 54 builds at
  v2.6.2 found the accounting exact everywhere -- `3*TAZ_final + CYN ==
  720` in all 54, every surviving triazine with exactly 3 aryl-ether bonds
  across 7,501 examined, zero bare -OH -- and a six-defect negative control
  was caught by 2-10 checks each. That is real evidence, but it is not a
  test and it does not run on a PR.) It is pure topology manipulation, so
  it can be tested deterministically in milliseconds against a synthetic
  `TopoCoord` carrying triazines at k=0,1,2,3 — assert atom counts, the
  residue census, cap placement, and that no unreacted bridge -OH
  survives.
- **An end-to-end example in CI.** A deliberately tiny build (a
  20-molecule, few-ps variant of example 0) run inside the container
  would cover `core/runtime.py` and `cure/curecontroller.py` — 1,056
  statements, both at 0% — in the only way that is honest, since faking
  the whole AmberTools/Gromacs tool chain to unit-test the orchestration
  is a large effort for less confidence.
- **Remaining zero-coverage modules**: `analysis/postsim.py` (203),
  `cli.py` (165), `analysis/analyze.py` (137), `utils/vmd_viz.py` (75),
  `utils/checkpoint.py` (54).
- **`analysis/plot.py` is at 34%** after the smoke-test pass. The
  diagnostics-log parsers (`diagnostics_graphs`, `_token_match`,
  `_parse_data`) are the part most likely to rot silently — they already
  broke once when modules were renamed — and they are testable against a
  small captured log fixture.
- **Coverage reporting in CI**, so the number is visible on a PR rather
  than something we remember to measure by hand.

## Release and distribution

- **The conda recipe drops `matplotlib`'s version floor.** `pyproject.toml`
  requires `matplotlib>=3.5`; the feedstock recipe lists a bare
  `matplotlib-base` with no lower bound, so conda can solve an environment
  pip would refuse. Found by grayskull's dependency analysis on feedstock
  PR #21, not by us -- `check-conda-sync.py` compared names with the version
  specifiers stripped and so could not see a floor at all (fixed 2026-09-09;
  it now compares floors, and reports this one). Practical risk is low
  because `analysis/plot.py`'s `_get_cmap` carries a pre-3.5 fallback, but
  the declaration should match. **Fold it into the hand-written feedstock PR
  that has to add `ycleptic >=2.4.1` anyway** rather than opening a separate
  one -- same file, two lines.

  Three other items in that same grayskull report are noise and should not be
  "fixed": `python {{ python_min }}` vs `python 3.10` is a Jinja rendering
  artifact, `rdkit >=2024.03` vs `>=2024.3` are the same version to both
  conda and pip, and `ambertools`/`graphviz` are deliberate conda-only
  additions with no PyPI equivalent.

- **The release preflight cannot tell whether the *previous* release
  actually shipped to conda-forge.** `scripts/check-conda-sync.py` compares
  `pyproject.toml`'s runtime deps against the feedstock recipe, which
  catches dependency drift the autotick bot cannot handle. It says nothing
  about whether the bot's last PR ever merged. That gap ran for four
  releases: `pip check` in the recipe's test section started failing when
  `ambertools` began pulling in distributions with unsatisfiable metadata,
  so the bot PRs for 2.2.0, 2.3.0, 2.3.1 and 2.4.0 all sat red and unmerged
  while conda-forge served 2.1.0 from June. Every one of those releases
  passed the preflight, because the deps genuinely did match. Nobody
  noticed until a bot email got read.

  Fixed on 2026-08-26 by dropping `pip check` from the recipe's test section
  (it fails on `ambertools`' bundled distributions, never on ours) and
  merging the 2.4.0 bump, which closed the other three. **The damage is
  permanent and visible**: conda-forge's version list for this package now
  reads 1.0.9 -> 2.1.0 -> 2.4.0, because 2.2.0, 2.3.0 and 2.3.1 were never
  built there and never will be. That is worth knowing as a diagnostic in
  its own right -- a published version list that skips releases the project
  actually made is the retroactive signature of this failure, in any
  package, without needing to have run any check at the time.

  The check is cheap: query `api.anaconda.org/package/conda-forge/htpolynet`
  for `latest_version` and compare it against the version being superseded.
  A mismatch does not have to block the release -- the fix is usually on the
  feedstock, not here -- but it must be loud, because the failure mode is
  silence. Consider also listing open PRs on the feedstock, since a red bot
  PR is the specific thing to look at.

  The deeper point is that publishing to conda-forge is the one leg of the
  release that completes *after* `release.sh` exits and on someone else's
  infrastructure, so it is the only one that can fail without anything here
  noticing.

- **External services still keyed to the old repo identity.** The Aug 2026
  transfer from `AbramsGroup/HTPolyNet` to `cameronabrams/htpolynet` moved
  the code but left every integration pointing at the old owner. Three
  broke and were fixed during the 2.2.0 release: the GHCR package path
  (docs referenced a package that had never been published under the new
  owner), PyPI trusted publishing (`invalid-publisher` — the claim no
  longer matched, so the v2.2.0 upload failed until the publisher was
  re-registered), and Read the Docs, whose project `repository.url` is
  **still** `https://github.com/AbramsGroup/htpolynet`. Builds succeed
  anyway because GitHub redirects the clone, but tag versions never sync —
  `/en/v2.2.0/` 404s and the API reports "No Version matches the given
  query". Fix the RTD project URL, then activate the tagged version. Worth
  keeping this list as the checklist if the repo ever moves again.

  A fourth instance turned up on 2026-08-26 and is fixed: the conda-forge
  recipe's `about:` block still gave `AbramsGroup/HTPolyNet` for `home` and
  `dev_url`, and a `doc_url` of `abramsgroup.github.io/HTPolyNet` that
  returns 404. Corrected in the same feedstock PR that unblocked the version
  bumps.


- **Mint a software DOI.** Enable the Zenodo GitHub integration for
  `cameronabrams/htpolynet`, then the next `scripts/release.sh` run
  archives the release automatically. Afterwards, add the concept DOI (not
  the version DOI) as a README badge and an `identifiers` entry in
  `CITATION.cff`. Note that enabling is not retroactive: releases before
  the toggle are not archived. Also note that adding a `.zenodo.json`
  would make Zenodo ignore `CITATION.cff` entirely — only worth doing if
  we need Zenodo-specific fields such as `grants` for funder linkage.
- **Per-minor Python classifiers.** `pyproject.toml` declares only
  `Programming Language :: Python :: 3`, so the PyPI Python badge reads an
  uninformative `python: 3`. Adding `:: 3.10` through `:: 3.13` would make
  it read `3.10 | 3.11 | 3.12 | 3.13`, which CI already verifies at both
  ends.

## Cure and repair

- **Filtering candidate crosslinks against disallowed topologies, and what that means
  for cyanate esters.**  Raised by Cameron 2026-09-21, who asked the right question:
  can "physically disallowable" even be defined for these systems?  Partly.  Three
  classes, in descending order of how firm the ground is.

  **What is filtered today.**
  - `TopoCoord.makes_shortcircuit` rejects a candidate whose two atoms sit in
    residues that are *already* bonded to each other, so one pair of residues cannot
    be doubly bonded.  Live for every chemistry.
  - `TopoCoord.pierces_ring` rejects a candidate whose segment passes through an
    existing covalent ring, and can dump a gro file of the piercing.  Live, and the
    one genuinely topological test in the code.
  - **The cycle-length filter is inert for cyanate esters.**
    `min_allowable_bondcycle_length` is read into `TC.min_bondcycle_length`, but the
    only code that consults it is the C-C bondchain pass, which exists for vinyl
    chain growth; `makes_bondcycle` itself is commented out in `bondtest`.  For an
    ether O-C bond, `chain_of` returns None and the whole check is skipped.  So
    example 6 and the study's bridge series have run with **no cycle filter at all**,
    and the config key silently does nothing there.  Worth saying out loud before
    designing more filters.

  **Class 1, firm: valence and connectivity.**  Already handled by `z` bookkeeping
  and the short-circuit test.  One gap arrives with the three-body reaction above: two
  -OCN groups of the *same* molecule entering the *same* triazine.  Not forbidden by
  valence, and Lin & Khare forbade the analogous intramolecular loop by fiat.  Make it
  a constraint of the triple search (one group per molecule per ring) with a config
  switch, not a post-hoc rejection.

  **Class 2, firm: catenation and threading.**  `pierces_ring` is one-directional: it
  asks whether a *new bond* pierces an *existing ring*.  Cyclotrimerization creates
  rings in situ, so the complementary test is needed -- does a *newly formed ring*
  encircle an existing strand?  Without it, trimerization can build catenanes that no
  chemistry would produce and no MD will undo.  This is a real, definable test and
  should land with the three-body work.

  **What the existing networks actually contain (htpolynet-study, 2026-09-21).**
  Measured on the v5 builds, 96 cured cells, 1440 bisphenol + 960 TAZ;
  `drops/htpolynet-study-cycle-sizes-20260921.md`, scripts and JSON in the study tree.
  - The atom-level shortest cycle through a cure bond is **quantized at 14n, with a
    floor of 28 in every cell, every conversion, every bridge.**  Nothing shorter
    occurs.  14 atoms is one junction-to-junction repeat (O + 4 aryl + X + 4 aryl + O,
    plus C-N-C through the triazine); a chemistry with a direct biphenyl link would
    give 13, so **14 is not a constant to hardcode**.
  - A 28-ring is two triazines joined by two *different* bisphenols.  Double-bridged
    triazine pairs are **1.4% of bridges, flat across chi = 0.90 to 0.99** (1540
    instances).  The flatness is itself odd -- a diffusion-limited cyclization rate
    should rise as ends get scarce -- and reads as a geometric consequence of random
    crosslinker insertion rather than a kinetic outcome.  Worth explaining before
    filtering these out.
  - Intramolecular loops, both arms of one bisphenol into one triazine, a 14-ring:
    **zero in 96 cells, zero in ~250,000 bridges.**
  - **And the reason is `makes_shortcircuit`, not geometry** -- which turned up
    something worth knowing.  That test keys on the `molecule` gro attribute, which is
    set at template build and replication and **never updated when a cure bond merges
    two molecules**.  So it asks only "are these two monomer instances already
    directly bonded to each other?", which forbids exactly the 14-ring and nothing
    else.  The staleness is load-bearing: were `molecule` the current connected
    component, the test would reject every cycle-closing bond and no system could gel
    (ours carry a cycle rank of ~430 on ~938 junctions).  Now noted in the code at
    both sites that read it, so it does not get "fixed".

  **Class 3, soft: strained cycles, and the threshold problem.**  A candidate bond
  whose two atoms are already close *through the bond graph* closes a covalent cycle.
  Small ones are unphysical -- a cycle threading a few rigid aromatic units cannot
  close without absurd strain -- but "small" needs a number, and inventing one is
  how a filter starts rejecting chemistry that is fine.  The clean generalization is
  a **graph-distance filter**: reject a candidate when the shortest existing bond path
  between its atoms is under N bonds, with N from config, computed by a BFS bounded at
  N on `Topology.bondlist`, which is cheap and chemistry-agnostic.  It subsumes the
  vinyl-only cycle logic conceptually and catches the same-molecule-twice case for
  free.
  - **The data says this filter is insurance, not a correction.**  The floor is a
    28-membered macrocycle -- large and unstrained -- and the one strained candidate,
    the 14-ring, is already excluded.  Because the spectrum is quantized, N has few
    distinct settings in this chemistry: 15 <= N <= 28 is a no-op, N = 29 removes the
    28-rings, N = 43 removes 28s and 42s.
  - **Count crosslinker repeats, not atoms** (the study's ask, and it is right).
    "Reject a bond closing a cycle spanning fewer than k junction-to-junction
    repeats" reads off the junction graph, is chemistry-independent, and avoids a
    per-chemistry magic number.
  - **Ship it off by default** (their second ask).  There is no evidence the 28-rings
    are wrong; intramolecular cyclization is real in step-growth cure and 1.4% is not
    an obviously wrong rate.  The filter's first use is to *measure* the effect of
    removing them, not to correct a defect.
  - **Re-measure for cyclotrimerization.**  These statistics come from the A2+B3
    route, where junctions are seeded at t=0 by random insertion rather than emerging
    where three ends meet, and the cycle spectrum inherits that.  They are a
    measurement of what htpolynet builds this way, not of polycyanurate.  Note also
    that under the three-body reaction the 14-ring is *not* covered by
    `makes_shortcircuit`: both arms of one bisphenol are the same residue, which that
    test asserts against rather than rejects.  The triple search must exclude it
    itself -- the class 1 constraint above.

  **Where the filters have to live.**  Not only in the pairwise search.  The
  single-step annealer chooses a whole connectivity at once, so these become
  constraints on its move set (reject the move, or price it in the cost), and the
  three-body search needs them as constraints on triple selection.  A filter
  expressed only as a post-hoc reject in `bondtest` will not transfer to either.
  Factor them as predicates over (candidate, current topology) that all three callers
  share.

  **Reporting.**  `bondtest` outcomes are already counted per `BTRC` code and logged;
  extend that enum rather than inventing a parallel mechanism, and log the
  per-reason histogram at each iteration.  A filter that silently removes half the
  candidates is worse than no filter.

- **A three-body cure reaction, for the real cyanate-ester chemistry.**  Asked for by
  Cameron 2026-09-21.  Cyanate esters cure by *cyclotrimerization*: three -O-C#N
  groups combine into one 1,3,5-triazine ring, an addition with no atom lost.
  htpolynet cannot express that today -- a reaction forms one bond between two
  reactants -- so example 6 and the study's whole bridge series use a surrogate: a
  pre-formed TAZ triazine monomer plus aryl-ether bonds (A2 + B3), with
  `postcure_repair: triazine_to_cyanate_cap` dismantling the rings that never
  completed back into -OCN caps.

  **Why it is worth doing.**  The surrogate is wrong in ways that cost real work:
  - Triazine rings exist before cure, so crosslink junctions are present at zero
    conversion, and the gel point of the model is not the gel point of the chemistry.
  - Incomplete rings are an artifact with no chemical counterpart, and the repair
    that cleans them up is where the +/-0.24 e charge bug lived (fixed in 2.9.0).
    With real trimerization an unreacted group simply stays -OCN, correct by
    construction, and repair is not needed at all for new builds.
  - Conversion becomes directly what the literature reports.  Today a bond conversion
    has to be mapped to OCN conversion (`chi_OCN ~ chi_bond^3`, exactly zero below f
    iterations -- see the memory note); with trimerization, one reaction event
    consumes three OCN groups and `chi_OCN` is measured, not inferred.

  **What already generalizes.**  `_yield_bonds_as_df` builds its bonds through a
  per-reactant `resid_mapper`, so template *building* is not limited to two
  reactants; `Molecule.generate` merges however many reactants it is given; a bond
  record is a pair of atoms, which stays true for a ring of three bonds.

  **The hard part, and the central new code.**  Per-bond template splicing cannot
  form a ring.  `map_from_templates` matches each new bond by its local context, but
  in the trimer template all three ring bonds exist, while the instance, partway
  through the batch, has only some of them -- so the template's bystanders cannot
  match the instance's, in either direction.  A ring needs a **product-wide splice**:
  map every atom of the three participating residues from the trimer template at
  once -- types, charges, and all angles, dihedrals and 1-4 pairs involving the new
  ring -- rather than three independent per-bond splices.  Treat that as a new
  function beside `map_from_templates`, not a modification of it; CURE's pairwise
  path must keep working unchanged.

  **Plan.**
  1. *Reaction model.*  Let a cure reaction declare three reactants and three bonds.
     Generalize `prepare_new_bonds`, whose residue-offset arithmetic is written for
     exactly two reactants, and the `len(R.reactants) == 2` guards.  Symmetry
     expansion needs care but is bounded: three copies of a dicyanate with two
     equivalent sites give four distinct site combinations, not eight.
  2. *The trimer template.*  One product per site combination, about 120 atoms for a
     real dicyanate: a one-time AM1-BCC cost of tens of minutes, then cached.  Ring
     closure is the question.  Preferred: build it as today by merging and placing the
     three pieces, then close the third bond under a type-6 restraint with the
     existing attenuation ladder before parameterizing -- the drag machinery already
     does exactly this.  Fallback: declare the trimer from SMILES as a
     `count: 0` constituent, which sidesteps the geometry but needs a way to name
     residues inside one generated molecule, which htpolynet has no syntax for.
  3. *Triple search.*  Enumerate candidate triples from the linkcell -- pairs within
     the radius, then triangles among them -- score by the sum of the three
     distances, and choose a disjoint set.  Greedy is the simple version; the
     set-packing problem is close enough to the single-step annealer above that the
     two should share code if both get built.  Forbid two groups of the same molecule
     in one ring, as the single-step entry forbids intramolecular loops.
  4. *Ring closure in the system.*  Three sites have to converge, which is strictly
     harder than bringing two together.  Restrain all three pairs, walk them in with
     the drag ladder, then form the ring in one batch and relax.
  5. *Accounting and checks.*  Conversion counted in OCN groups consumed, three per
     event.  `check_iterations_vs_functionality` assumes a residue carries its own
     crosslink sites; junctions now form in situ from three molecules, so its
     functionality model needs rework rather than a patch.  The per-molecule charge
     check applies unchanged and is worth watching, since a ring closure changes
     types and charges on three residues at once.
  6. *Docs, example, validation.*  A new shipped example, or an ex6 variant, built
     both ways: surrogate versus true trimerization at matched conversion, comparing
     density, Tg, and the fragment and junction statistics, against the study's BADCy
     results.  Builds belong to htpolynet-sweep.

  **Risks.**  Ring closure failing in a dense box is the likely one, and the
  three-way convergence has no precedent in the code.  Beyond that: template cost if
  a system has several distinct cyanate monomers (one trimer template per
  combination), and the fact that changing conversion semantics breaks comparability
  with every existing cyanate build, which is a documentation problem as much as a
  code one.  Out of scope: the dimer and oxazoline intermediates real cyanate cure
  passes through.

- **Single-step crosslinking (the Khare method), as an alternative to CURE, not a
  replacement.**  Agreed with Cameron 2026-09-21.  All bonds are chosen at once by
  combinatorial optimization on one static snapshot, before any bond exists, so
  conversion is 100% by construction; the network is then relaxed into place by one
  restraint-strengthening ladder applied to every bond simultaneously.  CURE stays
  exactly as it is and remains the default.

  **Sources.**  literature's spec, with every parameter attributed to the paper it
  was read from, is
  `~/.local/state/fleet/drops/literature-khare-crosslink-algorithm-20260921.md`.
  The load-bearing ones: `Khare1993Generation` (origin; the open-TSP framing),
  `Lin2009Molecular` (the only complete annealing parameter set, and the viability
  gate), `Khare2012Directed` (the restraint ladder with all k and b0),
  `Jang2015Comparison` (Cameron's own; the acceptance formula, and the only
  head-to-head against multi-step).

  **The method, as specified.**
  - *Stage 1, connectivity by simulated annealing.*  Cost is the total bond length
    over all chosen pairs; the "temperature" is a **length**, starting at 2.5 A;
    accept if the cost drops, else with `exp(-dL/L_C)`; multiply L_C by 0.95 whenever
    attempted moves reach 100x or accepted moves 10x the number of reactive pairs.
    The move changes only two bonds, so a trial costs O(1) -- that is what makes it
    cheap (under 15 min CPU in Lin 2009).  Intramolecular loops are forbidden.
  - *The viability gate, and its loop.*  Annealing returns a sequence, not short
    bonds: the target C-N bond is 1.364 A, van der Waals approach about 3.4 A, and
    most chosen pairs start much further apart.  So before any bond is created, two
    50 ps NVT runs at 300 K under harmonic distance restraints at 8 A then 3 A, then
    a re-scan requiring every pair under 4 A, and **on failure, re-anneal**.  Lin
    2009 needed up to three passes for their nanocomposite.  This loop is the
    difference between a method that works and one that explodes.
  - *Stage 2, directed diffusion.*  Ten stages of 20 ps NVT at 703 K with
    `E = k(b-b0)^2`, k in kcal/mol/A^2 and b0 in A: (1,10) (2,9) (4,8) (8,7) (16,6)
    (32,5) (64,4) (128,3) (256,2) (332.7,1.458).  k doubles each stage; the last is
    the GAFF C-N bond itself.  Then the topology and charges are updated, sacrificial
    H's deleted, and MD run for 4 ns at more than 150 K above Tg.

  **Decisions taken (Cameron, 2026-09-21).**
  1. **Objective:** sum of lengths by default, as the better-documented choice, with
     sum of squares available as a config option.  `Lin2009Molecular` and
     `Jang2015Comparison` use lengths, `Khare2021Atomistic` uses squares, and no
     paper reconciles the two.
  2. **Encoding:** a capacitated matching with a partner-swap move, rather than a
     literal TSP path.  A path whose consecutive pairs are all bonds gives every
     interior site two bonds -- right for a linear chain, wrong for an epoxy carbon
     that forms one.  htpolynet already knows each site's functionality from `z`, so
     the matching formulation handles arbitrary functionality.

     **This is our design decision, not a reading of Khare.**  Lin & Khare 2009
     (journal pp. 4320-4321) name the problem and never state the answer: they
     extended the polystyrene approach "to incorporate the modifications required to
     account for the complexities (e.g., **the number of connections**) of the
     unreacted epoxy monomer, cross-linker, and POSS molecules" -- and what those
     modifications were appears nowhere.  Their prose points both ways ("sequence of
     bonds between randomly picked reacting site pairs" reads as a matching, "a path
     ... that connects all of the potential reacting sites" as a Hamiltonian path),
     but the move settles it: "only the bonds at the two end points of the selected
     path segment get altered" is 2-opt segment reversal on a path, which on a
     matching would change every pair in the segment.  In Khare 1993 there was no
     tension, because polystyrene is linear and an interior monomer really is
     bivalent.  The efficiency argument survives the change: swapping the partners of
     two sites also alters exactly two bonds, so the O(1) delta that makes the
     annealing schedule affordable is kept.  Credit Khare for the method and say
     plainly that the network adaptation is ours.
  3. **Ladder:** ship the 10-stage 2012 table as the default and make it
     configurable.  `Khare2018Quantitative` and `Khare2021Atomistic` use 12 stages
     from a Supporting Information the library does not hold (queued, with Ketan
     Khare's 2013 thesis).  Revisit if it arrives.
  4. **Naming:** the config block is `single_step:`, not a person's name; the docs
     credit Khare and cite the papers.

  **What it reuses.**  Templates, box construction, densification and precure
  unchanged; reactive-site enumeration from the same `reactions:` declarations that
  CURE uses, minus the radius test, so bonds get the same template parameters;
  type-6 restraints between not-yet-bonded pairs (`TopoCoord.add_restraints`), which
  is the same mechanism as CURE's drag; batch bond formation through
  `update_topology_and_coordinates`, which already forms hundreds of bonds per CURE
  iteration; and the whole postcure path.  The second-shell templates added in 2.9.0
  matter here: at 100% conversion every multi-site residue ends fully substituted,
  which is exactly the context they cover.

  **What it needs that does not exist.**
  - `src/htpolynet/singlestep/anneal.py`: cost, move set, acceptance, schedule, with
    no MD in it, so it can be tested on a lattice with a known optimum.
  - A controller mirroring `CureController`: the viability-gate loop, the restraint
    ladder, its own restart state beside `cure_state.yaml`.
  - Restraint parameters set explicitly per stage.  `attenuate_bond_parameters`
    interpolates from each bond's initial distance toward equilibrium; this method
    prescribes (k, b0) per stage for every bond alike.
  - A `single_step:` schema section, and a `do_crosslink` dispatch in the workflow
    that runs either this or CURE.  `do_cure` is not touched.
  - One fix: `check_iterations_vs_functionality` warns when iterations < site count,
    which is meaningless at 100% conversion in a single step.  Skip it there.

  **Validation, and the result to design against.**  `Jang2015Comparison` built
  chemically and stoichiometrically identical systems both ways: thermal and
  mechanical properties and the N-to-N contour length distribution were insensitive
  to the choice, but the fragment molecular-weight distribution after cutting each
  crosslinker was **not** -- single-step gave fewer, larger fragments.  The two
  methods build qualitatively different networks that agree on the properties usually
  reported.  So the acceptance test is example 3 (DGEBA/PACM) built both ways,
  comparing density, thermal and mechanical properties, N-to-N contour length, and
  that fragment distribution.

  **But do not treat the fragment result as a target to hit.**  Because the papers
  never state their network adaptation (see decision 2), that difference belongs to
  whatever connectivity rule their implementation used, not to "single-step" as a
  class.  A capacitated matching may not reproduce it, and would not be wrong for
  failing to.  Read it as the reason to *measure* the fragment distribution both ways
  rather than assume the networks are interchangeable; agreement on density, thermal
  and mechanical properties is the part worth requiring.  Builds belong to htpolynet-sweep.  The per-molecule charge check applies
  unchanged and is worth watching: this forms every bond in one batch.

- **Chain-growth cure leaves whole chains charged.**  The v2.10.0 sweep's new
  molecule-charge check found example 1's polystyrene warning: 7 of 49 chains
  over 0.01 e, the largest -0.040 e on the longest (1410-atom) chain, with the
  system total exactly 0.  htpolynet-sweep measured the same scale in the two
  preceding sweeps (9 of 49, largest -0.022 and -0.025 e), so it predates the
  check, which merely exposed it.  It is about 2-3e-5 e per atom, i.e. a little
  charge per chain-growth bond, which points at the same system-wide
  `adjust_charges` as the item below.  Small enough not to matter much per
  atom, but it scales with chain length, so the longest chain is the worst
  off.  The fix is the per-molecule settle below; what is not yet known is
  whether the residual comes from the templates disagreeing with their
  instances (a `STY~C1-C2~STY` bond in a long chain is not quite the trimer
  template) or only from the system-wide spread.  Measure before fixing.

- **Cure still settles template overcharge system-wide.**
  `map_from_templates` ends each batch of bonds with one `adjust_charges`
  over every mapped atom in the system.  When the templates are consistent
  with each other this does nothing: since the second-shell templates, all 13
  of example 6's molecules after cure are neutral (before them, the
  range was -0.37 to +0.29 e).  But any template that does not describe
  its instance exactly (a mapped atom whose context the template lacks, or
  charges from a different conformer) leaks charge between molecules without
  a trace.  A per-molecule settle, like
  `topology_surgery.neutralize_touched_fragments` in repair, would contain
  it.  Two cases the second-shell rule leaves out: reacted atoms *directly*
  bonded to the bonding atom (it leaves those to the C=C chain machinery), and
  atoms that aren't declared symmetry-equivalent.

- **`cap_min_clearance` has never been calibrated against a working metric.**
  The 0.150 nm default was chosen in 2.6.0 against a clearance that was pinned
  at 0.136 nm by a bug, so it fired on 100 % of caps for a reason that had
  nothing to do with crowding, and the 0.22 nm value rejected before it was
  rejected on the same broken measurement. The metric is real as of the fix,
  but the only two data points on a real box -- the accept260 builds on
  picotte, `/ifs/groups/abramsGrp/cfa22/htpolynet/cyanate-bridge-series/accept260/{unbiased,biased}/`
  -- were taken with the bug present and are uninformative. The geometry now
  bounds the achievable clearance at 0.272 nm (cap carbon to the attachment
  oxygen's aryl carbon, antiparallel), with the chemically sensible ~120
  degree placement at 0.236 nm -- but the aryl carbon has since been excluded
  from the metric precisely so it does not saturate, so that bound no longer
  applies and there is no analytic ceiling to reason from at all. What is
  unknown was the fire rate in a box at polymer density, and that has now been
  measured (see the real-box block at the end of this entry). What is left is
  the recalibration itself: set the default so it flags the tail rather than
  the bulk, using `blind_min_clearance_nm`, `blind_median_clearance_nm` and
  `n_below_target`. Do **not** calibrate against `min_clearance_nm` -- it
  saturates at whatever target is set. A threshold that fires on everything is
  not a threshold. The `n_preferred_out_of_angle` escape hatch is closed: it
  came back at 6 %, so the 90--150 degree window is not biting real O-H
  vectors and is not the thing to change.

  One thing the synthetic sweep did settle, so nobody re-derives it: the
  clipping the aryl-carbon exclusion removes is density-dependent, and at melt
  density it was mild. Median clearance against random neighbourhoods at 5,
  15, 33 and 60 heavy atoms per nm^3 ran 0.316 / 0.230 / 0.193 / 0.177 nm with
  the aryl carbon excluded and 0.233 / 0.222 / 0.191 / 0.174 nm with it in.
  It pins hard at 0.233 nm -- the C-O-C geometry at 118 degrees -- only in
  sparse surroundings, where the caps are comfortable anyway; by 33/nm^3 real
  neighbours are usually closer than the geometric bound and the two agree to
  about 2 %. So the median was censored from above rather than constant, and
  the exclusion buys an uncensored statistic rather than rescuing a dead one.
  Worth knowing before treating a change in that median as physics.

  For the cyanate-ester system specifically, the study session measured the
  box: heavy-atom count is conversion-independent at 7560 (14,040 atoms less
  6,480 H, and the cure removes only hydrogen), post-repair mass 100,192
  g/mol, and measured densities of 1.171--1.198 g/cm^3 across the as-cured,
  ladder-cold-end and monomer-melt states bracket **~53 heavy atoms per
  nm^3**. Repair runs before the postcure anneal so the placement-time value
  is not directly measured, but the spread between those three is smaller than
  the uncertainty over which stage to attribute. That is inside the swept
  range, so nothing here extrapolates. Rerunning the sweep at 53/nm^3 gives
  median blind clearance 0.143 nm (p10 0.075, p90 0.212), median searched
  clearance 0.179 nm, and about 45 % of caps keeping the O-H direction, with
  the censored and uncensored medians differing by under 1 %. Three things
  follow.

  First, and this needs stating at commit granularity rather than as "the
  clearance fix", because the two commits move the numbers by wildly different
  amounts in the same release: 5a58d3d, dropping the attachment oxygen, turns
  `median_clearance_nm` from the constant 0.136 into a distribution, so every
  placement field moves enormously between the accept260 runs and whatever
  comes next. That is the released bug being fixed and is expected. 82ddd3e,
  dropping the aryl carbon, is the <1 % above. Read as one change, the large
  shift becomes evidence about the aryl-carbon exclusion, which it is not.

  Second, a `blind_median_clearance_nm` anywhere near 0.233 in a box believed
  to be at this density is unambiguous evidence the box is wrong, since that
  value is only legitimate in the sparse regime.

  Third, the default sits *above* the blind median: 0.150 against 0.143, which
  is the same fact as only ~45 % of caps keeping the O-H vector. So at melt
  density roughly half of all cap sites cannot reach `cap_min_clearance` along
  the preferred direction and the search is doing real work rather than
  rubber-stamping. That makes 0.150 demanding but not absurd -- which is more
  than it was when it was unreachable by construction -- and it suggests the
  recalibration may move it *down* rather than up. Not the direction either
  session guessed, and worth not being surprised by.

  Independent check on the scale, by the study session: uncorrelated Poisson
  at 53/nm^3 has a median nearest-neighbour distance of 0.1462 nm from a probe
  point, dropping to 0.1236 when minimized over both atoms of a -C#N group at
  0.116 nm separation. The swept 0.143 sits inside [0.124, 0.146], near the
  top, which is where the excluded hole around the attachment oxygen should
  put it. Script at
  `~/devtests/htpolynet/bridge-series/verify_blind_clearance_scale.py`.

  A trap in that check, flagged by the study session and worth disarming
  before someone rederives it: the disjoint-sphere floor,
  `(ln 2 / (8/3 pi rho))^(1/3)`, evaluates to 0.116 nm at rho = 53, which is
  also `cn_len`, the C-N bond length. That is arithmetic coincidence at this
  one density and nothing else. The floor tracks number density -- 0.127 at
  40/nm^3, 0.116 at 53, 0.111 at 60, 0.101 at 80 -- while the bond length is
  fixed. Anyone who writes down "the floor is the C-N bond length" has made a
  claim that is invisible here and wrong everywhere else, and has attributed a
  number-density result to cap geometry.

  Two caveats on those figures. The sweep places uncorrelated Poisson points
  with a hole around the oxygen, which is not a melt -- real packing is
  correlated and has excluded volume between the neighbours too -- so it gives
  a scale, not a prediction, and in particular its `n_below_target` of zero is
  certainly optimistic. And replicate scatter in this system is topological
  rather than numerical: the study session measures density sd rising from
  0.71 kg/m^3 at chi_OCN 0 to 3.25 at 0.746, a factor of 4.6, because at zero
  conversion every build is the same molecular liquid. Cap sites are network
  sites, so expect the scatter in `blind_median_clearance_nm` to grow with
  conversion the same way, and size any series meant to calibrate the default
  accordingly rather than assuming the high-conversion points are as tight as
  the low ones.

  **Measured on real boxes, 2026-08-27** (study session, job 22107390): 14
  independent cured BPA-cyanate-ester boxes at v2.6.1, chi_bond 0.740--0.901,
  1955 caps, ~53 heavy atoms/nm^3 -- the first placement numbers taken without
  the clearance bug present. `blind_median_clearance_nm` 0.1198 +/- 0.0061
  against the synthetic sweep's 0.143: real boxes are ~16 % *tighter*, in the
  direction the sweep's own caveat predicted, so the sweep's absolute values
  should not be used to set the default. `blind_min_clearance_nm` 0.006--0.048.
  722 of 1955 caps (37 %) would have been placed inside 0.10 nm blind.
  `n_preferred_out_of_angle` 109 (6 %). `n_below_target` 1 of 1955, against
  the sweep's 0 that was flagged "certainly optimistic".

  And a correction to this entry's own advice, which the numbers force:
  `min_clearance_nm` is **not** a tail statistic and must not be calibrated
  against. Across all 14 boxes it came in at 0.1503 +/- 0.0006 against a 0.150
  target -- the target seen from above, not a measurement. The search exits at
  the first direction reaching target, so whenever it succeeds for every cap
  the worst-placed cap is one that only just cleared, and the minimum is
  pinned to the threshold by construction. The one box below it (cb0773,
  0.1488) is the one box with `n_below_target` = 1. This is the searched-median
  defect one statistic further along; it was missed because the median's
  version of it was the one being written up. Docs corrected.

  So the recalibration is unblocked and wants the `blind_*` columns. Note that
  0.150 sits well above the measured blind median of 0.120, not marginally
  above 0.143 as the sweep suggested -- so the search runs for well over half
  of all caps -- yet it succeeds for 1954 of 1955. Demanding and reachable at
  once, which is a different situation from the one this entry was written
  against, and it weakens the earlier guess that recalibration would move the
  default down.

  **The placement change is density-neutral at high conversion, 2026-08-28**
  (study session; picotte 22107390 complete, 14/14, exit 0:0). v2.6.1 anchors
  at chi_OCN 0.7438 give 1197.59 kg/m^3 (n = 2) against v2.3.0's 1197.94 +/-
  3.25 at 0.7458 (n = 8): offset -0.34 on a combined SE of 2.38, i.e. 0.1
  sigma, bounding |offset| < 4.7 kg/m^3 at 2 sigma. So the net of 5a58d3d and
  82ddd3e does not move bulk density where it has been checked. Two limits.
  The bound is small against the 21 kg/m^3 rise across the series, so the two
  versions pool for curve shape -- but it is *not* small against the 3.0 kg/m^3
  basin depth, so it licenses nothing about the basin. And it does not
  calibrate `cap_min_clearance`: a placement change being invisible in bulk
  density is a much weaker statement than the threshold being right, and this
  item stays open exactly as written above.

- **The cap direction search stops at the first adequate direction, not the
  best one.** `_choose_cap_placement` breaks out as soon as a direction
  reaches `cap_min_clearance`, so a cap that could have had 0.25 nm of room
  may be left with 0.15. Two consequences. The placement is worse than it
  needs to be, for compute that is genuinely negligible -- the angle window
  admits about 43 % of a 48-direction Fibonacci spiral, so a full scan is
  ~21 vectorized distance evaluations against a neighbour list of order 100.
  And `median_clearance_nm` is partly a readout of the threshold rather than
  of the box, since the distribution is truncated from below at the target --
  as is `min_clearance_nm`, which is pinned *at* the target rather than merely
  pulled toward it, measured at 0.1503 +/- 0.0006 across 14 real boxes.
  Removing the early exit would give both fields their meaning back.

  Deliberately not changed for 2.6.1: taking the best direction moves every
  cap that needed a search, and doing that on the same release as the
  clearance fix means the next real build cannot attribute a change in the
  numbers to either one. That precondition is now satisfied -- the 14-box
  measurement of 2026-08-27 recorded above is the clean pre-change baseline,
  so this is ready to do. Note that changing it invalidates the achieved
  columns of anything built before it while leaving the `blind_*` columns
  comparable, which is another reason those are the ones to correlate on.

- **The shortfall from the cube law between `f` and many iterations is
  unexplained, and is now known to be narrow.** Below `f` iterations
  crosslinker conversion is exactly zero, a counting constraint; the four runs
  of the nine-iteration cohort land +0.6 to +3.5 % above the cube (attested at
  nine for two of the four; the other two are the same cohort at the same
  target with the count unrecorded). In between, measured ratios
  to the cube law are 0.06 and 0.46--0.50 at three iterations and 0.86 at
  four -- an eightfold spread at an *identical* iteration count, so the
  iteration count does not determine it. The band above that is now measured
  and shows no shortfall at all: 14 runs at `n`/`f` 1.67--2.67 sit at
  1.020 +/- 0.029 of the cube, with 3 of the 14 *below* it to -3.7 % against
  replicate scatter of 0.8--1.8 %. So the unexplained region is only
  `f` to about 1.7`f`, and across 1.7--2.7`f` the cube is a two-sided estimate
  rather than the floor the docs used to call it -- corrected there.

  What that does **not** settle is 3`f`, and the temptation to settle it by
  arithmetic should be resisted: applying the 1.7--2.7`f` scatter to the
  nine-iteration cohort gives P(below cube) = 0.245 and makes four-of-four
  above a one-in-three outcome, but that presupposes the deviation is
  regime-independent, which is the claim in question. The cohort on its own is
  +2.64 % sd 1.39 % (n = 4), giving P = 0.03 -- an order of magnitude apart --
  and at n = 4 the sd's 95 % CI is [0.79, 5.18], so the two cannot be
  distinguished. 3`f` is untested, not disproved. Settling it wants more runs
  at 3`f`, not a wider distribution borrowed from a lower one. Whatever
  mechanism explains the shortfall has to switch off within a factor of two
  of `f`.

  The one-bond-per-residue-per-iteration rule was simulated against this and
  moves the n=3 prediction only from 30.0 to 25.9 against 15 observed: a fifth
  of the gap. The leading untested guess is spatial anti-correlation --
  bonds forming preferentially on crosslinkers that already have room, leaving
  the rest to compete. This matters because it is the regime where a user's
  reported conversion is most wrong, and because a warning better than the
  current `iterations < f` one needs a mechanism to threshold on. Widening
  that warning to `n < 2f` was considered and rejected: no threshold on `n`
  alone can work.

  That conclusion survives a correction to its own evidence, and comes out
  stronger. The 8x pair was **confounded**: conv40 r3 and conv50 r1 are both
  n = 3 but sit at chi_bond 0.40 and 0.50, so they never isolated `n` either.
  The new low-conversion builds split it -- at fixed chi_bond the residual
  spread is 1.36x, not 8x, and the chi_bond 0.40 ratio-to-cube reproduces to
  three digits across v2.3.0 and v2.6.1 five months apart. So at fixed `n` the
  ratio moves 7--10x purely with chi_bond, which is a sharper refutation of an
  `n` threshold than the confounded pair was. What does *not* survive is the
  explanation: the "product of per-iteration bond fractions" account is
  superseded by plain chi_bond, which is simpler and reproduces across
  versions. Docs corrected.

  **Reinstated 2026-08-29**, and worth recording because it was withdrawn and
  reinstated inside two days. The "product of per-iteration bond fractions"
  account was not superseded -- it was under-evidenced, and the iteration logs
  now make it the best-evidenced thing in this entry. Bonds per iteration decay
  steeply (bpa-fl0500: 156, 145, 59; bpa-cb0900: 162, 148, 128, 92, 46, 29, 17,
  26), so `720 x chi_bond / n` is an average no build actually realises. What
  matters is the LAST iteration, the one that has to supply a triazine's third
  bond: bpa-fl0400 spent 9 bonds there against an average of 96 and came in at
  ratio 0.064, and conv40 r3 -- same conversion, v2.3.0, five months earlier --
  spent 10 on its own last iteration and also came in at 0.064.

  So chi_bond is a **proxy for the distribution, not a replacement for it**.
  The algebra above still holds for *average* bonds-per-iteration, which is
  exactly collinear with chi_bond at fixed n (r = 1.0000). The *last-iteration*
  count is a third variable at r = 0.980 with chi_bond -- separable in
  principle, not on four points. The withdrawal on 08-28 was made on the
  strength of a correlation, which is the move this entry already warns about
  one paragraph earlier.

  **Measured below the band, 2026-08-28** (study session; picotte 22132281,
  12 BPA builds at v2.6.1, cures final, ladders still running): at bond
  conversions of 0.40--0.73 the deviation from the cube is one-sided and
  large -- 11 of 12 below, mean -19.9 %, worsening monotonically as bond
  conversion falls, -93.6 % at 0.40. The single build inside the documented
  band (0.7306, +1.5 %) agrees with it, so nothing here disturbs the
  1.7--2.7`f` scoping.

  The open question this creates is which variable is doing the work. The
  documented band is *both* 1.7--2.7`f` and bond conversion 0.74--0.90; this
  array varies bond conversion and its iteration counts are not yet in hand
  (they live in `diagnostics.log`, returned only when a task ends -- expected
  early 2026-08-29). So the docs state the degradation against bond
  conversion, which is what was measured, and do not attribute it to `n`/`f`.
  **Resolved 2026-08-29**, counts returned with the array. Four builds sit at
  exactly n = 3 spanning 15x in ratio-to-cube (0.064 / 0.633 / 0.795 / 0.937 at
  chi_bond 0.40 / 0.50 / 0.55 / 0.58), monotone in chi_bond. So "iteration
  count does not determine the shortfall" is established directly, on four
  same-n runs, rather than inferred from the confounded pair. Over all 26
  builds r(ratio, chi_bond) = +0.813 against r(ratio, n/f) = +0.567, with
  r(chi_bond, n/f) = +0.905; counts run 3 to 8, n/f 1.00 to 2.67. Docs cite
  the four-run evidence now.

  **The deviation changes sign, 2026-09-03** (study session, 54 builds
  audited at v2.6.2): near `chi_bond` 0.900 the observed crosslinker count
  runs *above* the cube -- 180.60 against 175.39 predicted over 30 builds,
  +2.17 pp at t = +10.8 against a replicate sd of 2.63 triazines -- while at
  `chi_bond` <= 0.64 it runs below, 32.17 against 41.28 over 6 builds, -3.80
  pp. The compact way to say it is that the effective exponent
  `ln(chi_OCN)/ln(chi_bond)` is not 3 and is not constant: ~3.5 near
  `chi_bond` 0.55, ~3.0 near 0.73, ~2.6--2.7 near 0.90. No single power law
  fits, so the cube is a local estimate around 0.9 and nothing more. This is
  consistent with the 1.020 +/- 0.029 above, and puts a significance on it.

- **A config-time version of that warning, before any compute is spent.**
  Largely superseded by the completed-crosslinker check that shipped -- nothing
  needs a proxy for a quantity that is exactly known by the time the cure
  ends -- and worth keeping only for
  what the post-hoc check cannot do: warn *before* a multi-hour build rather
  than after it. That is a weaker claim on much worse evidence, so it is a
  separate feature and not a substitute. The reasoning that motivated it:
  chi_bond predicts the ratio-to-cube to within 1.36x, and unlike a bond
  distribution -- which does not exist until the cure runs -- it is knowable
  before the cure starts. So
  htpolynet could warn at config time that a run will yield far less
  crosslinker conversion than the cube law implies, which is the number a user
  is actually reasoning with. The existing `iterations < f` warning fires only
  after the fact and only in the exactly zero regime.

  **Do not key it on `desired_conversion`, and do not key it on chi_bond
  either.** chi_bond is a config-time proxy, not the governing quantity, and
  the collinearity is not an accident of sampling that more builds would fix.
  Within a fixed iteration count `bonds = 720 x chi_bond` for this system, so
  bonds-per-iteration `= 720 x chi_bond / n` is *exactly* proportional to
  chi_bond. No dataset at fixed `n` can separate them, however large. Breaking
  the proportionality needs `max_conversion_per_iteration` or
  `min_bonds_per_iteration` varied AT FIXED `desired_conversion` -- which is
  a deliberate experiment nobody has run, not something a wider sweep of
  conversions will deliver.

  There is a further reason to keep it qualitative: the quantity that appears
  to govern the shortfall is the number of bonds formed in the *last*
  iteration, and that is not knowable at config time at all -- it is an
  outcome of the cure, not a setting. Any config-time warning is necessarily
  keyed on a proxy for it.

  Those same two directives are also why `desired_conversion` is the wrong
  key in ordinary use: they set bonds-per-iteration independently of it, which
  is exactly what the warning at `cure/curecontroller.py:674` already tells a
  user to do to spend more iterations at the same conversion, and what this
  project plans in order to reach 9+ iterations. A threshold keyed on
  `desired_conversion` alone would misfire on precisely the runs that took
  htpolynet's own advice.

  So keep it qualitative until something separates the two. A warning saying
  "below about 0.7 the cube law overstates badly, and worsens as you go lower"
  needs no attribution and could ship soon. One that quotes a number needs both
  the separation and a sense of how much of the curve is BADCy-specific -- it
  is one chemistry, one force field, one cure protocol, n = 1 per target apart
  from a single duplicate.

- **`completion_bias` biases the `B` side only, and that is a convention,
  not a law.** The new `CURE.controls.completion_bias` ranks bond candidates
  by how many bonds their `B`-side residue already carries, because
  htpolynet's A2+B3 idiom puts the multifunctional crosslinker in the `B`
  position -- as example 6 does. A user who declares the crosslinker as `A`
  gets the bias pointed at their *bridges* instead: not a no-op, a different
  claim about which partly-reacted species is a reactive intermediate. That
  case is now warned about at the first cure iteration, by comparing the
  initial functionality of the two sides, so it is visible rather than silent
  -- but it is still not *served*. The generalization is small (rank on the
  `A` side, or on the sum of both) and each variant is a different physical
  claim, none of them run. Do it when someone has a chemistry that needs it,
  and make them say which side rather than inferring it from functionality:
  the warning's heuristic is good enough to flag a probable mistake and not
  good enough to silently redirect the ranking on.

- **Nobody has run example 6 with `completion_bias` on.** The ranking is
  unit-tested -- ordering, tie-breaks, missing residues, the fallback when a
  `.grx` predates the `nreactions` attribute -- but the acceptance criteria
  that matter are system-scale and need gmx and a multi-hour build:
  incomplete triazines should collapse from tens to ~1 at
  `desired_conversion: 0.90`, crosslinker conversion should rise from ~0.76
  to ~0.90, `TAZ + CYN/3` must still equal the initial triazine count
  exactly, and atom conservation must still be exact. If the incomplete
  count does *not* collapse, the sort key is not surviving as far as the
  truncation and that is where to look. Until someone runs it, the directive
  is documented as untested at scale.

- **Cap placement is greedy and never backtracks.** Transferred `-C#N`
  fragments are now placed against the local neighbourhood rather than blindly
  along the old O-H vector, which removes the catastrophic overlaps: in a
  synthetic box at polymer density, blind placement put 91 of 158 caps within
  0.10 nm of a neighbour and the worst at 0.008 nm, while the neighbourhood
  search puts none below 0.128 nm. What it does not do is relax. Each cap is
  placed once, in the clearest direction available *at that moment*, and later
  caps must work around it; in a genuinely over-packed box a residue of
  sub-target placements remains (17 of 158, 33 of 378 in the same synthetic
  test) and is reported rather than fixed. The next lever is minimizing
  incrementally as caps land, or placing in order of how constrained each site
  is rather than in match order. Do it if the reported clearance warnings turn
  out to correlate with builds that still die -- the instrumentation to decide
  that now exists, and did not before. One thing the instrumentation still
  does not report is how far a transferred cap actually travelled:
  `_greedy_match` can reach the globally nearest free oxygen after ten radius
  doublings, and one audited build transferred 72 of 177 caps, but
  `repair-summary.yaml` carries clearance statistics only, so an external
  audit could not bound the transfer distance. Report it alongside the
  clearances.

- **`bdf.loc[:abs_max]` takes one bond more than the limit.** In
  `curecontroller.py::_searchbonds`, the truncation that applies
  `max_conversion_per_iteration` uses `.loc` with a slice, which is
  inclusive of its endpoint, so an iteration limited to `n` bonds forms
  `n + 1`. Harmless in practice -- the limit is a throttle, not a
  correctness bound -- but it is off by one, and fixing it changes the
  trajectory of every existing config by one bond per throttled iteration.
  Worth doing at a version boundary where a small reproducibility break is
  already expected, not before.

- **An algebraic self-check that no build currently runs.** For an A2 monomer
  of `A` atoms in the 360-bisphenol/240-triazine example, the total atom count
  is `360*A + 2160 - 1440` regardless of conversion -- verified at 12600
  (bpa, A = 33), 9720 (bpo, A = 25) and 14760 (tmb, A = 39). It is
  conversion-independent and closed-form, so any atom dropped or duplicated
  anywhere in cure or repair breaks it. Cheap to assert at the end of repair.

- **`_distance_attenuation` runs an open-loop stage ladder.** In
  `cure/curecontroller.py:606`, `this_nstages = int(maxL/d['increment'])` is
  computed once from the initial maximum bond length, and the loop then runs
  that many stages whatever happens. `maxL` *is* recomputed inside the loop,
  but only to log it -- it never extends, shortens or aborts the ladder -- and
  `restore_bond_parameters(saveT)` runs each stage whether or not the bonds
  actually closed. So a stage schedule that fails to bring a bond in has no
  way to say so. Robustness only; it was ruled out as the cause of the
  halogen failure below (quadrupling the stage count made that *worse*).
  The densification entry below is the same shape of problem -- a number
  computed once and trusted thereafter.

- **Decided against: gating and extending `CURE.relax` on density
  convergence.** Raised by Cameron 2026-09-05 from the observation that
  pestifer has a gate system htpolynet does not; recommended by the study
  session, twice revised, then withdrawn 2026-09-07 once the experiment that
  tested it had run. Recorded here so it is not re-proposed from the same
  premises. **The verdict is "not now, on these grounds", not "never"** --
  both legs of it have stated limits, recorded below with the evidence.

  **The experiment.** Arm R16 raised `CURE.relax` NPT `nsteps` 2000 ->
  32000, i.e. above-`Tg` relaxation during cure from ~120 ps to ~1.9 ns, on
  two chemistries at matched conversion, n = 4:

      bpa   1187.52 -> 1192.19 kg/m3   +4.67 +/- 2.61   (+1.8 sigma)
      bpf   1257.65 -> 1255.44         -2.21 +/- 2.15   (-1.0 sigma)

  Sixteen times the relaxation closes **17 % of the density deficit, and not
  significantly**; bpf, the control that should not move, did not. The
  premise -- that the fixed relaxation window generates the deficit and a
  gate would recover it -- survives only as a minor contributor.

  **R16 is underpowered for the effect it found, and must not be recorded as
  settled beyond its power.** Resolving +4.7 kg/m3 at 3 sigma needs n ~ 11 at
  the observed scatter; this is n = 4. "17 % of the deficit at +1.8 sigma" is
  a point estimate, not a tight upper bound, and better statistics could move
  it. What it does rule out is a *large* effect of relaxation duration --
  which is what the gate would have had to deliver to be worth its cost.

  **The stronger reason, which does not depend on that result.** R16 ran a
  fixed multiple, never to convergence, so ask what convergence would cost.
  Restoring the start-of-cure crossing fraction needs ~4.4x the displacement
  of unreacted species, which under `sqrt(t)` is ~19x the time: ~37 ns of
  relax per build against a 39.2 ns production ladder, roughly doubling
  build cost. **And `sqrt(t)` does not hold.** Past the gel point the
  unreacted species are bonded into the network -- topologically
  constrained, not merely slow -- and no amount of NPT time restores a
  capture-radius displacement to a covalently tethered fragment. So an
  honestly implemented gate would hit its ceiling on every late-cure
  iteration and report "may not have settled", correctly, with nothing to be
  done about it. It detects a condition it cannot cure.

  **That limit is now closed, and the conclusion is stronger than it was.**
  The argument used to lean partly on gelation physics rather than
  measurement. 32 cure-only builds at the production window (study session,
  2026-09-08) measure the effective diffusion exponent directly:
  displacement ~ t^alpha, with alpha **0.32** at the start of cure -- already
  sub-Fickian, which is 0.50 -- falling to **0.14** at the end. At alpha 0.14,
  restoring the displacement a gate would need takes **3.6e4 times** the
  relaxation time. That is not expensive, it is impossible, and it does not
  depend on the earlier n = 4 result at all.

  **And the shipped window never satisfies Varshney's criterion, not even at
  the start.** Mean crossing fraction is **0.61 at iteration 1** and **0.036
  at the last** -- so the window is outside the criterion it was sized against
  for the whole cure, not merely late in it. This is the calibration for the
  mobility report that shipped in v2.7.0: its 25 % warning threshold sits
  between those two numbers, so a default build stays quiet early and warns
  through the second half. Whether 25 % is the right line is not yet
  established -- nothing measures what crossing fraction a *trustworthy* build
  needs, only what the current one delivers -- so the threshold is a placeholder
  chosen to be quiet at iteration 1, and should be revisited against a build
  whose bonds are known to be well sampled.

  **And a structural mismatch.** Pestifer uses `density_equilibrate` as a
  terminal, one-shot replacement for a hand-written NPT ladder at the end of
  a build -- somewhere you can afford to run until converged. The proposal
  put the same machinery in an inner loop that executes ~10 times per build.
  The economics are entirely different.

  What does *not* change is the mechanism the proposal was built to explain,
  which is still unexplained. If the box is under-relaxed while bonds are
  forming, bonds form between whichever atoms happen to be adjacent in a
  too-open configuration, and **the topology is then permanent** -- a later
  remelt relaxes coordinates but cannot rewire bonds. That still ties
  together three otherwise-awkward results from the study's 64-build series:
  post-cure annealing is refuted as a cause of the density deficit (the
  ladder opens 2 ns at `Tg`+112 K, erasing pre-ladder history), the uncured
  monomer melt matches experimental dilatometry to 0.19 % while the cured
  network is 2.25 % low, and cure shrinkage is strongly bridge-dependent
  (+0.033 to -0.004 ml/g) and anti-correlates with monomer van der Waals
  volume. R16 says more relax time is not the fix; it does not say the
  mechanism is wrong.

  **The production ladder does not need a gate either.** Two structures with
  opposite histories held 5 ns at 480 K converge from opposite sides and
  stall 4.25 kg/m3 apart with the ladder's own hold inside that bracket;
  within-hold drift on the glassy branch is -0.35 to +0.53 kg/m3 per 667 ps
  across four replicates; and replicate scatter is topological rather than
  equilibrative (sd 3.2 kg/m3 at high conversion against 0.7 at `chi_OCN` 0).

  **Reference values for the instrumentation that did ship in v2.7.0.** The
  relax and mobility numbers a user now sees in their own log mean nothing
  without something to compare them against, and these are the only measured
  ones we have. Window size: the default relax sequence ends in an NPT
  segment of 2000 steps and `relax-npt.mdp` runs `dt = 0.001`, so 2 ps per
  stage; at the default `nstages: 6` over ~10 iterations that is roughly
  **120 ps of above-`Tg` constant-pressure time across an entire cure**,
  against a 39.2 ns production ladder. Density at cure termination is still
  climbing, bpa +1.29 +/- 0.62 and bpf +1.48 +/- 0.28 kg/m3 per iteration
  over the second half of cure (2.1 and 5.3 sigma, n = 4 per chemistry) --
  though a rising density is *partly* expected, since crosslinking genuinely
  densifies, and that slope cannot separate real densification from
  incomplete relaxation. Reactive-species mobility over cure: bpa 18.1 ->
  3.7 A with crossing fraction 97 % -> 13 %, bpf 17.4 -> 4.4 A and
  97 % -> 22 %. All measured by the study session.

- **Published gating designs we have not evaluated.** Audited by the study
  session 2026-09-05. These gate *bond acceptance*, not density, so R16's
  refutation of the density gate leaves them untouched.

  - **Rejection on relaxation failure -- Moore et al. 2021**
    (*Macromolecules* 54:6275), on di(cyanate ester)s, which is *our*
    chemistry. Verbatim: "If a failure occurred during relaxation following
    bond formation, the bond was removed and the crosslink cycle was
    repeated with different bonds selected for crosslinking." That is bond
    *acceptance* conditioned on a measured relaxation outcome. The same
    paper caps bonds at **1 % of total per iteration** and terminates on
    **stall** -- no new bonds in 20 iterations -- rather than on a step
    count. Both are directly comparable to knobs we have:
    `max_conversion_per_iteration` defaults to 1.0, i.e. no cap at all, and
    we terminate on conversion or `max_iterations`.
  - **Arrhenius kinetics gating -- Schichtel & Chattopadhyay 2020**
    (*Comput. Mater. Sci.*): bonding probability from cure temperature,
    cutoff distance and an activation energy. Their sec 2.4 *demonstrates*
    the failure mode rather than merely flagging it -- a density trajectory
    they call "not physical", ending in a metastable configuration, from
    "accelerated reaction kinetics caused by the usage of higher
    probabilities."

  Rejection catches a bond that will not relax at all; Arrhenius sets the
  rate. Neither sees gradual under-relaxation, but neither costs what the
  density gate would have, and both act where R16 says the leverage is not
  -- on which bonds form, rather than on how long the box is given
  afterwards.

  **A calibration number worth having before anyone panics about our
  deficit.** Moore uses a fixed 40 ps NPT relaxation window at 800 K -- the
  same number Varshney chose in 2008, with the same stated rationale, to let
  unreacted species diffuse between reactions -- and reports a **4.77 %**
  density deficit in BADCy after 100 ns of post-cure equilibration at 293 K.
  Ours is 2.25 %. So the fixed-window architecture is universal from 2008 to
  a 2021 GPU workflow, the resulting deficit does not equilibrate away, and
  htpolynet is not an outlier. That makes this a field-wide design question
  rather than a bug, which raises the value of getting it right and lowers
  the urgency.

  Nobody here has read either paper; the quotations above are the study
  session's, and both citations should be checked against the originals
  before being repeated anywhere public. The same applies to Varshney 2008,
  which the v2.7.0 docs now name as the origin of the 40 ps relaxation
  window and its mobility criterion -- that attribution is the study
  session's too, and it is the one of the three that has already been
  published in our documentation. Verify it first if any of them.

## Configuration

- **Move config handling onto ycleptic.** Approved by Cameron 2026-09-08.
  `core/configuration.py` (85 lines) reads YAML/JSON and `.get()`s eleven
  known top-level sections into attributes with **no validation at all**: an
  unrecognized section or a misspelled key is silently discarded. Defaults
  live in three separate Python dicts applied by three hand-rolled mechanisms
  at three different depths -- `Runtime.runtime_defaults` (shallow, one level,
  `core/runtime.py:177`), `CureController.curedict_defaults` (two levels,
  `cure/curecontroller.py:299-304`), and `Runtime.default_edict` (ad hoc). And
  937 lines of hand-maintained RST document those defaults, which is why six
  of them had drifted from the code by the time anyone checked.

  `src/htpolynet/schema/base.yaml` now holds the schema for the sections
  ycleptic 2.3.0 can express. It is **inert**: nothing reads it, and
  `pyproject.toml` does not declare a ycleptic dependency. All seven shipped
  example configs validate against it, and a `desired_converson` typo that
  htpolynet accepts silently today is rejected with the list of valid keys.

  Staging, agreed with the ycleptic session:

  1. **ycleptic to conda-forge** -- gates shipping and nothing else, and is
     now the ONLY thing gating it.
     `api.anaconda.org/package/conda-forge/ycleptic` was 404 while htpolynet
     is 200, and a conda package's run-deps must exist in the channel, so
     declaring the dependency would make htpolynet's feedstock unbuildable and
     trip `scripts/release.sh`'s `check-conda-sync.py --strict` preflight.
     Submitted 2026-09-08 with Cameron's direct approval as
     conda-forge/staged-recipes#34763, "Adding ycleptic"; still open, not
     merged, and `api.anaconda.org/package/conda-forge/ycleptic` is still 404.

     **That PR pins 2.3.0, and 2.4.0 is now out on PyPI.** As written, merging
     it creates a feedstock one version behind, and 2.4.0 then arrives only
     through an autotick-bot PR -- the exact cycle whose stalling is the
     central OBSERVED failure in Cameron's conda-forge skill. Updating the
     staged recipe to 2.4.0 before it merges would skip that cycle entirely.
     It is ycleptic's PR, not ours, and this is recorded as an observation
     rather than an action.
  2. **Three ycleptic grammar additions**: `value_attributes:` + `key_text:`
     for free-key mappings, `value_type:` for scalar-valued ones, and
     `list_defaults: replace|append` per attribute. **Released as ycleptic
     2.4.0 on PyPI**, and base.yaml uses all three.
  3. Port the flat sections and generate their reference docs -- done in
     base.yaml, awaiting 1 to go live.

     **The first release carrying ycleptic (2.8.0) was cut ahead of
     conda-forge, deliberately.** An earlier version of this entry said a bot
     PR without `ycleptic` in the recipe would be green and would ship a conda
     package that fails at import, and so forbade `--skip-conda-check`. That
     premise was wrong, and checked 2026-09-12: the feedstock's test runs
     `htpolynet --help`, and importing `htpolynet.cli` loads ycleptic. So the
     bot's PR for a ycleptic-carrying release **fails its own test**, sits red,
     and cannot auto-merge. No broken conda package ships; conda-forge users
     simply stay on the last good version.

     What that leaves to do once ycleptic is on the channel:

     a. Edit the recipe on the autotick bot's open PR for the current version
        (push to the bot's fork, not to `origin` -- see the conda-forge skill)
        adding `ycleptic >=2.4.1` to `run:`, and `>=3.5` to the bare
        `matplotlib-base`. Its test then passes and automerge takes it.
     b. If the bot has not opened one, open that recipe PR by hand.

     Two cautions that still stand. The protection is the recipe's
     `htpolynet --help` test; if that test is ever removed, this reasoning
     stops holding and `--skip-conda-check` becomes dangerous again. And a red
     bot PR must not be merged by hand -- red is the thing keeping the broken
     package off the channel.
  4. Port `constituents` and `reactions`. **`constituents` is done**: it is a
     `value_attributes` free-key mapping, so a typo inside a molecule record
     is now rejected naming the molecule (`Attribute 'smilse' invalid ...
     under 'constituents[STY]'`), `count` defaults per value, a required
     `conformers.count` is enforced, and `conformers.generator.name` is
     constrained to obabel/gromacs. Every list with a non-empty default
     carries `list_defaults: replace`, verified to use a user's ladder
     verbatim rather than appending it to the default one.

     **`reactions` and `postcure_repair` are now ported too.** ycleptic 2.4.0
     added `value_attributes` on a *list*, which describes multi-key records
     and closed the `lwalk` gap, plus `value_type:` for scalar-valued free-key
     mappings. Every config section is now validated except
     `gromacs.mdrun_options`, which is deliberately open. Errors name the
     offending item -- `under 'reactions[2]'`, `of 'postcure_repair[0]'`.
     Today an unknown reaction key is discarded by `Reaction.__init__` with a
     debug-level message, so this replaces a silent failure with a loud one.

     Superseded, kept for the record -- these were the two blockers:

     - ycleptic describes list items with the tagged-task idiom -- `lwalk`
       takes `itemname = list(item.keys())[0]` -- and htpolynet's reaction
       records are multi-key. With no schema node for the item there is
       nowhere to hang `value_attributes` for `reactions[].atoms`, which
       would otherwise qualify. The ycleptic session has deliberately not
       changed `lwalk`; it is a real feature request.
     - `value_attributes` requires every value to be a **mapping**, which
       ruled out scalar-valued free-key mappings. **Closed**: `value_type:`
       was bundled into ycleptic 2.4.0, and `constituents[].reactive_atoms`
       and `constituents[].rename_atoms` now declare `value_type: str` and
       are validated. `reactions[].reactants` (`{1: BPA, 2: HIE}`) is the
       same shape and would work too, but it lives inside a `reactions` list
       item, so the `lwalk` limitation above blocks it regardless.

- **`required: true` on a list attribute is silently ineffective in ycleptic
  2.4.0.** An absent required list is filled with `[]` and no error is raised
  (`walkers.py:389-397`); scalars raise, and dicts raise via their required
  children. Verified with a minimal spec, and it is the exact class of
  do-nothing declaration `speccheck` exists to catch, which does not catch
  this one. In base.yaml it means `reactions[].bonds`, `bonds[].atoms` and the
  repair spec's `ring_carbon_atoms`, `ring_nitrogen_atoms` and
  `reactive_oxygen_atoms` are declared required and are not enforced: omitting
  one validates, then fails deeper in htpolynet. The declarations are kept
  because they state the intent correctly and will start working when ycleptic
  fixes it. Reported to the ycleptic session; not worked around here.

  **The version floor is a trap; read this before pinning anything.** The
  staged-recipes submission is necessarily at ycleptic **2.3.0**, because
  staged-recipes needs a tarball that already exists. Step 3's base.yaml is
  written against 2.3.0 and works there. But the free-key node step 4 needs
  will land in **2.4.0**, and 2.4.0 reaches conda-forge only via an
  autotick-bot PR *after* that release -- and the central OBSERVED failure in
  Cameron's conda-forge skill is exactly those bot PRs sitting red for months
  while conda-forge serves the stale version. So: pin `ycleptic>=2.3.0` when
  step 3 goes live, and **do not raise the floor to 2.4.0 until the anaconda
  API confirms 2.4.0 is actually on the channel** -- not when it is released,
  not when the bot PR opens. Raising it early recreates precisely the
  unbuildable-feedstock problem step 1 exists to avoid. The ycleptic session
  will report when the API says so.

  Two hazards are recorded in the base.yaml header rather than here, because
  that is where someone editing it will look: ycleptic's `dwalk` rejects any
  key not declared, so user-keyed sections must currently be bare `dict` and
  are unvalidated; and a bare list's default *concatenates* rather than
  replaces, so a base default of `[min, nvt, npt]` plus a user's `[nvt]`
  yields `['min','nvt','npt','nvt']` silently. Every list in base.yaml
  therefore declares no default, which gives correct replace semantics, with
  the real default left in Python and marked.

- **Documented `CURE.drag`/`CURE.relax` defaults do not match the code, and
  several keys are undocumented.** The tables give `drag.increment` 0.08 and
  `drag.limit` 0.3 where `curedict_defaults` has 0.0 for both -- sentinels
  that disable dragging via the `:316` guard, so the documented defaults would
  turn on a stage ladder the code turns off. (`relax.increment` was the same
  kind of mismatch and is now fixed in the code's favor: it defaults to the
  documented 0.08. `relax.nstages: 6` is still declared and never read.)
  Undocumented entirely:
  `drag.trigger_distance`, `drag.kb`, `drag.nstages`, `drag.cutoff_pad`,
  `relax.nstages`, `relax.cutoff_pad`, the whole `CURE.output` block,
  `densification.scale`, `densification.aspect_ratio` and
  `densification.initial_boxsize`. Generating the reference from base.yaml is
  what makes this class of error impossible; until then the tables need a pass.

- **`resolve_type_discrepancies` is accepted in two places and documented in
  one.** `core/runtime.py:326` falls back from
  `cfg.gaff.get('resolve_type_discrepancies',[])` to
  `cfg.resolve_type_discrepancies`, so both a top-level key and a
  `GAFF`-scoped one work, with GAFF winning. The
  docs describe only the GAFF-scoped form. base.yaml declares both, marking the
  top-level spelling deprecated; pick one before the schema goes live.

- **A parameterized molecule built from a symmetry-sibling reactant would
  hand antechamber an incomplete bond table.**  Found while fixing the
  bondless-MOL2 writes reported by the 2026-09-13 example sweep. Symmetry
  siblings are loaded from the parent's top+gro and carry no `mol2_bonds`.
  When topologies merge (`Topology.merge` via `_myconcat(... 'mol2_bonds')`),
  a reactant without that table contributes nothing to the product's, so the
  product's `{name}.mol2`, which `Molecule.parameterize` gives to antechamber
  for atom typing, would silently lack that reactant's internal bonds.
  `TopoCoord.write_mol2` now falls back to topology bonds only when the table
  is absent entirely, which does not cover this partial case. Not exercised by
  any shipped example: the fresh-cache 2026-09-13 sweep had zero bondless
  writes outside diagnostic `-prebonding` files, and every sibling product
  copied its parent's topology rather than being parameterized. The fix, if
  a case ever appears, is for `merge` to synthesize `mol2_bonds` from `bonds`
  for any operand that lacks the table, rather than concatenating a partial
  table. Bond orders would then need the same single-bond assumption.

  **It has appeared (2026-09-17).** Example 2's stereoisomer copies `GMAS-*`
  have 78 bonds in their topology but only 24 in their `.mol2`, in the
  2026-09-14 sweep's files too. It shows up only when a cache is *reused*:
  `bond_source_check` then logs "Gromacs/Mol2 bond inconsistency detected"
  once per copy. Sweeps always start from an empty cache, so they never see
  it. `generate_stereoisomers` deep-copies the parent's TopoCoord and writes
  only a `.gro`; where the 24-bond `.mol2` comes from is not yet traced. What
  a reused cache does with the wrong table (ring detection, system `.mol2`
  output) is also unknown. No cyanate-ester molecule has stereoisomer copies.

- **Two gaps left after the percolation check moved to cure sites
  (2026-09-15).**  First, a system that crosslinks through a *multi-residue*
  molecule -- example 2's GMA (two methacrylate HIE units), example 5's
  three-armed THT -- is reported as "not assessed". A real assessment needs a
  molecule-level count: sum cure sites per molecule and ask how many molecules
  reacted at all of them. Second, `CureController._completion_bias_counts`
  still ranks candidates by `residue_reaction_counts` over *all* atoms, so a
  bond made while building a molecule (HIE's C4) counts as progress toward
  completion. That only matters with `completion_bias` on, which no shipped
  example uses. It was left alone because fixing it changes which bonds a
  biased cure forms, not just a diagnostic; `cure_site_mask` is the tool when
  that is wanted.

## Simulation defaults

- **The density-convergence tolerance may be tight for example 6.**  In the
  v2.10.0 sweep, example 6 failed the gate for the first time in four sweeps:
  "sem 1.070 exceeds tolerance 1.000" after 6 extensions, with the last window
  at 1112.0 kg/m^3 against 1110.5 and 1111.6 in the two sweeps that converged
  with no extension at all.  Monomer charges were byte-identical to the
  previous sweep, so this is run-to-run noise rather than a change -- but the
  build went on into the cure from a box the gate had just called unsettled,
  and said so only in a warning.  Worth deciding whether the tolerance should
  scale with system size or window length, and whether failing it should be
  more than a warning.

- **The halogen constraint failure is fixed but not explained.** v2.7.0
  switched the per-iteration equilibration to `constraints = h-bonds` with
  `lincs_order = 8`, which took a fluorinated bisphenol from 0 of 7 builds to
  4 of 4, and the eight-bridge series from 7 failures to 32 of 32. What is
  still not understood is the *timing*. The constraint artifact is present
  from the very first NPT step -- a melt of pristine, uncrosslinked BAF
  monomers at production density shows step-0 pressure of -1.34e5 bar, with
  no cure involved at all -- yet the old builds survived densification,
  precure and 4-6 cure iterations before dying. So it was a predisposing
  cause plus a threshold in the network's tolerance for it, not a single
  trigger, and nobody knows what sets the threshold. The failing structures
  were not preserved, so the LINCS atom indices were never mapped to
  residues; anyone chasing this should trap and keep them first. The
  diagnostic that identified it is worth reusing: halving `dt` made the
  artifact ~4x *worse*, and 1/dt^2 scaling is the signature of a one-shot
  constraint start-up projection rather than a physical clash.

- **`mdp_to_dict` crashes on any mdp line containing two `=` signs.** In
  `external/gromacs.py:239` it does `k,v = l.split('=')`, which raises
  `ValueError: too many values to unpack` rather than reporting anything
  useful. No packaged template trips it today -- checked across all eleven --
  but any comment containing an equals sign does, including the obvious one
  someone would write next to a timestep. It also does not strip `;`
  comments, so the comment text becomes part of the value and is re-emitted
  attached to it; GROMACS tolerates that, but `mdp_get` on a commented
  numeric key would not survive a `float()`. Split on the first `=` only,
  and strip comments.

- **`postsim` inherits `npt.mdp` wholesale, so cure defaults govern the
  production measurement.** `analysis/postsim.py:81` checks out `npt` and
  renames it; `build_mdp()` overrides `ref_t`, `ref_p`, `nsteps`, velocity
  generation and `tcoupl` (to `v-rescale`) but never touches `dt`,
  `constraints`, `lincs_*` or `pcoupl`. So whatever is chosen above also sets
  the timestep and constraint scheme of every postsim run, and the barostat
  there is still **Berendsen**, which does not sample a correct NPT ensemble
  and is deprecated in modern GROMACS. Changing the barostat is a physics
  change to published numbers; changing it silently would be worse than
  leaving it, so it wants a release note and probably a config knob.

## Usability

- **The atom-serial column in `final.gro` wraps modulo 10000, and nothing
  says so.** That is the GROMACS `.gro` format, not an htpolynet bug -- atom
  10000 prints as `0` -- but every htpolynet system large enough to matter
  crosses it, and the natural analysis script joins `final.gro` to `final.top`
  on that column. In one 12,600-atom build 2,601 atoms (21%) have a `.gro`
  serial that differs from their `.top` index, so such a script is silently
  wrong for a fifth of the box with no error anywhere. Worth a sentence in the
  analysis docs saying to join on line order, not on the serial. Same family
  as the existing trap that `min_clearance_nm` is a substring of
  `blind_min_clearance_nm` while PyYAML emits keys alphabetically, so a naive
  grep of `repair-summary.yaml` picks the wrong one.

- **`gen-slurm-script` doesn't stage to scratch.** The emitted script
  runs in the submit directory. A cure run does heavy small-file I/O
  every iteration, so on a cluster whose home and group storage are NFS
  that is the wrong place. The cluster-correct pattern is to run in
  node-local or parallel scratch and copy results back, with a `trap` so
  partial results survive a timeout. It also cannot emit a job array,
  which is the right shape for sweeping a set of configs.
- **`-restart` is documented as "EXPERIMENTAL: broken at the cure
  stage".** A build that dies mid-cure currently has to start over — the
  worst possible time to lose work, since cure is the longest stage.
- **Make `input-check` a real config linter.** It currently reports only
  the initial atom count. Everything a new config gets wrong is checkable
  cheaply and statically, before the user spends hours discovering it:
  (a) that each `symmetry_equivalent_atoms` group really is topologically
  equivalent — RDKit canonical ranks with `breakTies=False` settle it in
  about ten lines, and a wrong group silently generates reaction templates
  the cure stage will never match; (b) that the A2+B3 site counts actually
  balance at the given monomer `count`s, and that `desired_conversion` is
  reachable given them; (c) that no `reactive_atoms` name collides with the
  bare element names `_reset_names_to_element` assigns to every unmapped
  atom, since the repair stage looks atoms up by name; (d) that the
  requested `charge_method` matches the provenance record of any cached
  parameterization of the same name -- the build itself now rejects a
  mismatch, but reporting it in the pre-flight is cheaper than discovering
  mid-run that half the molecules need re-parameterizing. Building the
  cyanate-ester bridge series meant doing (a) and (b) by hand in a throwaway
  RDKit script, which is exactly the work a user should not have to
  reinvent.
- **Generated topologies cannot be compared byte-wise, because ParmEd
  stamps them.** Every `.top` htpolynet writes opens with a ParmEd header
  recording the invoking user, the host, and the date:

      ;   File TAZ.top  was generated
      ;   By user: cfa (1000)
      ;   On host: panacea.chemeng.drexel.edu
      ;   At date: Fri. May  5 15:55:38 2026

  So two physically identical parameterizations never hash the same, and the
  obvious check -- "did these two builds produce the same parameters for the
  shared monomers?" -- returns a false difference. The calibration study hit
  this on its bridge series and nearly reported the series confounded on the
  strength of an md5; the fix on their side was to strip comment lines before
  comparing, which every user will have to reinvent. Note the header is
  itself provenance, but the wrong kind: it records *when and where* rather
  than *what directives*, which is what the `parm` record now covers, and it
  actively prevents the comparison you would want. The `.itp` files carry no
  such header and the `parm` records are JSON with sorted keys, so both
  already compare cleanly. Worth either normalizing the `.top` header away or
  shipping a comparison helper; this belongs with the manifest entry below,
  since both are about being able to answer "is this the same build?".

- **`-lib` is read-only, and nothing says so.** `pfs.checkout()` and
  `pfs.exists()` consult the `-lib` user library first, then the user cache,
  then the system library -- but `pfs.checkin()` writes unconditionally to the
  user cache, and `UserLibrary` has no `checkin()` method at all. So pointing
  a run at `-lib somewhere` redirects where it *reads* parameterizations from
  while leaving where it *writes* them untouched. That asymmetry is invisible
  from the flag, from `--help`, and from the docs, and it is a natural thing
  to get wrong: a user who passes `-lib` for provenance reasonably concludes
  their products are being contained there. Either `checkin()` should prefer
  the user library when one is configured, or the flag and docs should say
  plainly that it governs lookup only. Verified empirically 2026-08-23: with
  `userlibrary` set, `pfs.checkin()` still landed the file in the user cache.

  What makes this more than a naming problem is that the writes are silent
  and cumulative. Nothing logs them, nothing warns, and a run that is
  quietly depositing molecules into a shared library looks identical to one
  that is not -- the only way to discover it is to audit mtimes, which
  nobody does unprompted. The calibration study ran three sweeps believing
  `-lib` contained it, and would have added roughly seventy entries to
  `~/.htpolynet` across its remaining bridges; it was caught by reading
  `checkin()`, not by anything the tool said. Whatever the fix, a run should
  be able to say where its parameterizations went -- which is the same
  requirement as the build manifest below, approached from the write side
  rather than the read side. The two may well be one item: *where did this
  come from, and what produced it*.

- **No way to run without writing to the user library.** `pfs.checkin()`
  declines to *replace* an entry unless `--force-checkin` is given, but it
  always *writes* one the library does not yet hold, so any run adds every
  molecule name it produces to `~/.htpolynet`. The flag's help text said
  "force check-in of generated parameter files to the system library", which
  reads as though check-in happens only with the flag, and it named the wrong
  library; that wording is fixed, but the missing capability is real. The
  calibration study wanted exactly this: a way to develop against a config
  without its intermediate products accumulating in a library shared with
  other work. Today the only lever is pointing `HTPOLYNET_CACHE` elsewhere,
  which is a blunt instrument because it also hides the entries you *do* want
  to reuse. A `--no-checkin` flag threaded through
  `Runtime._checkin_parameterization()` would cover it in a few lines.

- **Two different defaults for `charge_method`.** `AMBERTOOLS_DEFAULTS` in
  `external/ambertools.py` says `bcc`, which is what a direct
  `GAFFParameterize()` call with no directives gets; `Runtime.runtime_defaults`
  says `gas`, which is what every actual build gets, because
  `_apply_runtime_defaults()` fills it in before AmberTools is ever reached.
  Both were already there and the provenance record is consistent either way
  -- whichever default applies is the one recorded -- so this is a
  readability trap rather than a live bug. It is still worth collapsing to
  one value, and the answer is probably `gas`, since changing what builds
  default to would silently change everyone's charges, which is the exact
  class of harm the record was added to prevent.

- **Nothing durable records that a build reused parameterizations of
  unverified provenance.** The parameterization stage now warns per molecule
  and again in a block at the end of the stage, but both live only in the
  log, and a log is the first thing discarded. Someone reading a result six
  months later -- or a reviewer asking what a published network was actually
  parameterized with -- has no artifact to check. htpolynet writes no build
  manifest today; `profile.json` holds timings only, and the diagnostic log
  is the log. A small `build-manifest.json` in the project directory listing
  each molecule, its origin (`newly parameterized` / `previously
  parameterized`), and its provenance record where one exists would carry
  that, and would be worth more than this one flag: it is also the natural
  home for the config hash, the htpolynet and AmberTools versions, and the
  seed once seeds exist. Raised by the calibration study, which caught the
  original cache bug from a wall-clock anomaly rather than from any log line.
  See also the `-lib` entry above: "a run should be able to say where its
  parameterizations went" and "a result should be able to say what produced
  it" are the same requirement from two directions, and a manifest that
  recorded the check-in destination would answer both.

- **A cached parameterization is not checked against the input structure
  it was built from.** The provenance record added alongside the cache now
  covers `charge_method`, `net_charge` and `atom_type` -- everything in the
  AmberTools invocation -- but not the structure antechamber consumed. So
  editing `lib/molecules/inputs/TAZ.mol2` and re-running still reuses the
  parameterization of the *old* geometry, silently, exactly as the charge
  method used to. Hashing the input would close it, and for a monomer that
  is easy: `Molecule.parameterize()` has the input file in the working
  directory at the moment it runs. Two things stopped it going in with the
  rest: (a) a molecule built by a reaction has no stable input to hash --
  `generate()` writes its mol2 from the merged reactant TopoCoord, whose
  coordinates vary run to run, so hashing it would make every generated
  molecule a permanent cache miss; the hash would have to be recorded only
  for `origin == 'unparameterized'` monomers and compared only when both
  sides carry one. (b) `Runtime._cached_parameterization_mismatch()` runs
  before `generate()` checks the input structure out, and `pfs` has no way
  to resolve a library file to an absolute path without copying it into the
  working directory -- `checkout()` always copies. A `pfs.locate(filename)`
  returning the resolved source path across user library, user cache and
  system library is the missing piece, and is worth having on its own.

- **No seed control anywhere, so a build cannot be reproduced.** Three
  independent sources of randomness are all unseeded: `random.sample` for
  conformer selection (`core/runtime.py`), `np.random.random()` for the
  per-bond probability test (`cure/curecontroller.py`), and every `.mdp`
  sets `gen-vel = yes` without `gen-seed`, so Gromacs picks a pseudo-random
  one per run. Two runs of one config therefore diverge — example 6 gave 57
  vs 60 incomplete triazines on two machines at the same 0.90 conversion.
  That is convenient in one direction (independent replicas of a
  quenched-disorder ensemble come free by re-running) but it means a build
  reported in a paper cannot be reproduced exactly, and a failure seen once
  may not reappear. The fix is a top-level `seed:` that feeds all three:
  seed `random` and `np.random` at startup and write `gen-seed` into every
  generated mdp, with replicas then requested by varying it rather than by
  relying on entropy.
- **Drop the pre-3.5 matplotlib fallback** in `analysis/plot.py`'s
  `_get_cmap()` once `matplotlib>=3.6` is a safe floor; the
  `matplotlib.colormaps` registry is then always present.

- **The 18 postsim TODO placeholders are the last unfinished tutorial work.**
  Every one of the seven `postsim.rst` pages carries `**TODO:** insert the
  rho_v_ns.png plots`, `**TODO:** insert tg.png` and `**TODO:** report the
  fractional free volume`, and `6-cyanate-ester/monomers.rst` still wants a
  `pics/TAZ.png`. These are unfinished rather than stale: nothing in them is
  wrong, the numbers and figures were simply never produced. Filling them
  needs real compute -- `htpolynet postsim` runs, which the 2026-09-08 example
  sweep did not do, since it ran `htpolynet run` only.

- **Tutorial timings are hardware-dependent and nothing said so until now.**
  The stage tables in examples 4, 5 and 6 were rewritten on 2026-09-09 from
  measured 24-core runs and now say "on 24 cores". Before that, example 5's
  table quoted ~12 hours from a 16-core run of an earlier configuration of
  that example, against 3h45m measured now -- and its iteration count was 15
  against 9 measured, which is a real behavior change rather than hardware.
  Whatever replaces this should either state the machine or stop quoting
  absolute times.

- **`htpolynet analyze` has no MSD stage, and one would have to be written
  carefully.** A user asked for an MSD workflow on a DGEBA/PACM network
  (2026-09-12), and the only thing to point them at was `gmx msd` directly.
  The trap for whoever adds one: `gmx msd` defaults to `-rmpbc`, which tries to
  make each molecule whole every frame, and for a periodic network that is
  impossible. It does not fail -- it prints "inconsistent shifts" and silently
  inflates the MSD, by about 2.4x at 100 ps on example 3 (0.379 against 0.154
  nm^2). The correct recipe is `trjconv -pbc nojump` and then
  `gmx msd -normpbc`. Avoid `-nopbc`, which segfaulted `gmx msd` in GROMACS
  2025.4. `-mol` is meaningless, since the whole network is one molecule. The
  documented recipe is on the "Analyzing Trajectories of Periodic Networks"
  page.

  Relatedly, the existing `freevolume` stage passes no `-normpbc`, so it also
  prints the warnings. For free volume they are harmless -- 0.211 +/- 0.007
  with default settings against 0.210 +/- 0.006 with `-normpbc` -- but they are
  noise a user will reasonably worry about. Adding `-normpbc` to its defaults
  would silence them; it is not done yet because `Analyze`'s option dict
  expects key/value pairs, and a bare boolean flag needs checking there.

## Example depot

- **Example 6's postcure anneal peaks too close to *T*:sub:`g` to relax the
  network, and the fix is not simply "run it longer".** Measured on four
  independent BPA builds: the postcure NPT plateau gives 1.1712 ± 0.0029
  g/cm³ where the same systems melted and slowly re-cooled give 1.1983 ±
  0.0037 g/cm³ -- 2.31% apart, with the plateau 2.6-2.8% below experiment
  and the re-cooled value within 1%. The plateau is an under-relaxed
  structure.

  The mechanism, and the reason the obvious fix is the wrong one: the
  example anneals at a 500 K peak, and the measured *T*:sub:`g` for this
  system is 487.9 K. That is 12.1 K above the glass transition, for 80 ps.
  The Tg ladder that produced the relaxed value spent 15,785 ps above
  *T*:sub:`g`, peaking 112 K above it -- 197x the time, in the regime where
  a crosslinked network actually moves. Lengthening an anneal that sits 12 K
  above *T*:sub:`g` extends a process that barely moves anything, so the
  lever is the peak temperature, not the duration.

  **But raising the peak is unlikely to close the gap entirely.** Two
  structures of the same network were held for 5 ns at a 480 K setpoint
  that thermostatted to 477.1 K -- about 11 K *below* the measured
  *T*:sub:`g`, and far longer than the 80 ps the example spends near it.
  Relaxation there is fast at first and then stops: about 87% of an
  initial 34.8 kg/m³ density gap -- the gap measured over the first 200 ps
  -- closes within 3 ns, and the remaining 4.25 kg/m³ persists at
  7.2 sigma, with a last-half trend of +0.030 kg/m³/ns, i.e. not decaying.
  (State the baseline and the window whenever this number is quoted: the
  same data give 83.5-87.8% across the defensible choices of each, so a
  bare percentage is not a fact.)  So heat and time buy most of the
  relaxation quickly and then buy nothing, and a residual difference
  survives that a longer anneal at this scale does not appear able to
  remove.  Whatever replaces the current anneal should
  therefore be argued as a large improvement, not as a fix; an entry
  promising equilibration would oversell what was measured.

  Note what that hold also says about *T*:sub:`g` itself: **relaxation does
  not switch off at the glass transition**, it slows continuously through
  it.  A transition 20-34 K wide by experiment means 11 K below
  *T*:sub:`g` is inside the transition, not deep in the glass, so fast
  early relaxation there is expected rather than anomalous.  This does not
  soften the argument above -- the 197x figure is unchanged and a peak well
  clear of the transition is still the lever -- but it does mean *T*:sub:`g`
  is a soft kinetic boundary for this material, which any protocol
  expressed relative to it has to accommodate.

  **Untested.** Nobody has yet built with a hotter anneal and compared it
  against the melt-and-recool value, which is the experiment that would
  settle it; the argument above is mechanism, not measurement. That is why
  the example is documented rather than changed. Raising the peak also costs
  wall-clock in a tutorial that should stay small, so the test needs to
  report both whether it works and what it costs.

  The general problem is worth more than the instance: **any protocol whose
  relaxation depends on being above *T*:sub:`g` is fragile when *T*:sub:`g`
  is unknown at configuration time**, which it always is -- you cannot know
  it before building the thing. A fixed absolute peak temperature is a guess
  that happens to be right or wrong per chemistry, silently. Something
  expressed relative to an estimated or measured *T*:sub:`g`, or a postcure
  that measures *T*:sub:`g` and then anneals accordingly, would be robust
  where an absolute number is not. That is a real design change, not a
  config edit.

  Raised by the calibration study, which had already published the
  pessimistic density before catching it -- the failure mode is reporting a
  protocol artifact as a force-field result.

- **Example 5 is still not laptop-scale.** The `68d81fb` retune cut the
  monomer pool to make HTPB/IPDI "fit in a reasonable wall time on a
  laptop", but a validation run on 16 CPU cores took **10h17m** — about
  ten times any other example (the next longest, PACM/DGEBA, was 2h15m).
  The cost is structural: densification starts at 10 kg/m³ in a 27 nm box
  and needs all 50 NPT repeats to compress geometrically (a steady 2.9%
  box-side reduction per repeat) up to ~860 kg/m³, and the cure stage then
  does drag/relax on long chains. Worth revisiting whether the low
  starting density is really necessary, or whether the example can start
  denser with a shorter densification.


- **More cyanate-ester variants.** Example 6 covers bisphenol-A
  dicyanate. The same reaction and repair machinery carries over to other
  bisphenol bridges with only a SMILES swap in `constituents` — bisphenol
  F (methylene), bisphenol E (methylethylidene), hexafluorobisphenol A,
  thioether, sulfone, dicyclopentadiene. A homologous series would
  exercise the repair stage across chemistries and give the tutorials a
  structure-property story.
