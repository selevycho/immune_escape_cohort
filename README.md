# Immune Escape Convergence — full cohort

The code behind forty synthetic tumour-normal pairs, and the analysis
that measures what each stage of a neoantigen pipeline recovers from
them.

Each pair puts somatic mutations from a TCGA breast or ovarian tumour
into panel sequencing from a healthy 1000 Genomes individual. Neither
half is simulated; the pairing is. Because every mutation was placed at a
known position with a known allele fraction, every tool downstream can be
scored against what is actually in the file rather than against another
tool's opinion.

**This repository is not meant to be run.** It depends on a cluster, a
workspace layout and thirty gigabytes of reference data specific to where
it was built. It is here so the analysis can be read and the numbers
traced back to the code that produced them.

For a version that does run — five samples, eight stages, fetching its
own data — see
[immune-escape-demo](https://github.com/selevycho/immune-escape-demo).

---

## The data

The forty pairs are published at **[DOI to be inserted]**: eighty BAMs
with indexes, a truth table listing all 1 611 placed mutations, and the
manifest pairing each backbone with its donor. 26 GB.

---

## What is here

```
scripts/simulation/     the pipeline itself, one file per stage
scripts/01_cohort_tcga/ building the cohort: TCGA download, liftover,
                        panel construction, sample selection
scripts/02_verify/      every check and sweep reported in the results
scripts/03_figures/     the figures, and the plotting library they share
results/                the tables those scripts produce
panel/                  350 genes, 23.81 Mb, 5 373 intervals
```

Each script carries its reasoning at the top — what it does, what was
tried first, and why the current version is shaped the way it is. Several
record mistakes that cost days to find, which is the part most worth
reading.

---

## What the numbers say

```
1 611 mutations placed, 1 572 present in the reads

Mutect2         74.5% of substitutions, 82.7% of indels
                no false calls across the cohort
OptiType        forty distinct genotypes, 78.2% correct
                against laboratory-typed controls
mhcflurry       1 830 strong binders from 49 863 peptides
NetMHCpan       agrees on sample ranking at r = 0.975
LOHHLA          fifteen simulated losses detected,
                none reported that was not made
```

Three findings a consensus-based benchmark could not have produced:

**The missing quarter is largely a filter setting.** Loosening
Mutect2's filters recovers 222 more true mutations without a single false
call.

**The MHC is invisible.** None of the 37 substitutions placed inside it
were recovered — not filtered out, never considered — and on files with
no mutations at all the caller reports variants there at eight times the
density of anywhere else.

**HLA loss detection runs out of positions, not reads.** 338 of 376
locus attempts never reached the statistical test, because the two
alleles of a locus rarely differ at enough covered positions on panel
data.

---

## Reading order

Start with `scripts/simulation/step2_inject.sh` — it explains how a
mutation is placed into real reads and what realignment does to it.

Then `scripts/02_verify/check_mutect2.py` for how recall is scored, and
why it is scored against verified mutations rather than against the truth
set as written.

`scripts/02_verify/check_lohhla.py` is the most heavily commented, for
the reason that it took the longest to get right.

---

## Licence

MIT for the code. The tools it calls are under their own terms and none
is redistributed here: LOHHLA under the Francis Crick Institute's
academic terms, netMHCpan under DTU's, IMGT/HLA under CC BY-NoDerivs.
