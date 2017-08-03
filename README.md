## Regulatory Element Locus Intersection (RELI) Analysis

A one-paragraph description of what the tool does.

## Installation

RELI requires a C++11 compiler (_e.g._ GNU CC 4.7 or higher) and `libgsl` and
`libgslcblas` from the [GNU Scientific Library][gsl].

A GNU-style `Makefile` is provided. You can build the RELI binary with just

    make

or, optionally, download sample data and run a test analysis with

    make test

### Manually downloading sample data

If you don't have `curl` available, you can manually download and extract the
sample datasets from

    https://tf.cchmc.org/external/RELI/data.tar.bz2

such that the decompressed data is inside a `data` subdirectory, within the
`RELI_public` repository you cloned above. A `.zip`-format archive is also
provided, in case for some reason you don't have `bzip2` available.

You can run the sample analysis by changing into the `example` directory and
running `example_run.sh` in a terminal like so:

    user@[/path/to/repo]$ cd example
    user@[/path/to/repo]$ ./example_run.sh

## Command-line options

_Required options are in **bold text**_

| Option                | Explanation
|-----------------------|------------------------------------------
| **`-snp FILE`**       | Phenotype snp file in 4 column bed format
| `-ld FILE`            | (optional) Phenotype linkage disequilibrium structure for snps, default: no ld file
| **`-index FILE`**     | ChIP-seq index file
| **`-data DIR`**       | Specify directory where ChIP-seq data are stored
| **`-target STRING`**  | Target label of ChIP-seq experiment to be tested from index file
| **`-build FILE`**     | Genome build file
| **`-null FILE`**      | Null model file
| **`-dbsnp FILE`**     | dbSNP table file
| **`-out DIR`**        | Specify output directory name under currentg working folder.
| `-match`              | (optional) Boolean switch to turn on minor allele frequency based matching, default: off
| `-rep NUMBER`         | (optional) Number of permutation/simulation to be performed, default: 2000
| `-corr NUMBER`        | (optional) Bonferroni correction multiplier for multiple test, default: 1
| `-phenotype STRING`   | (optional) User-provided phenotype name, default: "."
| `-ancestry STRING`    | (optional) User provided ancestry name, default: "."

## How to cite

_fill this in whenever the publication details are available_

## Authors

| Name              | Email              | Institution                    |
|-------------------|--------------------|--------------------------------|
| Dr. Xiaoting Chen | chenxt@mail.uc.edu | Cincinnati Children's Hospital

[md]: https://help.github.com/categories/writing-on-github/
[gsl]: https://www.gnu.org/software/gsl/
