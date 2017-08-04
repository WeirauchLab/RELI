## Regulatory Element Locus Intersection (RELI) Analysis

A one-paragraph description of what the tool does.

## Installation

RELI requires a C++11 compiler (_e.g._ GNU CC 4.7 or higher) and `libgsl` and
`libgslcblas` from the [GNU Scientific Library][gsl].

A GNU-style `Makefile` is provided. You can build the RELI binary with just

    make

In order to run a test analysis, **you need to download the example data**
either manually (see the next section) or just type

    make test

which will download and validate the sample datasets automatically, then run
a the test analysis in `example/example_run.sh`.

### Manually downloading sample data

If you have problems with `make test` (perhaps you don't have `curl`
available), you can manually download and extract the sample datasets from

> <https://tf.cchmc.org/external/RELI/data.tar.bz2>

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

## Adding custom data sets

### Adding ChIP-seq data

To add an additional ChIP-seq dataset, create an entry in the ChIP-seq index
file (`data/ChIPseq.index`) with the following tab-delimited format:

    label ⇥  source ⇥  Cell ⇥  TF ⇥  Cell label ⇥  PMID ⇥  Group ⇥  EBV Status ⇥  Species

where `label` corresponds to the filename, which you should deposit in the
`data/ChIP-seq` directory (in BED 4 column format).

### Adding a new genome build or species

To use a different genome build, use the UCSC [`fetchChromSizes` utility][fcs]
(usage information [here][fcsusage]) to download chromosome information for that
build. You may wish to prune lines representing unmapped chromosome information
(_e.g._, `chrN_glXXXXXX_random` and `chrUn_glXXXXXX`) from the downloaded data
file.

Be advised, however, that the null model included with the data was generated
for _Homo sapiens_ at build hg19; using a later "hg" build may invalidate this
model.

Please contact us via email (or [file an issue][ghi] against the public GitHub
repository) for additional details, or if you need support for a different
organism.

## How to cite

_fill this in whenever the publication details are available_

## Authors

| Name              | Email              | Institution                    |
|-------------------|--------------------|--------------------------------|
| Dr. Xiaoting Chen | chenxt@mail.uc.edu | Cincinnati Children's Hospital

## Credits

Project avatar based on Wikimedia Commons [Chromosome_18.svg][wmc]

[md]: https://help.github.com/categories/writing-on-github/
[gsl]: https://www.gnu.org/software/gsl/
[wmc]: https://commons.wikimedia.org/wiki/File:Chromosome_18.svg
[fcs]: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
[fcsusage]: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
[ghi]: https://github.com/WeirauchLab/RELI/issues