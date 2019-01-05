## Regulatory Element Locus Intersection (RELI) Analysis

RELI (Regulatory Element Locus Intersection) is an algorithm for discovering
transcription factors (TFs) that bind a significant number of loci associated
with a given disease or phenotype (e.g., through a Genome Wide Association
study, or GWAS).

The major data components are

1. An input set of disease or phenotype-associated genetic variants (RS IDs)

2. An internal “library” consisting of many ChIP-seq dataset peaks (in the form
   of .bed files)

3. An internal file containing information on genetic variant allele
   frequencies, etc.

To assess the significance of the intersection between the input disease
variants and a given TF ChIP-seq dataset, RELI performs simulations, generating
a null distribution used for P-value calculations.

The output of RELI is a series of statistics based upon the significance of the
overlap between the input genetic variants and the selected ChIP-seq dataset.

Additional details on RELI and the associated findings can be found in its
[accompanying publication][pubmed].

## Installation

RELI requires a C++11 compiler (_e.g._ GNU CC 4.7 or higher) and `libgsl` and
`libgslcblas` from the [GNU Scientific Library][gsl].

You may download the latest release as a compressed archive from GitHub, or
clone the repository with Git:

    # GitHub
    git clone https://github.com/WeirauchLab/RELI.git

    # Weirauch Lab GitLab
    git clone https://tfwebdev.research.cchmc.org/gitlab/ches2d/RELI_public.git

### Installation with Docker

This is the recommended method if you have some familiarity with Docker, as it
does not require you to download a compiler or any of the dependencies
necessary to build RELI from source.

You will need to install the appropriate [Docker client][dockercli] for your
OS. Please see the [official docs][dockerinst] for help with that.

The most straightforward way to get started is to simply use the pre-built
image available on [Docker Hub][hubi]:

```
docker run -it --rm weirauchlab/reli RELI --help
```

Or, if you've cloned the `RELI_public` repository (see above), you can locally
build the CentOS 7-based Docker container and compile RELI from source as
follows:

```
cd /path/to/cloned/repo
docker build -t reli .

# test to see if it works
docker run -it --rm reli RELI --help
```

### Installation on Linux

A GNU-style `Makefile` is provided in the repository. With GSL installed
system-wide, you can build the RELI binary with just

    make

then run `./RELI` with no arguments to verify that you have a working binary
(you should get a help screen).


In order to run a test analysis, **you need to download the sample data**
either manually (see the next section) or just type

    make test

which will download and validate the sample datasets automatically, then invoke
`example/example_run.sh` to invoke RELI on the sample data.

**This test analysis requires around 10 GB of RAM** to finish successfully; 16 GB
is recommended.

#### Toolchain or libraries in non-standard locations

The included `Makefile` will respect `CFLAGS` and `LDFLAGS` if set in the
environment, for example, if you have a locally-built GSL that is installed in
a non-standard place (such as in your home directory):

    CFLAGS=-I/path/to/include LDFLAGS=-L/path/to/lib make

If `g++` is not available in your `PATH` (or it has a different name), you will
likely want to modify the Makefile directly, beginning around line 33 with the
`CC` variable.

### Other platforms

RELI has also been verified to build and run on the following platforms (in
addition to GNU/Linux):

* Windows with [Cygwin][cyg] and GCC 5.4.0 (ensure the `gcc-g++`, `make`,
  `gsl`, and `libgsl-devel`, and `curl` packages are installed, at a minimum)

* Mac OS X 10.11.6 (El Capitan) with LLVM 8.0.0 (provided by the
  [Xcode Command Line Tools][xct]) and GSL installed from [MacPorts][mp] 

On Windows, make sure you run `make` (or the `example/example_run.sh` script)
from within the Cygwin shell, _not_ the Windows Command Prompt or PowerShell.
You may need to lightly modify the CDT build toolchain settings if your
installation of Cygwin is not at `C:\Cygwin64`.

[Eclipse CDT][cdt] project settings files are also included for both of the
above toolchains. Just create a copy (or symlink) of the appropriate one
called `.cproject`, then choose _File_ &rarr; _Import..._ &rarr; _Existing
Projects into Workspace_ and browse to where you cloned the repository.

### Manually downloading sample data

If you have problems with `make test` (perhaps you don't have `curl`
available), you can manually download and extract the sample datasets from

> <https://tf.cchmc.org/external/RELI/RELI_public_data.tar.bz2>

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

To use a different genome build, use the UCSC [`fetchChromSizes`][fcs] utility
(usage information [here][fcsusage]) to download chromosome information for that
build. You may wish to prune lines representing unmapped chromosome information
(_e.g._, `chrN_glXXXXXX_random` and `chrUn_glXXXXXX`) from the downloaded data
file.

Be advised, however, that the null model included with the data was generated
for _Homo sapiens_ at build hg19; using a later "hg" build may invalidate this
model.

If you need support for a different organism, please contact us via email for
additional details (see "Feedback" section, below), or [file an issue][ghi]
against the public GitHub repository.

## How to cite

[_Transcription factors operate across disease loci, with EBNA2 implicated in
autoimmunity._][pubmed]

Harley JB, Chen X, Pujato M, Miller D, Maddox A, Forney C, Magnusen AF, Lynch A,
Chetal K, Yukawa M, Barski A, Salomonis N, Kaufman KM, Kottyan LC, Weirauch MT.

_Nat Genet._ 2018 Apr 16. doi: [10.1038/s41588-018-0102-3][doi]. [Epub ahead of
print]

PMID: 29662164

## Feedback

Please report any issues with RELI (or feature suggestions) in our
[GitHub issue tracker][ghi].

With other questions, you may contact [Dr. Chen][xiaoting] (the primary author
of RELI) or [Dr. Weirauch][matt] via email.

## Contributors

| Name                 | Institution                    | Remarks
|----------------------|--------------------------------|----------------------
| Dr. Xiaoting Chen    | Cincinnati Children's Hospital | _primary author_

## Credits

Project avatar based on Wikimedia Commons [Chromosome_18.svg][wmc]

[md]: https://help.github.com/categories/writing-on-github/
[gsl]: https://www.gnu.org/software/gsl/
[dockercli]: https://hub.docker.com/search/?type=edition&offering=community
[dockerinst]: https://docs.docker.com/install/
[hubi]: https://hub.docker.com/r/weirauchlab/reli
[wmc]: https://commons.wikimedia.org/wiki/File:Chromosome_18.svg
[fcs]: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
[fcsusage]: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
[pubmed]: https://www.ncbi.nlm.nih.gov/pubmed/29662164
[doi]: https://dx.doi.org/10.1038/s41588-018-0102-3
[ghi]: https://github.com/WeirauchLab/RELI/issues
[mp]: https://www.macports.org/
[cdt]: https://www.eclipse.org/cdt/
[cyg]: http://cygwin.org/
[xct]: https://stackoverflow.com/a/32894314
[xiaoting]: mailto:Xiaoting.Chen%20-at-%20cchmc.org?subject=RELI%20feedback&body=(don't%20forget%20to%20change%20%22-at-%22%20to%20%22%40%22%20before%20sending!)
[matt]: mailto:Matthew.Weirauch%20-at-%20cchmc.org?subject=RELI%20feedback&body=(don't%20forget%20to%20change%20%22-at-%22%20to%20%22%40%22%20before%20sending!)
