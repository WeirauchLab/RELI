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

_I recommend explaining whether options are required or optional here, like
the example I used for `-match` (not sure if that's really optional or not
but it's just an example)._

_It's possible to add another column to the table if you want; the complete
reference for "GitHub-flavored Markdown" is [here][md]._

| Option                       | Explanation                                 |
|------------------------------|---------------------------------------------|
| `-snp <snp_file.snp>`        | description of this option                  |
| `-ld <ld_file.ld>`           |                                             |
| `-index path/to/index`       |                                             |
| `-data path/to/datafile`     |                                             |
| `-target path/to/target`     |                                             |
| `-build path/to/hg19.txt`    |                                             |
| `-null path/to/null/model`   |                                             |
| `-dbsnp path/to/snptable`    |                                             |
| `-out path/output/dir`       |                                             |
| `-match`                     | _(optional)_ description of this option     |
| `-rep <num>`                 |                                             |
| `-corr <num>`                |                                             |
| `-phenotype path/to/pheno`   |                                             |
| `-ancestry path/to/ancestry` |                                             |

## How to cite

_fill this in whenever the publication details are available_

## Authors

| Name              | Email              | Institution                    |
|-------------------|--------------------|--------------------------------|
| Dr. Xiaoting Chen | chenxt@mail.uc.edu | Cincinnati Children's Hospital

[md]: https://help.github.com/categories/writing-on-github/
[gsl]: https://www.gnu.org/software/gsl/
