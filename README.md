# RELI analysis

A one-paragraph description of what the tool does.

## Installation

RELI requires a C++11 compiler (_e.g._ GNU CC 4.7 or higher) and `libgsl` and
`libgslcblas` from the [GNU Scientific Library][gsl] library.

A build script (for Windows) is included; simply clone this repository

    git clone https://tfwebdev.research.cchmc.org/gitlab/ern6xv/RELI_public.git

and run

    src\compile.bat

from the Windows Command Prompt.

_FIXME: This obviously won't work on Windows unless you have something like
Cygwin installed; we either need to say that Cygwin or MinGW (or similar) is
required, or provide a Visual Studio `.vcxproj` file which will build using
the latest Visual Studio Community Edition._

### Downloading sample data

Download and extract the sample datasets from

    https://tf.cchmc.org/external/RELI/data.zip

such that the decompressed data is inside a `data` subdirectory, within the `RELI_public` repository you cloned above.

Finally, change into the `example` directory and run `Example_run.bat` from the Command Prompt:

    c:\path\to\repo> cd example
    c:\path\to\repo> Example_run


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
