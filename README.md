# fastqtiler

This is a program for read quality visualization of Illumina sequencing runs.

Installation
============

* Clone the repository:

		git clone https://github.com/har-wradim/fastqtiler.git

or download the files as [zip](https://github.com/har-wradim/fastqtiler/archive/master.zip)

* Optionally copy the scripts somewhere to your PATH:

		sudo cp fastqtiler/*.pl /usr/local/bin

* For the `fastqtiler` program to work ensure that [ImageMagick](http://www.imagemagick.org/script/index.php) is installed and `convert` is available in the PATH.

Launching the program
=====================

Depending on how you installed the program you can either launch it as `fastqtiler.pl`, `./fastqtiler.pl` or `perl fastqtiler.pl`.

The general form of the command looks like `fastqtiler.pl [parameters] [files]`

Parameters
==========

| **Short** | **Long**      |               **Description**                |     **Default value**     |
|-----------|---------------|----------------------------------------------|:-------------------------:|
| `-k`      | `--kit`       | sequencing **k**it/flowcell type (see below) | none                      |
| `-n`      | `--name`      | **n**ame for output files                    | tiles                     |
| `-o`      | `--outdir`    | **o**utput directory                         | .                         |
| `-t`      | `--threads`   | max number of **t**hreads                    | 2                         |
| `-l`      | `--len`       | read **l**ength                              | taken from the first read |
| `-f`      | `--frac`      | **f**raction of the reads to analyze         | 0.01                      |
| `-s`      | `--scale`     | **s**cale factor (pixel size in FC units)    | 1000                      |
| `-g`      | `--goodcolor` | color for the high-quality bases             | magenta                   |
| `-b`      | `--badcolor`  | color for the low-quality bases              | green                     |
| `-i`      | `--lanemin`   | lane number reads start at                   | 1                         |
| `-j`      | `--lanemax`   | lane number reads end   at                   | kit-specific              |
| `-h`      | `--help`      | output the help page                         |                           |
| `-v`      | `--version`   | same as above                                |                           |

The `--kit` options is mandatory. The supported kits are as follows:

| **Kit**                         | **Full value**               | **Alternative shortcuts**                                                  |
|---------------------------------|------------------------------|----------------------------------------------------------------------------|
| HiSeq 2000 or 2500, High Output | `hiseq2000/2500:high-output` | `hiseq2000:high-output`, `hiseq2500:high-output`, `hs2000:ho`, `hs2500:ho` |
| HiSeq 2000 or 2500, Rapid Run   | `hiseq2000/2500:rapid-run`   | `hiseq2000:rapid-run`, `hiseq2500:rapid-run`, `hs2000:rr`, `hs2500:rr`     |
| HiSeq 3000 or 4000              | `hiseq3000/4000`             | `hiseq3000`, `hiseq4000`, `hs3000`, `hs4000`                               |
| HiSeq X                         | `hiseqX`                     | `hsX`                                                                      |
| MiSeq                           | `miseq`                      | `ms`                                                                       |
| MiSeq Micro                     | `miseq:micro`                | `ms:micro`                                                                 |
| MiSeq Nano                      | `miseq:nano`                 | `ms:nano`                                                                  |
| NextSeq 500 or 550, High Output | `nextseq500/550:high-output` | `nextseq500:high-output`, `nextseq550:high-output`, `ns500:ho`, `ns550:ho` |
| NextSeq 500 or 550, Mid Output  | `nextseq500/550:mid-output`  | `nextseq500:mid-output`, `nextseq550:mid-output`, `ns500:mo`, `ns550:mo`   |

If you don't know the kit used you can always run the `fastqtester.pl` script providing it your file(s) (one exemplary file would suffice). It will also output some other useful information

Input files
===========

Either gzipped, bzipped or raw fastq, but not a mix of different types.

Output files
============

As output you will get a folder with the images corresponding to the quality scores at each one of the sequencing cycles and an html report. Here is how a [good run](http://alephreish.github.io/fastqtiler/good.html) looks like and here is a [bad one](http://alephreish.github.io/fastqtiler/bad.html).
