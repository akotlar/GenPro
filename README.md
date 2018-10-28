GenPro
======

Make personal protein databases using next generation sequencing data.

## Installation

Installation on a Unix/linux distribution (including OS X) is straightforward
and enumerated below. To simplify this further (and offer a solution that
should work on most platforms), we have included a docker file that will setup
a docker container with all necessary dependencies installed.

### With [Docker](https://www.docker.com/)

First, [install the community edition docker application](https://www.docker.com/community-edition)
for your system.

Next, clone the github repo, which has the latest docker image for genpro and
use the docker image to build a container with all dependencies (and GenPro
installed).

```
git clone https://github.com/wingolab-org/GenPro.git
cd GenPro
docker build -t genpro .
```

There are a few options for running GenPro using the docker container. The
simplest way is to login to the container using the command below. For this to
be useful, we need to allow the docker container to write to the host system.
This is done with the `-v` option, as in `-v /your/machine:/docker/container`.
You will need to set the host directory accordingly and note that in the example
below it is set to the working directory where you launch the `genpro` docker
container.

```
docker run -it -v $(pwd):/GenProData genpro bash
```

### Installing on a unix/linux

The approach is:
1. Install [local::lib](https://metacpan.org/pod/local::lib), which allows 
installation of GenPro (and other perl packages) without using `sudo`.
2. Install  [App::cpanminus](https://metacpan.org/pod/App::cpanminus), which
simplifies and automates the installation of perl packages.
3. Clone the GenPro repository.
4. Install the prebuilt GenPro package in the repository using `cpanm`.

Install [`local::lib`](https://metacpan.org/pod/local::lib) by downloading the
latest tarball and unpack it. See the
[bootstrapping](https://metacpan.org/pod/local::lib) section of the
documentation.

Example:
```
curl -O https://cpan.metacpan.org/authors/id/H/HA/HAARG/local-lib-2.000019.tar.gz
tar xzvf local-lib-2.000019.tar.gz
cd local-lib-2.000019
perl Makefile.PL --bootstrap
make test && make install
```

Install [`App::cpanminus`](https://metacpan.org/pod/App::cpanminus).

Example:
```
curl -L http://cpanmin.us | perl - App::cpanminus
```

Clone the GenPro repository with [`git`](https://git-scm.com). Install GenPro with 
`cpanm`, which will automatically install any needed dependencies.

Example:
```
git clone https://github.com/wingolab/GenPro.git
cd GenPro
cpanm GenPro.tar.gz
```

An alternative approach is clone the repository, unpack the GenPro.tar.gz
tarball and install it manually. Unless `local::lib` is installed, this 
approach will require using `sudo`.

Example:
```
git clone https://github.com/wingolab/GenPro.git
cd GenPro
tar xzvf GenPro.tar.gz
cd GenPro
perl Makefile.PL
make test
make install    
```

## Usage

### Memory requirements

The two programs that require the most memory are `GenPro_create_db.pl` and
`GenPro_make_perprotdb1.pl`. For `GenPro_create_db.pl`, the memory 
requirement scales with genome size and gene density.
For `GenPro_make_perprotdb1.pl`, the size of the reference protein database and
the number of samples raise the memory requirement.

The table below gives representative memory consumption. The WGS samples were 
obtained from 1000genomes phase1 (hg19) `vcf` files and converted to the `snp` 
format for GenPro using `vcfToSnp`.

```
--------------------------------------------------------------------------------
program                   | perl    | memory use | processed
--------------------------------------------------------------------------------
GenPro_create_db.pl       | v5.16.3 | 18.8 Gb    | hg38, chromosome 1
--------------------------------------------------------------------------------
GenPro_make_refprotdb.pl  | v5.16.3 | 795.7 Mb   | hg38, chromosome 1
--------------------------------------------------------------------------------
GenPro_make_perprotdb1.pl | v5.16.3 | 21.5 Gb    | 50 HapMap WGS, phase 1
--------------------------------------------------------------------------------
GenPro_make_perprotdb2.pl | v5.16.3 | 603 Mb     | HapMap NA06994, phase 1
--------------------------------------------------------------------------------
```

### 1. Download genomic data for a particular organism.

```
GenPro_download_ucsc_data.pl -d hg38 -g hg38
```

The above command will perform a dry-run download of hg38 (genome and annotated
gene coordinates). It relies on `rsync` being installed, which should be present
on unix, linux, and OS X by default. In the example, the data will be downloaded
into `hg38` directory, which may be created if it did not already exist.
`GenPro_download_ucsc_data.pl` will download knownGenes track and the genome of
the organism by default. By default, `GenPro_download_ucsc_data.pl` is set to a
dry-run (i.e., no download). Use the `--act` or `-a` switch to "act", i.e.,
download the data. Take care when using it since it since it will download from
a remote server.


### 2. Generate a binary index of the genome for the organism.

- Use `GenPro_create_db.pl`.
- There are helper scripts `sh/runall_create_db.sh` that work with SGE to
  build all chromosomes on a cluster with SGE. For example:
```
qsub -v USER -v PATH -cwd -t 1-26 runall_create_db.sh <genome>
```
- Alternatively, iterate over all the chromosomes. For example:
```
for ((i=1;i<27;i++)); do
  GenPro_create_db.pl -g hg38 -c $i \
    --genedir hg38/gene             \
    --geneset knownGene             \
    -f hg38/chr -o hg38/idx;
done;
```

### 3.  Make reference proteins for the organism.
- `GenPro_make_refprotdb.pl` uses the indexed binary annotation and creates
  all proteins for a given chromosome.
- A helper script `runnall_make_refprot.sh` automates this for SGE. For example,
```
qsub -v USER -v PATH -cwd -t 1-26 runnall_make_refprot.sh \
  <genome>              \
  <binary genome index> \
  <output dir>
```

### 4.  Make the personal proteins for the sample.
- This is a 2-step process that uses `GenPro_make_perprotdb1.pl` and
  `GenPro_make_perprotdb2.pl`.
- `GenPro_make_perprotdb1.pl` creates a per chromosome database of all relevant
  (i.e., nonsense/missense variants) for each sample in the snp file.
- `GenPro_make_perprotdb2.pl` creates a finished personal protein database for
  each sample. It provides two outputs:
    1. A json-encoded file that enumerates the variant protein information
    2. A fasta file with all full-length reference and variant proteins, which
    may be used as input for a proteomics search program.
-  Creating the final db is on a per sample basis, but it can take quite a
  while. It depends on how many proteins have multiple variants. It is especially
  slow for proteins with >10 variants. For example, if you have a protein with
  20 substitutions, there are 20! permutations. By design, all 20! proteins
  will be considered and only the proteins that contribute unique peptides will
  be retained. Any protein with more than 20 sites will have all variants
  inserted into the reference protein without performing any permutation.
- *To begin*, you will need to have genotype calls in the `snp` file format. To 
convert `vcf` to `snp` format you will need 
[bcftools](https://samtools.github.io/bcftools/) installed. The helper script, 
`bin/vcfToSnp`, calls `bcftools` internally so `bcftools` will need to be in 
your path.
