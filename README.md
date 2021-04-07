# SpacePHARER: CRISPR Spacer Phage-Host pAiRs findER

SpacePHARER is a modular toolkit for sensitive phage-host interaction identification using CRISPR spacers. SpacePHARER adapts the fast homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2) to sensitively query short spacer sequences. It introduces a novel approach of aggregating sets of spacer-based hits to discover phage-host matches. SpacePHARER is GPLv3-licensed open source software implemented in C++ and available for Linux and macOS. The software is designed to run efficiently on multiple cores.

[Zhang, R., Mirdita, M., Levy Karin, E., Norroy, C., Galiez, C., & Söding, J. (2020). SpacePHARER: Sensitive identification of phages from CRISPR spacers in prokaryotic hosts. bioRxiv.](https://doi.org/10.1101/2020.05.15.090266)

<p align="center"><img src="https://github.com/soedinglab/spacepharer/blob/master/.github/SpacePHARER.png" height="250"/></p>

## Installation

SpacePHARER can be used by compiling from source (see below) or downloading a statically compiled version. It requires a 64-bit system. We recommend using a system with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux).

    # install from bioconda
    conda install -c conda-forge -c bioconda spacepharer
    # pull docker container
    docker pull soedinglab/spacepharer
    # static Linux AVX2 build
    wget https://mmseqs.com/spacepharer/spacepharer-linux-avx2.tar.gz; tar xvzf spacepharer-linux-avx2.tar.gz; export PATH=$(pwd)/spacepharer/bin/:$PATH
    # static Linux SSE4.1 build
    wget https://mmseqs.com/spacepharer/spacepharer-linux-sse41.tar.gz; tar xvzf spacepharer-linux-sse41.tar.gz; export PATH=$(pwd)/spacepharer/bin/:$PATH
    # static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
    wget https://mmseqs.com/spacepharer/spacepharer-osx-universal.tar.gz; tar xvzf spacepharer-osx-universal.tar.gz; export PATH=$(pwd)/spacepharer/bin/:$PATH

Precompiled binaries for other architectures (ARM64, PPC64LE) and very old AMD/Intel CPUs (SSE2 only) are available at [https://mmseqs.com/spacepharer](https://mmseqs.com/spacepharer).

### Compile from source

Compiling SpacePHARER from source has the advantage of system-specific optimizations, which should improve its performance. To compile SpacePHARER `git`, `g++` (4.9 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the SpacePHARER binary will be located in the `build/bin` directory.

    git clone https://github.com/soedinglab/spacepharer.git
    cd spacepharer
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j
    make install
    export PATH=$(pwd)/spacepharer/bin/:$PATH

:exclamation: If you want to compile SpacePHARER on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and SpacePHARER will not be able to run multithreaded. Adjust the `cmake` call above to:

    CC="$(brew --prefix)/bin/gcc-10" CXX="$(brew --prefix)/bin/g++-10" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## Input

SpacePHARER will conduct a similarity search between six-frame translated CRISPR spacer sequences and sets of phage **ORFs** (open reading frames), combine multiple evidences (**hits**) found between the two sets and predict prokaryote-phage pairs (**matches**) with strictly controlled **FDR** (false discovery rate). The starting point are FASTA files of nucleotide sequences (`.fasta` or `.fasta.gz`). Spacers should be provided in multiple FASTA files, each containing spacers from one genome. For spacers, SpacePHARER also accepts output files from the following common CRISPR array analysis tools: [PILER-CR](https://www.drive5.com/pilercr/), [CRT](http://www.room220.com/crt/), [MinCED](https://github.com/ctSkennerton/minced) (derived from CRT format) and [CRISPRDetect](http://crispr.otago.ac.nz/CRISPRDetect/predict_crispr_array.html). Phage genomes are supplied as separate FASTA files (one genome per file), or can be downloaded with `downloaddb`(see below).

## Running SpacePHARER

### Main Modules

* `easy-predict`      Predict phage-host matches from multiFASTA and common spacer files (PILER-CR, CRISPRDetect and CRT)
* `downloaddb`        Download spacers or GenBank phage genomes and create sequence database
* `createsetdb`       Create sequence database from FASTA input
* `predictmatch`      Predict host-phage matches
* `parsespacer`       Parse a file containing CRISPR array in supported formats (CRT,PILER-CR and CRISPRDetect)

### Important parameters

    --reverse-fragments      reverse AA fragments (ORFs) to generate control setDB
    --fdr                    false discovery rate cut-off to determine S_comb threshold of predictions
    --fmt                    output format for predictmatch. (0: short (only matches); 1: long (matches and hits); 2: long with nucleotide alignment)


### Quick start

To start, you need to create a database of the phage genomes `targetSetDB` and control sequences `targetSetDB_rev`, against which no true match is expected. Here we create control sequences by reversing the extracted ORFs.

    spacepharer createsetdb examples/GCA*.fna.gz targetSetDB tmpFolder
    spacepharer createsetdb examples/GCA*.fna.gz targetSetDB_rev tmpFolder --reverse-fragments 1

Alternatively, you can use `downloaddb` to download a list of phage genomes. Here, the reversed control sequences are automatically created.

    # GenBank_phage_2018_09 is a set of nearly 8000 phage genomes
    spacepharer downloaddb GenBank_phage_2018_09 targetSetDB tmpFolder
      
    # Alternatively you can pass a list of URLs to downloadgenome
    spacepharer downloaddb examples/genome_list.tsv targetSetDB tmpFolder

The `easy-predict` workflow directly returns a tab-separated (`.tsv`) file containing phage-host predictions from (multiple) FASTA or supported CRISPR array files (CRT/MinCED, PILER-CR or CRISPRDetect) queries.

    spacepharer easy-predict examples/*.fas targetSetDB predictions.tsv tmpFolder

### Creating databases

Before search, query or target sequences contained in FASTA files need to be converted to database format by calling `createsetdb`. This command first creates a sequence DB, then extracts and translates all putative protein fragments (ORFs), and finally generates associated metadata. For spacer sequences, setting the parameter `--extractorf-spacer 1` is important to properly extract putative protein fragments from the spacers, which are usually a short partial ORF, not necessarily in frame.

    spacepharer createsetdb Query1.fasta [...QueryN.fasta] querySetDB tmpFolder --extractorf-spacer 1
    spacepharer createsetdb Target1.fasta [...TargetN.fasta] targetSetDB tmpFolder

You will also need to generate a control target set DB to allow SpacePHARER to calibrate the cutoff for reporting matches. SpacePHARER enables generating such control by reversing the protein fragments of your provided target DB using the parameter ```--reverse-fragments 1```:

    spacepharer createsetdb Target1.fasta [...TargetN.fasta] controlSetDB tmpFolder --reverse-fragments 1
      
      
#### Downloading query CRISPR spacer sets

As an alternative to creating query setDB, you can use `downloaddb` to download a comprehensive set of CRISPR spacers. The query setDB will be automatically created in the provided path.

    # spacers_shmakov_et_al_2017 is a set of more than 30000 CRISPR spacer sets (Shmarkov et al., 2017)
    spacepharer downloaddb spacers_shmakov_et_al_2017 querySetDB tmpFolder
    
    # spacers_dion_et_al_2021 is a set of more than 490000 CRISPR spacer sets (Dion et al., 2021)
    spacepharer downloaddb spacers_dion_et_al_2021 querySetDB tmpFolder

#### Downloading target genomes

As an alternative to creating target and control setDB, the `downloaddb` module will download the provided list of URLs to phage genomes or a predefined list of phage genomes and create a target setDB in the provided path.

    spacepharer downloaddb GenBank_phage_2018_09 targetSetDB tmpFolder
      
    # Generating control sequences can be disabled if a different set will be used
    spacepharer downloaddb GenBank_phage_2018_09 targetSetDB tmpFolder --reverse-setdb 0

A list of predefined spacer or phage catalogues can be shown by executing `downloaddb` without additional parameters:

    spacepharer downloaddb

A file containing URLs can also be supplied to `downloaddb`:

    spacepharer downloaddb urls.txt targetSetDB tmpFolder
    # urls.txt content
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/836/905/GCA_000836905.1_ViralProj14035/GCA_000836905.1_ViralProj14035_genomic.fna.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/845/445/GCA_000845445.1_ViralProj14409/GCA_000845445.1_ViralProj14409_genomic.fna.gz

#### Adding taxonomic labels

If input spacers are supplied together with taxonomic identifiers, lowest common ancestors (LCA) of the phages are computed for each spacer. If input genomes are supplied together with taxonomic identifiers, the LCA of the host are computed for each phage.
Also, the [SpacePHARER output](#the-spacepharer-output) will contain taxonomic information for each match.

Databases download from the predefined entries in `downloaddb` come with taxonomic information already included. For custom databases, extra steps have to be taken:

#### Taxonomic labels in `createsetdb`

When calling `createsetdb` you can supply a tab-separated list of file names to NCBI taxonomy identifiers with the `--tax-mapping-file` parameter:

    # create phage database with taxonomic labels
    spacepharer createsetdb Target1.fasta [...TargetN.fasta] targetSetDB tmpFolder --tax-mapping-file targets.tsv
    # targets.tsv content (\t is a tab character)
    Target1.fasta\t10665
    ...
    TargetN.fasta\t31754

#### Taxonomic labels in `downloaddb URLFILE`

If an file containing URLs is supplied to `downloaddb` taxonomic labels can be provided in the second column of the URL file:

A file containing URLs can also be supplied to `downloaddb`:

    spacepharer downloaddb urls.txt targetSetDB tmpFolder
    # urls.txt content (\t is a tab character)
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/836/905/GCA_000836905.1_ViralProj14035/GCA_000836905.1_ViralProj14035_genomic.fna.gz\t10679
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/845/445/GCA_000845445.1_ViralProj14409/GCA_000845445.1_ViralProj14409_genomic.fna.gz\t244310

### Parsing spacer files

If you wish to provide spacer files (CRT/MinCED, PILER-CR or CRISPRDetect) as query, `parsespacer` will extract spacer sequences from each file and create a sequence DB. Then use `createsetdb` to generate associated metadata.

    spacepharer parsespacer spacerFile1.txt [...spacerFileN.txt] queryDB 
    spacepharer createsetdb queryDB querySetDB tmpFolder --extractorf-spacer 1

#### Sample commands for running spacer extraction tools

Below are sample commands for common spacer extraction tools whose format is accepted by SpacePHARER. If you use any of these tools, it is advisable to follow any updates to their commands on their user manuals.

PILER-CR: Use of `-noinfo` is mandatory. Otherwise, the format cannot be read.

    pilercr -noinfo -quiet -in prok.fasta -out prok.txt

CRT:

    java -cp CRT1.2-CLI.jar crt prok.fasta prok.txt

MinCED:

    minced prok.fasta prok.txt

CRISPRDetect is only available as a [web server](http://crispr.otago.ac.nz/CRISPRDetect/predict_crispr_array.html).

### Searching and predicting matches

The `predictmatch` workflow gives more control of the execution of the prediction. Here a seperate control sequence set DB `controlSetDB` can be used. For example, we can assume that any spacer hit towards a eukoryota-targeting virus is a false positive:

    spacepharer downloaddb GenBank_phage_2018_09 targetSetDB tmp --reverse-setdb 0
    spacepharer downloaddb GenBank_eukvir_2018_09 controlSetDB tmp --reverse-setdb 0
    spacepharer predictmatch querySetDB targetSetDB controlSetDB outputFileName.tsv tmpFolder

### The SpacePHARER output

Upon completion, SpacePHARER outputs a tab-separated text file (`.tsv`). Each prokaryotic-phage match spans two or more lines:

    #prok_acc  phage_acc   S_comb      num_hits
    >spacer_acc      phage_acc   p_bh    spacer_start      spacer_end  phage_start phage_end   5'_PAM|3'_PAM    5'_PAM|3'_PAM(reverse strand)
    *NUCL_SEQ_ALN_SPACER*
    *NUCL_SEQ_ALN_PHAGE*

The first line starts with `#`: prokaryotic accession, phage accession, combined score and number of hits in the match.

Each following line describes an individual hit: spacer accession, phage accession, p besthit, spacer start and end, phage start and end, possible 5’ PAM|3’ PAM, possible 5’ PAM|3’ PAM on the reverse strand.

Optionally, the aligned spacer and phage sequences can be printed in two additional lines following each hit line, using `--fmt 2`

`--fmt 0` will output a short-format, if you wish to only see the match line.

If the phage database was created with taxonomic labels, a result file with the suffix `_lca.tsv` is also created. The first column of this file contains the spacer accession and the remaining columns are described in the [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-output-and-tsv).

If the spacer database was created with taxonomic labels, a result file with the suffix `_lca_per_target.tsv` is also created. The first column of this file contains the phage genome file name and the remaining columns are described in the [same MMseqs2 wiki entry](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-output-and-tsv). Additionally, each match in the base `.tsv` will also contain these columns.

### Removing temporary files

During the workflow execution, SpacePHARER will keep all intermediate outputs in `tmpFolder`, passing the `--remove-tmp-files` parameter will clear out the `tmpFolder` after workflows have finished.

## Hardware requirements

SpacePHARER will scale its memory consumption based on the available main memory of the machine. SpacePHARER needs a 64-bit CPU with at least the SSE2 instruction set to run.
