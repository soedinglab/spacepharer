# SpacePHARER: CRISPR Spacer Phage-Host pAiRs findER

SpacePHARER is a modular toolkit for sensitive phage-host interaction identification using CRISPR spacers. SpacePHARER adapts the fast homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2) to sensitively query short spacer sequences. It introduces a novel approach of aggregating sets of spacer-based hits to discover phage-host matches. SpacePHARER is GPLv3-licensed open source software implemented in C++ and available for Linux and macOS. The software is designed to run efficiently on multiple cores.

<p align="center"><img src="https://github.com/soedinglab/spacepharer/blob/master/.github/SpacePHARER.png" height="250"/></p>

## Installation

SpacePHARER can be used by compiling from source (see below) or downloading a statically compiled version. It requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).

     # install from bioconda
     conda install -c conda-forge -c bioconda spacepharer
     # pull docker container
     docker pull soedinglab/spacepharer
     # static SSE4.1 build
     wget https://mmseqs.com/spacepharer/spacepharer-linux-sse41.tar.gz; tar xvzf spacepharer-linux-sse41.tar.gz; export PATH=$(pwd)/spacepharer/bin/:$PATH
     # static AVX2 build
     wget https://mmseqs.com/spacepharer/spacepharer-linux-avx2.tar.gz; tar xvzf spacepharer-linux-avx2.tar.gz; export PATH=$(pwd)/spacepharer/bin/:$PATH

### Compile from source

Compiling SpacePHARER from source has the advantage of system-specific optimizations, which should improve its performance. To compile SpacePHARER `git`, `g++` (4.8 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the SpacePHARER binary will be located in the `build/bin` directory.

      git clone https://github.com/soedinglab/spacepharer.git
      cd spacepharer
      mkdir build
      cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j
      make install
      export PATH=$(pwd)/spacepharer/bin/:$PATH

:exclamation: If you want to compile SpacePHARER on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and SpacePHARER will not be able to run multithreaded. Adjust the `cmake` call above to:

      CC="$(brew --prefix)/bin/gcc-9" CXX="$(brew --prefix)/bin/g++-9" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## Input

SpacePHARER will conduct a similarity search between six-frame translated CRISPR spacer sequences and sets of phage **ORFs** (open reading frames), combine multiple evidences (**hits**) found between the two sets and predict prokaryote-phage pairs (**matches**) with strictly controlled **FDR** (false discovery rate). The starting point are FASTA files of nucleotide sequences (`.fasta` or `.fasta.gz`). Spacers should be provided in multiple FASTA files, each containing spacers from one genome. For spacers, SpacePHARER also accepts output files from the following common CRISPR array analysis tools: [PILER-CR](https://www.drive5.com/pilercr/), [CRT](http://www.room220.com/crt/), [MinCED](https://github.com/ctSkennerton/minced) (derived from CRT format) and [CRISPRDetect](http://crispr.otago.ac.nz/CRISPRDetect/predict_crispr_array.html). Phage genomes are supplied as separate FASTA files(one genome per file), or can be downloaded by `downloadgenome`(see below).

## Running SpacePHARER

### Main Modules

      easy-predict      Predict phage-host matches from multiFASTA and common spacer files (PILER-CR, CRISPRDetect and CRT)
      downloadgenome    Download GenBank (phage) genomes and create sequence database
      createsetdb       Create sequence database from FASTA input
      predictmatch      Predict host-phage matches
      parsespacer       Parse a file containing CRISPR array in supported formats (CRT,PILER-CR and CRISPRDetect)

### Important parameters

     --reverse-fragments      reverse AA fragments (ORFs) to generate control setDB
     --fdr                    false discovery rate cut-off to determine S_comb threshold of predictions
     --fmt                    output format for predictmatch. (0: short (only matches); 1: long (matches and hits); 2: long with nucleotide alignment)

### Quick start

To start, you need to create a database of the phage genomes `targetSetDB` and control sequences `targetSetDB_rev`, against which no true match is expected. Here we create control sequences by reversing the extracted ORFs.

      spacepharer createsetdb examples/GCA*.fna.gz targetSetDB tmpFolder
      spacepharer createsetdb examples/GCA*.fna.gz targetSetDB_rev tmpFolder --reverse-fragments 1

Alternatively, you can use `downloadgenome` to download a list of phage genomes. Here, the reversed control sequences are automatically created.

      # GenBank_phage_2018_09 is a set of nearly 8000 phage genomes
      spacepharer downloadgenome GenBank_phage_2018_09 targetSetDB tmpFolder
      
      # Alternatively you can pass a list of URLs to downloadgenome
      spacepharer downloadgenome examples/genome_list.tsv targetSetDB tmpFolder

The `easy-predict` workflow directly returns a tab-separated (.tsv) file containing phage-host predictions from (multiple) FASTA or supported CRISPR array files (CRT/MinCED, PILER-CR or CRISPRDetect) queries.

      spacepharer easy-predict examples/*.fas targetSetDB predictions.tsv tmpFolder

### Creating databases

Before search, query or target sequences contained in FASTA files need to be converted to database format by calling `createsetdb`. This command first creates a sequence DB, then extracts and translates all putative protein fragments (ORFs), and finally generates associated metadata. For spacer sequences, setting the parameter `--extractorf-spacer 1` is important to properly extract putative protein fragments from the spacers, which are usually a short partial ORF, not necessarily in frame.

      spacepharer createsetdb Query1.fasta [...QueryN.fasta] querySetDB tmpFolder --extractorf-spacer 1
      spacepharer createsetdb Target1.fasta [...TargetN.fasta] targetSetDB tmpFolder

You will also need to generate a control target set DB to allow SpacePHARER to calibrate the cutoff for reporting matches. SpacePHARER enables generating such control by reversing the protein fragments of your provided target DB using the parameter ```--reverse-fragments 1```:

      spacepharer createsetdb Target1.fasta [...TargetN.fasta] controlSetDB tmpFolder --reverse-fragments 1

### Downloading target genomes

As an alternative to creating target and control setDB, the `downloadgenome` module will download the provided list of URLs to phage genomes or a predefined list of phage genomes and create a target setDB in the provided path.

      spacepharer downloadgenome GenBank_phage_2018_09 targetSetDB tmpFolder
      
      # Generating control sequences can be disabled if a different set will be used
      spacepharer downloadgenome GenBank_phage_2018_09 targetSetDB tmpFolder --reverse-setdb 0

A list of predefined phage catalogues can be shown by executing `downloadgenome` without additional parameters:

      spacepharer downloadgenome

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

      ./minced prok.fasta prok.txt

CRISPRDetect:

CRISPRDetect is available as web server tool [here](http://crispr.otago.ac.nz/CRISPRDetect/predict_crispr_array.html).

### Searching and predicting matches

The `predictmatch` workflow gives more control of the execution of the prediction. Here a seperate control sequence set DB `controlSetDB` can be used. For example, we can assume that any spacer hit towards a eukoryota-targeting virus is a false positive:

    spacepharer downloadgenome GenBank_phage_2018_09 targetSetDB tmp --reverse-setdb 0
    spacepharer downloadgenome GenBank_eukvir_2018_09 controlSetDB tmp --reverse-setdb 0
    spacepharer predictmatch querySetDB targetSetDB controlSetDB outputFileName.tsv tmpFolder

### The SpacePHARER output

Upon completion, SpacePHARER outputs a tab-separated text file (.tsv). Each prokaryotic-phage match spans two or more lines:

    #prok_acc  phage_acc   S_comb      num_hits
    >spacer_acc      phage_acc   p_bh    spacer_start      spacer_end  phage_start phage_end   putative_5'_PAM|putative_3'_PAM
    *NUCL_SEQ_ALN_SPACER*
    *NUCL_SEQ_ALN_PHAGE*

The first line starts with `#`: prokaryotic accession, phage accession, combined score and number of hits in the match.

Each following line describes an individual hit: spacer accession, phage accession, p besthit, spacer start and end, phage start and end, putative 5’ PAM|putative 3’ PAM.

Optionally, the aligned spacer and phage sequences can be printed in two additional lines following each hit line, using `--fmt 2`

`--fmt 0` will output a short-format, if you wish to only see the match line.

### Removing temporary files

During the workflow execution, SpacePHARER will keep all intermediate outputs in `tmpFolder`, passing the `--remove-tmp-files` parameter will clear out the `tmpFolder` after workflows have finished.

## Hardware requirements

SpacePHARER will scale its memory consumption based on the available main memory of the machine. SpacePHARER needs a CPU with at least the SSE4.1 instruction set to run.
