# VODKA2
<!-- GETTING STARTED -->
## Getting Started

This is the VODKA2.0 version that was used to analyze data in E. Achouri, SA. Felt, M. Hackbart, NS. Rivera-Espinal and CB. Lopez in 'VODKA2: A fast and accurate method to detect non-standard viral genomes from large RNA-seq datasets' (submitted)

### Prerequisites

* perl 5
* bowtie2-v2.4.1 or higher
* blast-v2.11.0 or higher

<!-- USAGE EXAMPLES -->
## Usage

### Generate VODKA2 database

 1. DB for cbVG analysis

    perl VODKA2_genomefa_to_newfasta_CB_v2.pl <genomefa> <nb_from_right> <read_length> <newfasta_name>

* \<genomefa\>  Strandard virus reference genome in FASTA format (1 line sequence)
* \<nb_from_right\> Number of nt from right end of the sequence (use genome size for a full length analysis)
* \<read_length\> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
* \<newfasta_name\> Name of the multifata file to be generated. Should be as <virus>.<nb_from_right>.<read_length>.fasta

 2. DB for delVG analysis

    perl VODKA2_genomefa_to_newfasta_DEL_v2.pl <genomefa> <nb_from_right> <read_length> <newfasta_name> <gap_size>

* \<genomefa\>  Strandard virus reference genome in FASTA format (1 line sequence)
* \<nb_from_right\> Number of nt from right end of the sequence (use genome size for a full length analysis)
* \<read_length\> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
* \<newfasta_name\> Name of the multifata file to be generated. Should be as <virus>.DEL.<nb_from_right>.<read_length>.fasta
* \<gap_size\> e.g. 10

 3. Build bowtie2 index

    bowtie2-build --large-index <newfasta_name> <b2_idx>

* \<newfasta_name\> Name of the multifasta DB (output from previous step)
* \<bt2_idx\>  Name of VODKA2 DB index for Bowtie2


### Run VODAK2 analysis setup script

 1. Setup analysis folder and scripts

    bash VODKA2_analysis_setup.sh [-h] -f <files.txt> -d <db_bt2_idx> -v <ref.fasta> -p <project> [options]

* -f \<files.txt\>  TXT file containing the list of fastq files
* -d \<bt2_idx\>  Name of VODKA2 DB index for Bowtie2
* -v \<ref.fasta\> Virus reference genome in FASTA format
* -p \<project\>  Project name used to label output folders and files


 2. Run analysis

    bash cbVG_VODKA2_analysis_<project>.sh
or
    bash delVG_VODKA2_analysis_<project>.sh


<!-- CONTACT -->
## Contact

Emna Achouri - emna@wustl.edu

Carolina López - clopezzalaquett@wustl.edu

Project Link: [https://github.com/lopezlab-washu/VODKA2](https://github.com/lopezlab-washu/VODKA2)

<p align="right">(<a href="#top">back to top</a>)</p>

