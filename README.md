# VODKA2
<!-- GETTING STARTED -->
## Getting Started

This is the VODKA2.0 version that was used to analyze data in E. Achouri, S. Felt, M. Hackbart and CB. Lopez in 'VODKA2: An accurate method to detect copy-back and deletion viral genomes from next-generation sequencing data' (submitted)

### Prerequisites

* blast-v2.11.0 or higher
* bowtie2-v2.4.1
* perl 5


<!-- USAGE EXAMPLES -->
## Usage

1- Generate VODKA2 multifata copy-back viral genome database

    perl VODKA2_genomefa_to_newfasta_CB_v2.0b.pl <genomefa> <nb_from_right> <read_length> <newfasta_name>

* \<genomefa\>  Strandard virus reference genome in FASTA format (1 line sequence)
* \<nb_from_right\> Number of nt from right end of the sequence (use genome size for a full length analysis)
* \<read_length\> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
* \<newfasta_name\> Name of the multifata file to be generated. Should be as <virus>.<nb_from_right>.<read_length>.fasta

2- Build bowtie2 index

    bowtie2-build --large-index <newfasta_name> <b2_idx>

* \<newfasta_name\> Name of the multifata file to be generated. Should be as <virus>.<nb_from_right>.<read_length>.fasta
* \<bt2_idx\>  Name of VODKA2 DB index for Bowtie2

3- Run VODAK2 analysis setup script

    bash VODKA2demo_analysis_setup.sh [-h] -f <files.txt> -d <bt2_idx> -p <project> [options]

* -f \<files.txt\>  TXT file containing the list of fastq files (1 per sample, with the relative path)
* -d \<bt2_idx\>  Name of VODKA2 DB index for Bowtie2 (make sure already exists within VODKA installation folder)
* -p \<project\>  Project name used to label output folders and files

options:
  * -i \<folder\> VODKA2-v2.0b installation folder (default: ./VODKA2-v2.0b-master)
  * -s \<suffix\> Suffix of FASTQ files (anything after samplename, default: _nonHost.fq.gz)
  * -l \<on|off\> LSF on or off (default: on)
  * -h  Show this help text and exit

4- Run VODKA2 analysis

    bash cbVG_VODKA2_analysis_<project>.sh

* replace \<project\> with the project name specified in step 3.

5- Generate extra report by sample

    bash vodka_extra-info_v2.sh -s <samplename> -p <project> -t CB -n <N> -a <genome_annotation.txt>

* \<samplename\>  Sample name
* \<project\> Project name as specified in steps 3. and 4.
* \<N\> Number of nucletides to be used around break and rejoin positions for species aggregation
* \<genome_annotation.txt\> Reference genome annotation file

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

Emna Achouri - emna@wustl.edu

Carolina López - clopezzalaquett@wustl.edu

Project Link: [https://github.com/lopezlab-washu/VODKA2](https://github.com/lopezlab-washu/VODKA2)

<p align="right">(<a href="#top">back to top</a>)</p>

