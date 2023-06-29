# VODKA2
<!-- GETTING STARTED -->
## Getting Started

This is the VODKA2.0 version that was used to analyze data in E. Achouri, SA. Felt, M. Hackbart, NS. Rivera-Espinal and CB. Lopez in 'VODKA2: A fast and accurate method to detect non-standard viral genomes from large RNA-seq datasets' (submitted)

### Prerequisites

* perl 5
* bowtie2-v2.4.1 or higher
* blast-v2.11.0 or higher
* R-v4.3+

or

* docker (docker container image is available on docker hub repository)

### Installation

Option 1: Insatll all required softwares mentionned above and download VODKA2 scripts

Option 2: Run the VODKA2 docker container
    <code>docker pull emna20/vodka2</code>

<!-- USAGE EXAMPLES -->
## Usage

### Generate VODKA2 database

 1. DB for cbVG analysis

    <code>perl VODKA2_genomefa_to_newfasta_CB_v2.pl <genomefa> <nb_from_right> <read_length> <newfasta_name></code>

* \<genomefa\>  Strandard virus reference genome in FASTA format (1 line sequence)
* \<nb_from_right\> Number of nt from right end of the sequence (use genome size for a full length analysis)
* \<read_length\> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
* \<newfasta_name\> Name of the multifata file to be generated. Should be as <virus>.<nb_from_right>.<read_length>.fasta

 2. DB for delVG analysis

    <code>perl VODKA2_genomefa_to_newfasta_DEL_v2.pl <genomefa> <nb_from_right> <read_length> <newfasta_name> <gap_size></code>

* \<genomefa\>  Strandard virus reference genome in FASTA format (1 line sequence)
* \<nb_from_right\> Number of nt from right end of the sequence (use genome size for a full length analysis)
* \<read_length\> Read length of samples intended to be analysed using the current db (can be greater value than the actual sample read length)
* \<newfasta_name\> Name of the multifata file to be generated. Should be as <virus>.DEL.<nb_from_right>.<read_length>.fasta
* \<gap_size\> e.g. 10

 3. Build bowtie2 index

    <code>bowtie2-build --large-index <newfasta_name> <b2_idx></code>

* \<newfasta_name\> Name of the multifasta DB (output from previous step)
* \<bt2_idx\>  Name of VODKA2 DB index for Bowtie2

    <code>bowtie2-build \<genomefa\> <b2_idx></code>

* \<genomefa\> Strandard virus reference genome in FASTA format (output from previous step)
* \<bt2_idx\>  Name of Bowtie2 index (keep the extension .fasta)


### Run VODAK2 analysis setup script

 1. Setup analysis folder and scripts

    <code>bash VODKA2_analysis_setup.sh [-h] -f <files.txt> -d <db_bt2_idx> -v <ref.fasta> -p <project> [options]</code>

* -f \<files.txt\>  TXT file containing the list of fastq files
* -d \<bt2_idx\>  Name of VODKA2 DB index for Bowtie2
* -v \<ref.fasta\> Virus reference genome in FASTA format
* -p \<project\>  Project name used to label output folders and files


 2. Run analysis

    <code>bash cbVG_VODKA2_analysis_\<project\>.sh</code><br/>
or<br/>
    <code>bash delVG_VODKA2_analysis_\<project\>.sh</code><br/>

<!-- TEST RUN -->
## TEST RUN

Download test data from this repository<br/>
* Reference genome:
    * RSVKC731482geneG.fasta
* Simuated Illumina Seq data:
    * RSV_geneG_artificial_CB_R1.fastq
    * RSV_geneG_artificial_CB_R1.fastq
<br/>

1. Start the VODKA2 Docker container
```
   docker pull emna20/vodka2:v1
```
Once the VODKA2 container has started, make sure Blast commands are in your PATH:
```
   export PATH="/blast/bin:$PATH"
``` 

2. Generate the VODKA2 cbVG database (DB)
```
   perl /VODKA2/genomefa_to_newfasta_cb_v2.pl RSVKC731482geneG.fasta 966 150 RSVKC731482geneG.966.150.fasta
```

3. Generate bowtie2 index for VODKA2 DB and ref genome
```
   bowtie2-build --large-index RSVKC731482geneG.966.150.fasta RSVKC731482geneG.966.150
```
```
   bowtie2-build RSVKC731482geneG.fasta RSVKC731482geneG.fasta
```

4. Setup the VODKA2 analysis script

Create a file with the fastq files:
```
   ls RSV_geneG_artificial_R*.fastq > samples.txt
```

Run the setup script:
```
   bash /VODKA2/VODKA2_analysis_setup.sh -f samples.txt -d RSVKC731482geneG.966.15 -v RSVKC731482geneG.fasta -p TEST
```

5. Run the analysis
```
   ./cbDVG_analysis_TEST.sh
```

6. Check results

    <b>TEST/ content should be:</b>
    
    <code>vodka2_output_CB/
    TEST_CB_dvg/</code>    (these folders contain intermediate outputs)<br/>
    <code>RSV_R1.vodka2.all-info_CB.N5.txt
    RSV_R2.vodka2.all-info_CB.N5.txt
    RSV_R1.vodka2.all-info_CB.N5_mode.txt
    RSV_R2.vodka2.all-info_CB.N5_mode.txt</code>    (these are the VODKA2 result tables that should be used for downstream analysis)<br/>
    <code>RSV_R1.vodka2.all-info_CB.N5_mode_plot.tiff
    RSV_R2.vodka2.all-info_CB.N5_mode_plot.tiff</code>     (plots using data from result tables)<br/>
    <code>TEST_CB_dvg_summary.txt</code>     (summary results table)<br/>
    
<br/>
<!-- CONTACT -->
## Contact

Emna Achouri - emna@wustl.edu

Carolina López - clopezzalaquett@wustl.edu

Project Link: [https://github.com/lopezlab-washu/VODKA2](https://github.com/lopezlab-washu/VODKA2)

<p align="right">(<a href="#top">back to top</a>)</p>

