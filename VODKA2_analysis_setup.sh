#!/bin/bash

# make sure to add VODKA2 scripts folder to $PATH

set -e

usage=$'\n\t\033[1m*** DVG ANALYSIS SETUP ***\033[0m\n
Setup VODKA and BLAST anaysis.\n
Usage:\n\t\033[1m'$(basename "$0")$'\033[0m [-h] \033[1m-f SAMPLES.txt -d VODKA_DVG_index -v REF_VIRUS -p PROJECT\033[0m [-n NUMBER]
where:
\t \033[1m-f <file.txt>\t\tTXT file containing the list of fastq files\033[0m (1 per sample, with the relative path)
\t \033[1m-d <vodka_index>\tName of VODKA DB index for Bowtie2\033[0m (including path to folder)
\t \033[1m-v <virus_fasta>\tVirus reference genome in fasta format\033[0m (including path to folder)
\t \033[1m-p <project>\t\tProject name to label files\033[0m
(option) -n <number>\t\tnumber of nt shift allowed for Blast B/R (default: 5) 
\t -h  \t\t\tshow this help text and exit\n'

# default values
number=5

options=':hf:d:v:p:g:r:n:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    f) samplesfile=$OPTARG;;
    v) virusfas=$OPTARG;;
    d) vodkabt2=$OPTARG;;
    p) project=$OPTARG;;
    n) number=$OPTARG;;
    :) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m missing argument for $OPTARG\n"; exit 1;;
   \?) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m illegal option -$OPTARG\n"; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$samplesfile" ] || [ ! "$vodkabt2" ] || [ ! "$project" ] ; then
  echo "$usage"
  echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m options -f, -d, and -p must be provided\n"
  exit 1
fi

# project folder
vodka_folder=$(realpath ${project}_vodka_results)

# virus ref and vodka2 db
vodkabt2_name=$(basename ${vodkabt2})
virusfas_name=$(basename ${virusfas})

# suffix of VODKA result filenames
vodka_suffix="${vodkabt2_name}_RESULTS.txt"

# extract DVGtype, genome length and read length from VODKA index name
gl=$(echo ${vodkabt2_name} | rev | cut -f 2 -d "." | rev)
rl=$(echo ${vodkabt2_name} | rev | cut -f 1 -d "." | rev)
ty=$(echo ${vodkabt2_name} | rev | cut -f 3 -d "." | rev)

# setup variables and script name according to DVG type
if [ ${ty} == "DEL" ]; then
  analysis_script=delDVG_analysis_${project}.sh
  DVGtype=DEL
else
  analysis_script=cbDVG_analysis_${project}.sh
  DVGtype=CB
fi

echo
echo -e "\033[1m\t* Setup ${DVGtype} DVG analysis *\033[0m"
echo


# print the list of samples
echo -e "\033[4mList of samples:\033[0m"
echo
while read i; do
  samplename=$(basename "${i}" | cut -d "_" -f 1)_$(basename "${i}" | sed -e "s/.*\(R[12]\).*/\1/")
  echo ${samplename}
done < ${samplesfile}

# check sample list is ok
echo
echo -n "Please confirm the samplenames listed above are correct "
read -p "(y/n): "
if [[ $REPLY =~ ^[Yy]$ ]]; then
  echo "thank you."
else
  echo
  echo -e "\033[1mPlease check samplename parsing method!\033[0m Exit."
  echo
  exit 1
fi

# create results folder on Acive
mkdir -p ${vodka_folder}/vodka_output_${DVGtype} ${vodka_folder}/${project}_${DVGtype}_dvg
cp ${virusfas} ${vodka_folder}/${project}_${DVGtype}_dvg/.
if [ ${ty} == "DEL" ]; then
  cp *annotation.txt ${vodka_folder}/.
fi

## SETUP BLAST AND SPECIES ANALYSIS SCRIPTS ###
echo "cd ${vodka_folder}/${project}_${DVGtype}_dvg" > vodka_blast_${project}_${DVGtype}.sh
cmd="makeblastdb -in \"${virusfas}\" -out ${virusfas_name} -dbtype nucl"
echo "${cmd}" >> vodka_blast_${project}_${DVGtype}.sh
echo "cd ${vodka_folder}" > vodka_species_${project}_${DVGtype}.sh

echo
echo -e "\033[4mResults folder:\033[0m ${vodka_folder}"
echo
echo "(content)"
ls -1 ${vodka_folder}
echo

## VODKA ANALYSIS SETUP ##
while read i; do
        samplename=$(basename "${i}" | cut -d "_" -f 1)_$(basename ${i} | sed -e "s/.*\(R[12]\).*/\1/")

	# create sample analysis folders
	alignment="${vodka_folder}/vodka_output_${DVGtype}/${samplename}_vodka_output/alignment"
	results="${vodka_folder}/vodka_output_${DVGtype}/${samplename}_vodka_output/results"
	mkdir -p ${alignment} ${results}

	# VODKA STEP 1: Alignment against VODKA DVG database
	cmd1="bowtie2 -p 30 -x \"${vodkabt2}\" --local -U \"${i}\" -S ${alignment}/${samplename}.${vodkabt2_name}.sam --no-sq --no-unal --mp 0,0"
	echo "${cmd1}" >> vodka_step1_${project}_${DVGtype}.sh

	# VODKA STEP 2: Extract reads matching DVG ref
	cmd2="perl search1.pl ${alignment}/${samplename}.${vodkabt2_name}.sam ${alignment}/${samplename}.${vodkabt2_name}.search_output.txt ${gl} ${rl}"
	echo "${cmd2}" >> vodka_step2_${project}_${DVGtype}.sh

	# VODKA STEP 3: generate DVG candidates fasta file (for step4) and prepare results files (for step5)
	cmd3="perl organize.pl ${vodka_folder}/vodka_output_${DVGtype}/${samplename}_vodka_output ${vodkabt2_name} ${samplename}"
	echo "${cmd3}" >> vodka_step3_${project}_${DVGtype}.sh

	## BLAST VERIFICATION ##
	echo "cp ${results}/${vodka_suffix} ${samplename}_${vodka_suffix}" >> vodka_blast_${project}_${DVGtype}.sh
	echo "bash VODKA2_blast.sh -r ${virusfas_name} -v ${samplename}_${vodka_suffix} -s _${vodka_suffix} -t ${DVGtype}" >> vodka_blast_${project}_${DVGtype}.sh

	## EXTENDED REPORT WITH SPECIES ##
        if [ ${DVGtype} == "DEL" ]; then
          echo "bash VODKA2_report.sh -s ${samplename} -p ${project} -t ${DVGtype} -n ${number} -a *annotation.txt" >> vodka_species_${project}_${DVGtype}.sh
        else
	  echo "bash VODKA2_report.sh -s ${samplename} -p ${project} -t ${DVGtype} -n ${number}" >> vodka_species_${project}_${DVGtype}.sh
        fi

done < ${samplesfile}

### SETUP SUMMARY REPORT SCRIPT ###
echo "cd ${vodka_folder}" > vodka_summary_${project}_${DVGtype}.sh
echo "echo -e \"sample\tnb_junctions\tnb_species\tnb_${DVGtype}_reads\" > ${project}_${DVGtype}_dvg_summary.txt" >> vodka_summary_${project}_${DVGtype}.sh
while read i; do
  samplename=$(basename "${i}" | cut -d "_" -f 1)_$(basename "${i}" | sed -e "s/.*\(R[12]\).*/\1/")
  if [ ${DVGtype} == "DEL" ]; then x=23; else x=19; fi
  samplereport=${vodka_folder}/${samplename}.vodka2.all-info_${DVGtype}.N${number}.txt
  echo "if [ -f \"${samplereport}\" ]; then
    nb_junc=\$( tail -n +2 "${samplereport}" | cut -f 1 | sort | uniq | wc -l )
    nb_spec=\$( tail -n +2 "${samplereport}" | cut -f ${x} | sort | uniq | wc -l )
    nb_read=\$( tail -n +2 "${samplereport}" | wc -l )
  else
    nb_junc=0
    nb_spec=0
    nb_read=0
  fi
  echo -e \"${samplename}\t\${nb_junc}\t\${nb_spec}\t\${nb_read}\" >> ${project}_${DVGtype}_dvg_summary.txt" >> vodka_summary_${project}_${DVGtype}.sh
done < ${samplesfile} >> vodka_summary_${project}_${DVGtype}.sh

cat vodka_step*_${project}_${DVGtype}.sh vodka_blast_${project}_${DVGtype}.sh vodka_species_${project}_${DVGtype}.sh vodka_summary_${project}_${DVGtype}.sh >> ${analysis_script}
chmod 777 ${analysis_script}

echo
echo "Your script is ready!"
echo
echo -e "Run analysis: \033[1m./${analysis_script}\033[0m"
echo

