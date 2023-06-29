#!/bin/bash

# make sure to add VODKA2 scripts folder to $PATH

set -e

usage=$'\n\t\033[1m*** VODKA2 nsVG ANALYSIS SETUP ***\033[0m\n
Usage:\n\t\033[1m'$(basename "$0")$'\033[0m [-h] \033[1m-f SAMPLES.txt -d VODKA2_DB -v REF_VIRUS -p PROJECT\033[0m [-n NUMBER]
where:
\t \033[1m-f <file.txt>\t\tTXT file containing the list of fastq files\033[0m (1 per sample, with the relative path)
\t \033[1m-d <vodka2_index>\tName of VODKA2 DB index for Bowtie2\033[0m (including path to folder)
\t \033[1m-v <virus_fasta>\tVirus reference genome in fasta format\033[0m (including path to folder)
\t \033[1m-p <project>\t\tProject name to label files\033[0m
(option) -n <number>\t\tnumber of nt shift allowed for species aggregation (default: 5) 
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

# VODKA2 scripts folder
VODKA2scripts="/VODKA2"

# project folder
vodka2_folder="${PWD}/${project}_vodka2_results"

# virus ref and vodka2 db
vodkabt2_name=$(basename ${vodkabt2})
virusfas_name=$(basename ${virusfas})

# suffix of VODKA2 result filenames
vodka2_suffix="${vodkabt2_name}_RESULTS.txt"

# extract DVGtype, genome length and read length from VODKA2 index name
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

# create results folder
mkdir -p ${vodka2_folder}/vodka2_output_${DVGtype} ${vodka2_folder}/${project}_${DVGtype}_dvg
cp ${virusfas} ${vodka2_folder}/${project}_${DVGtype}_dvg/.
if [ ${ty} == "DEL" ]; then
  cp *annotation.txt ${vodka2_folder}/.
fi

## SETUP FOLDERS AND REF INDEX ##
cmd="bowtie2-build ${virusfas_name} ${virusfas_name}"
echo "${cmd}" >> vodka2_step1_${project}_${DVGtype}.sh
echo "cd ${vodka2_folder}/${project}_${DVGtype}_dvg" > vodka2_step5_${project}_${DVGtype}.sh
cmd="makeblastdb -in ${virusfas_name} -out ${virusfas_name} -dbtype nucl"
echo "${cmd}" >> vodka2_step5_${project}_${DVGtype}.sh
echo "cd ${vodka2_folder}" > vodka2_step6_${project}_${DVGtype}.sh
echo "cd ${vodka2_folder}" > vodka2_step7_${project}_${DVGtype}.sh

echo
echo -e "\033[4mResults folder:\033[0m ${vodka2_folder}"
echo
echo "(content)"
ls -1 ${vodka2_folder}
echo

## VODKA2 ANALYSIS SETUP ##
while read i; do
        samplename=$(basename "${i}" | cut -d "_" -f 1)_$(basename ${i} | sed -e "s/.*\(R[12]\).*/\1/")

	# create sample analysis folders
	alignment="${vodka2_folder}/vodka2_output_${DVGtype}/${samplename}_vodka2_output/alignment"
	results="${vodka2_folder}/vodka2_output_${DVGtype}/${samplename}_vodka2_output/results"
        nonVirus="${vodka2_folder}/vodka2_output_${DVGtype}/${samplename}_vodka2_output/nonVirus"
	mkdir -p ${alignment} ${results} ${nonVirus}

        # STEP 1 : Host reads removal
        cmd1="bowtie2 -p 12 -x ${virusfas_name} -U ${i} --un-gz ${nonVirus}/${samplename}_nonViral.fq.gz -S ${samplename}_bowtie2virus.sam --no-sq --no-unal"
        echo "${cmd1}" >> vodka2_step1_${project}_${DVGtype}.sh

	# STEP 2: Alignment against VODKA2 DVG database
	cmd2="bowtie2 -p 12 -x \"${vodkabt2}\" --local -U ${nonVirus}/${samplename}_nonViral.fq.gz -S ${alignment}/${samplename}.${vodkabt2_name}.sam --no-sq --no-unal --mp 0,0"
	echo "${cmd2}" >> vodka2_step2_${project}_${DVGtype}.sh

	# STEP 3: Extract reads matching DVG ref
	cmd3="perl ${VODKA2scripts}/search1.pl ${alignment}/${samplename}.${vodkabt2_name}.sam ${alignment}/${samplename}.${vodkabt2_name}.search_output.txt ${gl} ${rl}"
	echo "${cmd3}" >> vodka2_step3_${project}_${DVGtype}.sh

	# STEP 4: Generate DVG candidates fasta file (for step5) and prepare results files (for step6)
	cmd4="perl ${VODKA2scripts}/organize.pl ${vodka2_folder}/vodka2_output_${DVGtype}/${samplename}_vodka2_output ${vodkabt2_name} ${samplename}"
	echo "${cmd4}" >> vodka2_step4_${project}_${DVGtype}.sh

	# STEP 5: BLAST VERIFICATION
	cmd5="cp ${results}/${vodka2_suffix} ${samplename}_${vodka2_suffix} && \
        bash ${VODKA2scripts}/VODKA2_blast.sh -r ${virusfas_name} -v ${samplename}_${vodka2_suffix} -s _${vodka2_suffix} -t ${DVGtype}"
        echo "${cmd5}" >> vodka2_step5_${project}_${DVGtype}.sh

	# STEP 6: Extended report
        if [ ${DVGtype} == "DEL" ]; then
          cmd6="bash ${VODKA2scripts}/VODKA2_report.sh -s ${samplename} -p ${project} -t ${DVGtype} -n ${number} -a *annotation.txt"
        else
	  cmd6="bash ${VODKA2scripts}/VODKA2_report.sh -s ${samplename} -p ${project} -t ${DVGtype} -n ${number}"
        fi
        echo ${cmd6} >> vodka2_step6_${project}_${DVGtype}.sh

        # STEP 7: Species ID (mode) and plot
        report=${samplename}.vodka2.all-info_${DVGtype}.N${number}.txt
        cmd7="Rscript ${VODKA2scripts}/VODKA2_species_plot.R ${report} ${gl} ${samplename} ${DVGtype}"
        echo "${cmd7}" >> vodka2_step7_${project}_${DVGtype}.sh 

done < ${samplesfile}

### SETUP SUMMARY REPORT SCRIPT ###
echo "cd ${vodka2_folder}" > vodka2_summary_${project}_${DVGtype}.sh
echo "echo -e \"sample\tnb_junctions\tnb_species\tnb_${DVGtype}_reads\" > ${project}_${DVGtype}_dvg_summary.txt" >> vodka2_summary_${project}_${DVGtype}.sh
while read i; do
  samplename=$(basename "${i}" | cut -d "_" -f 1)_$(basename "${i}" | sed -e "s/.*\(R[12]\).*/\1/")
  if [ ${DVGtype} == "DEL" ]; then x=19; else x=20; fi
  samplereport=${vodka2_folder}/${samplename}.vodka2.all-info_${DVGtype}.N${number}_mode.txt
  echo "if [ -f \"${samplereport}\" ]; then
    nb_junc=\$( tail -n +2 "${samplereport}" | cut -f 1 | sort | uniq | wc -l )
    nb_spec=\$( tail -n +2 "${samplereport}" | cut -f ${x} | sort | uniq | wc -l )
    nb_read=\$( tail -n +2 "${samplereport}" | wc -l )
  else
    nb_junc=0
    nb_spec=0
    nb_read=0
  fi
  echo -e \"${samplename}\t\${nb_junc}\t\${nb_spec}\t\${nb_read}\" >> ${project}_${DVGtype}_dvg_summary.txt" >> vodka2_summary_${project}_${DVGtype}.sh
done < ${samplesfile} >> vodka2_summary_${project}_${DVGtype}.sh

cat vodka2_step*_${project}_${DVGtype}.sh vodka2_summary_${project}_${DVGtype}.sh >> ${analysis_script}
chmod 777 ${analysis_script}

echo
echo "Your script is ready!"
echo
echo -e "Run analysis: \033[1m./${analysis_script}\033[0m"
echo

