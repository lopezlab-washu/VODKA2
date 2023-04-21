#!/bin/bash

set -e

usage=$'\n\t\033[1m*** delVG extra info report ***\033[0m\n
Extracts read orientation information from SAM file (FLAG).\n
Usage:\n\t\033[1m'$(basename "$0")$'\033[0m [-h] \033[1m-s SAMPLENAME -p PROJECT -t DEL|CB [-a GENOME_ANNOTATION -n NUMBER]\033[0m
where:
\t \033[1m-s  samplename \033[0m
\t \033[1m-p  project \033[0m
\t \033[1m-t  dvg type (DEL or CB)\033[0m
\t \033[1m-a  genome annotation file\033[0m (only for DEL)
\t \033[1m-n  shift nt within the same species group (default 5)
\t -h  show this help text and exit\n'

options=':hs:p:t:a:n:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    s) SAMPLENAME=$OPTARG;;
    p) PROJECT=$OPTARG;;
    a) ANNOT=$OPTARG;;
    t) DVGTYPE=$OPTARG;;
    n) N=$OPTARG;;
    :) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m missing argument for $OPTARG\n"; exit 1;;
   \?) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m illegal option -$OPTARG\n"; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$SAMPLENAME" ] || [ ! "$PROJECT" ] || [ ! "$DVGTYPE" ]; then
  echo "$usage"
  echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m options -s, -p and -t must be provided\n"
  exit 1
fi

if [ ! "$N" ]; then
  N=5
fi

# Extract genome length from vodka output file
# project folder
vodka_folder="${PWD}/${PROJECT}_vodka_results"
gl=$(ls ${vodka_folder}/vodka_output_${DVGTYPE}/*/*/*search_output.txt | head -1 | rev | cut -d "." -f 4 | rev)

echo
echo -e "\t* ${DVGTYPE}VG extra info report *"
echo

# extract DVG read ids from vodka output
vodka_results="${vodka_folder}/${PROJECT}_${DVGTYPE}_dvg/${SAMPLENAME}.2ranges.confirm.vodka.txt"
if [[ ! $(wc -l <${vodka_results}) -ge 2 ]]; then
  echo -e "Sample ${SAMPLENAME} does not have any ${DVGTYPE}VG"
  echo -e "Not making report. Done.\n"
  exit 0
fi
echo -e "Extracting DVG fastq read ids\n"
cut -f 2 ${vodka_results} | grep -v READ_ID > ${SAMPLENAME}.readids.txt

# extract FLAG info from SAM file
sam="${vodka_folder}/vodka_output_${DVGTYPE}/${SAMPLENAME}_vodka_output/alignment/${SAMPLENAME}.*.sam"
echo -e "Extracting read orientation\n"
join --nocheck-order <(sort ${SAMPLENAME}.readids.txt) <(sort ${sam}) | cut -f 1,2 -d " " > ${SAMPLENAME}.readstrands.txt

# write results
echo -e "Formatting output\n"
join -1 2 <(tail -n+2 ${vodka_results} | sort -k2) <(sort ${SAMPLENAME}.readstrands.txt) | awk -v OFS='\t' '{swap=$1;$1=$2;$2=swap;print $0}' > ${SAMPLENAME}.results_${DVGTYPE}.txt

# get gene positions
if [ $DVGTYPE == "DEL" ]; then
  if [ ! -z ${ANNOT} ]; then
    # get gene positions and dvg/del sizes
    echo -e "Getting the gene positions\n"
    declare -A dict
    while read i; do
      # dict key is region end position
      key=$(echo ${i} | cut -d " " -f 3 )
      # data is start position and gene name
      gene=$(echo ${i} | cut -d ' ' -f 1 )
      dict[$key]=${gene}
    done < ${ANNOT}

    # get genes list
    genes=$(grep -Ev 'name|leader|igr-|trailer' ${ANNOT} | cut -f 1)

    # store list of region 'end' positions
    keys=$(echo ${!dict[@]} | tr ' ' '\n' | sort -nr | grep -v 'end' | tr '\n' ' ' )

    # prepare output file (header line)
    echo -e "$(head -1 ${vodka_results})\tSAM_FLAG\tLENGTH\tBREAK_POSITION\tREJOIN_POSITION\tB_REGION\tB_GENE\tR_REGION\tR_GENE\tSAME_REGION\tDELETION_SIZE\tPERC_STANDARD\tSPECIES" > ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.N${N}.txt
  else
    echo -e "$(head -1 ${vodka_results})\tSAM_FLAG\tLENGTH\tBREAK_POSITION\tREJOIN_POSITION\tDELETION_SIZE\tPERC_STANDARD\tSPECIES" > ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.N${N}.txt
  fi

  # loop through dvg reads
  while read i; do
    # store Break and Rejoin positions
    B=$(echo ${i} | cut -f 1 -d "_" )
    R=$(echo ${i} | cut -f 2 -d "_" | cut -f 1 -d "(" )
    # add delVG length, deletion size and percentage of delVG compared to standard virus size
    DVGsize=$(( ${gl} - ${R} + ${B} + 1 ))
    delsize=$(( ${R} - ${B} - 1 ))
    perc=$(awk "BEGIN{print ${DVGsize}/${gl}*100}" )

    if [ ! -z ${ANNOT} ]; then
      # loop through genes 'end' position
      for k in ${keys}; do
        # check B is <= current end
        if (( ${B} <= ${k} )); then GB=${dict[$k]}; fi
      done
      for k in ${keys}; do
        # check R is <= current end
        if (( ${R} <= ${k} )); then GR=${dict[$k]}; fi
      done
      # check whether B and R are in the same gene
      if [ ${GB} == ${GR} ]; then same="YES"; else same="NO"; fi
      # check B/R within a gene
      if echo ${genes} | grep -w -q ${GB} ; then Bgene="YES"; else Bgene="NO"; fi
      if echo ${genes} | grep -w -q ${GR} ; then Rgene="YES"; else Rgene="NO"; fi
      # write to file
      echo -e "${i}\t${DVGsize}\t${B}\t${R}\t${GB}\t${Bgene}\t${GR}\t${Rgene}\t${same}\t${delsize}\t${perc}" >> ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp
    else
      echo -e "${i}\t${DVGsize}\t${B}\t${R}\t${delsize}\t${perc}" >> ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp
    fi
    done < ${SAMPLENAME}.results_${DVGTYPE}.txt

elif [ $DVGTYPE == "CB" ]; then
  # prepare output file (header line)
  echo -e "$(head -1 ${vodka_results})\tSAM_FLAG\tLENGTH\tBREAK_POSITION\tREJOIN_POSITION\tLOOP_SIZE\tSTEM_SIZE\tPERC_STEM\tSPECIES" > ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.N${N}.txt
  # get dvg/stem/loop sizes
  while read i ; do
    # store Break and Rejoin positions
    B=$(echo ${i} | cut -f 1 -d "_" )
    R=$(echo ${i} | cut -f 2 -d "_" | cut -f 1 -d "(" )
    # delVG size is (gl - Break + Rejoin +1)
    DVGsize=$(( ${gl} - ${B} + 1 + ${gl} - ${R} + 1 ))
    loopsize=$(( ${R} - ${B} ))
    stemsize=$(( $(( ${gl} - ${R} + 1 )) * 2 ))
    perstem=$(awk "BEGIN{print ${stemsize}/${DVGsize}*100}" )
    echo -e "${i}\t${DVGsize}\t${B}\t${R}\t${loopsize}\t${stemsize}\t${perstem}" >> ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp
  done < ${SAMPLENAME}.results_${DVGTYPE}.txt
fi

# get species occurence/break/rejoin/size/ (species is defined by SIZE and B/R) and sort by SIZE and B
sort -k13n -k14n ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp > ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp.sorted
# initialize variables: occurence, break and rejoin
# flag variable is used for flagging equivalent B and R in order to aggregate them
# o1/O: junction occurence, b1/B/r1/R: BREAK and REJOIN positions, s1/S: DVG size
b1=0; r1=0; s1=0; species=(); flag=0
while read i; do
  S=$(echo ${i} | cut -d ' ' -f 13)
  B=$(echo ${i} | cut -d ' ' -f 14)
  R=$(echo ${i} | cut -d ' ' -f 15)
  if [[ ${s1} == ${S} ]]; then # same size
    # start new species group
    if [[ ${flag} == 0 ]]; then # prev species is new
      if (( ${B} <= $(( ${b1} + ${N} )) )); then # current species is similar to prev
        # store new B/R/S and species list
        b1=${B}; r1=${R}; s1=${S}; species+=("${i}")
        # flag updated species
        flag=1
      else # current species is different
        # print prev species info to file
        OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
          echo -e "${j}\t${b1}_${r1}"
        done; IFS=$OLDIFS
        # store current species info
        b1=${B}; r1=${R}; s1=${S}; species=("${i}")
      fi
    # update previous species group
    elif (( ${B} <= $(( ${b1} + ${N} )) )); then # current species is similar to prev
      # store new B/R/S and sum of occurences
      b1=${B}; r1=${R}; s1=${S}; species+=("${i}")
    else # current species is different
      # print prev species info to file
      OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
        echo -e "${j}\t${b1}_${r1}"
      done; IFS=$OLDIFS
      # store current species info
      b1=${B}; r1=${R}; s1=${S}; species=("${i}")
      # flag new species
      flag=0
    fi
  else # current species is different
    #print prev species info to file
    OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
      echo -e "${j}\t${b1}_${r1}"
    done; IFS=$OLDIFS   
    # store current species info
    b1=${B}; r1=${R}; s1=${S}; species=("${i}")
    flag=0
  fi
done < ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp.sorted >> ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.N${N}.txt

# Don't forget the last species (not printed at the end of the loop)
lastrow=$(tail -n 1 ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp.sorted)
B=$(echo ${lastrow} | cut -d ' ' -f 14)
R=$(echo ${lastrow} | cut -d ' ' -f 15)
# print to file
OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
  echo -e "${j}\t${B}_${R}" >> ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.N${N}.txt
done; IFS=$OLDIFS

rm ${SAMPLENAME}.readids.txt ${SAMPLENAME}.readstrands.txt ${SAMPLENAME}.results_${DVGTYPE}.txt ${SAMPLENAME}.vodka2.all-info_${DVGTYPE}.tmp*
echo -e "Done.\n"
