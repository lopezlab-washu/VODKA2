#!/bin/bash

# Make Blast DB previous to running this sript
# makeblastdb -in ref.fasta -dbtype nucl -parse_seqids

set -e

usage=$'\n\t\033[1m*** BLAST VALIDATION SCRIPT v4 ***\033[0m\n
Runs blast to validate the DVG junction regions predicted by VODKA.\n
Usage:\n\t\033[1m'$(basename "$0")$'\033[0m [-h] \033[1m-r REF_GENOME -t DVG_TYPE -v VODKA_RESULTS -s SUFFIX\033[0m [-e EVALUE] [-n NUMBER]
where:
\t \033[1m-r  reference genome for BLAST\033[0m
\t \033[1m-t  type of DVG\033[0m (CB or DEL)
\t \033[1m-v  VODKA results table\033[0m (*_RESULTS.txt)
\t \033[1m-s  suffix of VODKA output file\033[0m (anything after samplename)
(option) -e  E-value threshold (for BLAST, default: 0.001)
(option) -n  add this number to BREAK and REJOIN position (default: 5)
\t -h  show this help text and exit\n'

DVGTYPE="CB"
EVALUE=0.001
N=5

options=':hr:v:s:e:n:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    r) REFDB=$OPTARG;;
    v) SAMPLETXT=$OPTARG;;
    s) SUFFIX=$OPTARG;;
    e) EVALUE=$OPTARG;;
    n) N=$OPTARG;;
    t) DVGTYPE=$OPTARG;;
    :) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m missing argument for $OPTARG\n"; exit 1;;
   \?) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m illegal option -$OPTARG\n"; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$REFDB" ] || [ ! "$DVGTYPE" ] || [ ! "$SAMPLETXT" ] || [ ! "$SUFFIX" ]; then
  echo "$usage"
  echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m options -r, -t, -v and -s must be provided\n"
  exit 1
fi

echo
echo -e "\t* VODKA ${DVGTYPE} DVG validation analysis *"
echo

# extract DVG junction region sequences from vodka output
echo -e "Extracting DVG fasta sequences\n"
samplefaB=$(basename ${SAMPLETXT} .txt)_B.fa
samplefaR=$(basename ${SAMPLETXT} .txt)_R.fa
sample=$(basename ${SAMPLETXT} ${SUFFIX})
tail -n +2 ${SAMPLETXT} | cut -f 1,2,10 | tr -d "<>" | sed 's/[ATCGN]*::://' | sed 's/\t/-/' | sed 's/^/>/g' | tr "\t" "\n" > ${fastaB}
tail -n +2 ${SAMPLETXT} | cut -f 1,2,10 | tr -d "<>" | sed 's/:::.*//' | sed 's/\t/-/' | sed 's/^/>/g' | tr "\t" "\n" > ${fastaR}

# genome length
gl=$(echo ${SUFFIX} | rev | cut -d "." -f 3 | rev)

# run BLAST
echo -e "Blasting against ref ${REFDB}"
echo "E-Value threshold is ${EVALUE}"
# upstream B
blastn -max_hsps 1 -db ${REFDB} -query ${samplefaB} -out ${sample}_B.blast -outfmt "6 qseqid qstart qend sseqid sstart send sstrand" -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -evalue ${EVALUE} -perc_identity 0.1 || exit 1
# downstream R
blastn -max_hsps 1 -db ${REFDB} -query ${samplefaR} -out ${sample}_R.blast -outfmt "6 qseqid qstart qend sseqid sstart send sstrand" -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -evalue ${EVALUE} -perc_identity 0.1 || exit 1
cat ${sample}_B.blast ${sample}_R.blast > ${sample}.blast
echo "Done!"
echo

# set n variable according to DVG type (used for strand orientation checking)
if [ ${DVGTYPE} = "DEL" ]; then n=1;
elif [ ${DVGTYPE} = "CB" ]; then n=2;
fi

echo -e "Checking BLAST results..."
# sort blast output
sort ${sample}.blast > ${sample}.blast.sorted
# count number of ranges (e.g. DVG with 3 ranges " 3 B_R(B_R)")
cut -f 1 ${sample}.blast.sorted | uniq -c > ${sample}.blast.counts
# store list of ranges values (1,2,3,4,17)
if [[ ! $(less ${sample}.blast.counts | grep -Eoh "^ +[0-9]+ " | tr -d " " | sort -n | uniq | grep -v ^1$) ]]; then
  # print report
  n=$(( $(wc -l ${SAMPLETXT} | cut -d " " -f 1) - 1 ))
  echo " * Found 0 ${DVGTYPE} DVG \"junctions\" (supported by 0 DVG reads) *"
  echo "VODKA detected ${n} DVG reads."
  echo "Blast reported 2 ranges for 0 DVG."

  # delete intermediate files
  rm ${sample}.blast.*
  rm ${samplefa} ${SAMPLETXT}

  echo
  echo "Finished."
  echo

  exit 0
fi
less ${sample}.blast.counts | grep -Eoh "^ +[0-9]+ " | tr -d " " | sort -n | uniq | grep -v ^1$ > ${sample}.blast.uniq.ranges
# keep only DVG with at least 2 ranges
less ${sample}.blast.counts | sed  "s/^ *//" | sort -n -k 1 | grep -v "^1 " > ${sample}.blast.junctions
# split DVG list according to number of ranges (e.g. DVG with 3 ranges is saved into file sample.blast.junctions.3)
while read i ; do grep -E "^${i} " ${sample}.blast.junctions | cut -d " " -f 2 > ${sample}.blast.junctions.${i} ; done < ${sample}.blast.uniq.ranges 
# set n variable according to DVG type (used for strand orientation checking)
if [ ${DVGTYPE} = "CB" ]; then n=1; elif [ ${DVGTYPE} = "DEL" ]; then n=2; fi
# keep only DVG with 2 ranges
sort ${sample}.blast.junctions.2 > ${sample}.blast.junctions.2.sorted
join -j1 -o "0,1.7" ${sample}.blast.sorted ${sample}.blast.junctions.2.sorted > ${sample}.blast.junctions.2.strands
# check ref strand orientation (sstrand values for range1,range2 are plus,minus or minus,plus for CB (n=2) and plus,plus or minus,minus for DEL (n=1))
uniq -c ${sample}.blast.junctions.2.strands | sed  "s/^ *//" | grep "^${n} " | cut -d " " -f 2 | uniq > ${sample}.blast.junctions.2.checked
# extract reads with the 2 ranges and correct strand orientation from blast output
join -j1 ${sample}.blast.sorted ${sample}.blast.junctions.2.checked > ${sample}.blast.junctions.2.data

# Keep only DVG with Break and Rejoin position included within the alignment positions for the 2 ranges
echo -e "Getting blast B/R positions"
echo -e "N value is ${N}"
# DELETION DVGs
if  [ ${DVGTYPE} = "DEL" ]; then
  i=0
  while read li ; do
    # BREAK +/- N
    BREAK=$(echo ${li} | cut -d "_" -f 1)
    BREAKminusN=$(( ${BREAK} - ${N} ));
    BREAKplusN=$(( ${BREAK} + ${N} ));
    # REJOIN +/- N
    REJOIN=$(echo ${li} | cut -d "_" -f 2 | cut -d "(" -f 1)
    REJOINminusN=$(( ${REJOIN} - ${N} ));
    REJOINplusN=$(( ${REJOIN} + ${N} ));
    if [[ $(( $i % 2 )) == 0 ]]; then
      # range 1
      range1start=$(echo ${li} | cut -d " " -f 5);
      range1end=$(echo ${li} | cut -d " " -f 6);
    else
      # range 2
      range2start=$(echo ${li} | cut -d " " -f 5);
      range2end=$(echo ${li} | cut -d " " -f 6);
      if (( ${range1end} < ${range2start} )); then
        # if B is included in range 1 and R is in range 2
        if (( ${BREAKminusN} <= ${range1end} && ${range1end} <= ${BREAKplusN} && ${REJOINminusN} <= ${range2start} && ${range2start} <= ${REJOINplusN} )) ; then
          # print to file vodka B_R 
          echo ${li} | cut -d " " -f 1
        fi
      elif (( ${range2end} < ${range1start} )); then
        # if B is included in range 2 and R is in range 1
        if (( ${BREAKminusN} <= ${range2end} && ${range2end} <= ${BREAKplusN} && ${REJOINminusN} <= ${range1start} && ${range1start} <= ${REJOINplusN} )) ; then
          # print to file vodka B_R
          echo ${li} | cut -d " " -f 1
        fi
      fi
    fi
    i=$((i+1))
  done < ${sample}.blast.junctions.2.data > ${sample}.blast.junctions.2.confirmed

# COPYBACK DVGs
elif  [ ${DVGTYPE} = "CB" ]; then
  i=0
  while read li ; do
    # BREAK +/- N
    BREAK=$(echo ${li} | cut -d "_" -f 1)
    BREAKminusN=$(( ${BREAK} - ${N} ));
    BREAKplusN=$(( ${BREAK} + ${N} ));
    # REJOIN +/- N
    REJOIN=$(echo ${li} | cut -d "_" -f 2 | cut -d "(" -f 1)
    REJOINminusN=$(( ${REJOIN} - ${N} ));
    REJOINplusN=$(( ${REJOIN} + ${N} ));
    if [[ $(( $i % 2 )) == 0 ]]; then
      # range 1
      r1s=$(echo ${li} | cut -d " " -f 5);
      r1e=$(echo ${li} | cut -d " " -f 6);
    else
      # range 2
      r2s=$(echo ${li} | cut -d " " -f 5);
      r2e=$(echo ${li} | cut -d " " -f 6);
      range1start=$(( ${r1s} < ${r1e} ? ${r1s} : ${r1e} ));
      range2start=$(( ${r2s} < ${r2e} ? ${r2s} : ${r2e} ));
      if (( ${range1start} < ${range2start} )); then
        # if B is included in range 1 and R is in range 2
        if (( ${BREAKminusN} <= ${range1start} && ${range1start} <= ${BREAKplusN} && ${REJOINminusN} <= ${range2start} && ${range2start} <= ${REJOINplusN} )) ; then
          # print to file vodka B_R 
          echo ${li} | cut -d " " -f 1
        fi
      elif (( ${range2start} < ${range1start} )); then
        # if B is included in range 2 and R is in range 1
        if (( ${BREAKminusN} <= ${range2start} && ${range2start} <= ${BREAKplusN} && ${REJOINminusN} <= ${range1start} && ${range1start} <= ${REJOINplusN} )) ; then
          # print to file vodka B_R
          echo ${li} | cut -d " " -f 1
        fi
      fi
    fi
    i=$((i+1))
  done < ${sample}.blast.junctions.2.data > ${sample}.blast.junctions.2.confirmed
fi

echo
echo -e "Formatting results"
# extract and report confirmed DVG info from vodka and blast outputs
head -1 $(basename ${SAMPLETXT}) > ${sample}.2ranges.confirm.vodka.txt
while read i ; do dvg=$(echo ${i} | tr "-" "\t"); grep "^${dvg}" $(basename ${SAMPLETXT}) ; done < ${sample}.blast.junctions.2.confirmed >> ${sample}.2ranges.confirm.vodka.txt
echo -e "6\tqseqid\tqstart\tqend\tsseqid\tsstart\tsend\tsstrand" > ${sample}.2ranges.confirm.blast.txt
while read i ; do grep "^${i}" ${sample}.blast ; done < ${sample}.blast.junctions.2.confirmed >> ${sample}.2ranges.confirm.blast.txt
# report the list of confirmed DVG junctions
echo "A_C(A_C)" > ${sample}.2ranges.confirm.list.txt
cut -d "-" -f 1 ${sample}.blast.junctions.2.confirmed | uniq >> ${sample}.2ranges.confirm.list.txt

# print report
n=$(( $(wc -l ${SAMPLETXT} | cut -d " " -f 1) - 1 ))
r2=$(wc -l ${sample}.blast.junctions.2 | cut -d " " -f 1)
vodkajunc=$(( $(wc -l ${sample}.2ranges.confirm.list.txt | cut -d " " -f 1) - 1 ))
reads=$(( $(wc -l ${sample}.2ranges.confirm.vodka.txt | cut -d " "  -f 1) - 1 ))
echo " * Found ${vodkajunc} ${DVGTYPE} DVG \"junctions\" (supported by ${reads} DVG reads) *"
echo "VODKA detected ${n} DVG reads."
echo "Blast reported 2 ranges for ${r2} DVG."

# delete intermediate files
rm ${sample}.blast.*
rm ${samplefa} ${SAMPLETXT}

echo
echo "Finished."
echo
