#!/bin/bash

version="version : HoLDuP 0.05"

usage() {
	
	
	#error code
	#if "-h" option is set, status of usage =0, else it's 1 (=error)
	status=$1
	
	echo -e "$version\n\nUsage: $0 -f <sample design> -c <programs config>\n
 \tRequired arguments :\n                 
                    -f <sample design (you can adapt the file \"sample_design.txt\" supplied with the pipeline)>\n
                        \t\t-> example of supplied file (1st column fastq file, 2nd the condition, 3rd mate indication, 4th the library type) :\n
                      
                     \n\t\t\t\t/home/path/to/B67T11.R1.fastq.gz	prostate_tumor	mate1	fr-firststrand
                     \n\t\t\t\t/home/path/to/B67T11.R2.fastq.gz	prostate_tumor	mate2	fr-firststrand
                     \n\t\t\t\t/home/path/to/B67T14.R1.fastq.gz	prostate_normal	mate1	fr-firststrand
                     \n\t\t\t\t/home/path/to/B67T14.R2.fastq.gz	prostate_normal	mate2	fr-firststrand\n
                   
                                            \t********* warnings *********
                   
                              \t\t-> keep the order of mate1 & mate2 (don't put mate2 before mate1)
                               \t\t-> between each column, we have tabulations (not spaces)
                   
                                          \t***************************\n
                     \n\t\t  =========================================================================================\n 
                                                         
                    -c <programs config (adapt the \"programs_config.txt\" supplied with the pipeline)>\n
                    \t-> text file containing the location of the programs (or their bin)\n
                    \n\t\t  =========================================================================================\n                                  
 \toptional arguments :\n
 
 
                    -h < show Usage >\n                
                    
                    -v < show version >\n
                    \n\t\t  =========================================================================================\n                                           
\tResults :\n

                  -  4 gtf files of intergenic & antisense contigs (class 1 & 3)
                  
                      \t-> contigs_sans_sens.anti.class3w.quantile.0.2.gtf
                      \t-> contigs_sans_sens.inter.class3w.quantile.0.2.gtf
                      \t-> contigs_sans_sens.anti.class1w.quantile.0.2.gtf
                      \t-> contigs_sans_sens.inter.class1w.quantile.0.2.gtf
                      
                      remarks : 
                      
                        * class 3 contigs : contigs with RPKM expression above the value of a quantile threshold (we use 0.2)
                      
                        * class 1 contigs : subset of class 3 ; contigs with RPKM expression above the value of a quantile threshold + presence of EST(s) + presence of junction(s)          
\n" 1>&2; exit $status;}
                  
[[ $# -eq 0 ]] && usage 1

while getopts ":f:c:hv" opt; do
  case $opt in
  

      f)
             #design of the samples
             sample_design=$OPTARG
      
      
             ;;
      
      
      c)
             #take the config file
             programs_config=$OPTARG
     
             ;;
             
      h)
             #give the Usage
             usage 0
             
             
            ;;             
             
	  v)
             #give the version
			 echo -e "$version"
			 
			 exit 0

			 ;;
      

      
      #invalid options (options not in the list)
      ######################
      
      
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    
  esac
done

if [ ! -f "$sample_design" ] || [ ! -f "$programs_config" ]; then
      
     echo -e "\nplease, supply the config files ! !!\n"
     
     usage 1
      
fi

echo -e "\nthe design of the samples is in : $sample_design\n"


########## parse the programs config file #########

#retrieve samtools
process_samtools=$(grep "samtools=" $programs_config|grep -v "^#")
eval $process_samtools

#retrieve bedtools path
process_bedtools_path=$(grep "BEDTOOLS_PATH=" $programs_config|grep -v "^#")
eval $process_bedtools_path

#retrieve blastall path
process_blastall_path=$(grep "BLASTALL_PATH=" $programs_config|grep -v "^#")
eval $process_blastall_path

#retrieve htseq-count
process_htseq_count=$(grep "HTSEQCOUNT=" $programs_config|grep -v "^#")
eval $process_htseq_count

#retrieve tophat
process_TOPHAT_BIN=$(grep "TOPHAT_BIN=" $programs_config|grep -v "^#")
eval $process_TOPHAT_BIN

#retrieve bedtools path
process_bowtie2_path=$(grep "BOWTIE2_PATH=" $programs_config|grep -v "^#")
eval $process_bowtie2_path

#retrieve Rscript
process_Rscript=$(grep "r_script=" $programs_config|grep -v "^#")
eval $process_Rscript

#retrieve output directory
process_output_dir=$(grep "WORK_DIR=" $programs_config|grep -v "^#")
eval $process_output_dir

#retrieve project name
process_PROJECT_NAME=$(grep "PROJECT_NAME=" $programs_config|grep -v "^#")
eval $process_PROJECT_NAME

#retrieve annotation in bed format
process_ANNOTATION_BED_FILE=$(grep "ANNOTATION_BED_FILE=" $programs_config|grep -v "^#")
eval $process_ANNOTATION_BED_FILE

#retrieve annotation in gtf format
process_ANNOTATION_GTF_FILE=$(grep "ANNOTATION_GTF_FILE=" $programs_config|grep -v "^#")
eval $process_ANNOTATION_GTF_FILE

#retrieve ESTs file in bed format
process_EST_BED_FILE=$(grep "EST_BED_FILE=" $programs_config|grep -v "^#")
eval $process_EST_BED_FILE

#retrieve genome in fasta format
process_FASTA_FILE=$(grep "FASTA_FILE=" $programs_config|grep -v "^#")
eval $process_FASTA_FILE

#retrieve bowtie index
process_BOWTIE_INDEX=$(grep "BOWTIE_INDEX=" $programs_config|grep -v "^#")
eval $process_BOWTIE_INDEX

#retrieve Rscript
process_r_script=$(grep "r_script=" $programs_config|grep -v "^#")
eval $process_r_script

#retrieve proteins database
process_CPC_DATA=$(grep "CPC_DATA=" $programs_config|grep -v "^#")
eval $process_CPC_DATA

#retrieve PYTHONPATH
process_python_path=$(grep "python_path=" $programs_config|grep -v "^#")
eval $process_python_path

#retrieve the full path of some scripts
PATH_SCRIPTS=$( cd "$(dirname "$0")" ; pwd -P )/HoLDuP_scripts
echo -e "HolDuP scripts are in $PATH_SCRIPTS\n"
chmod 755 $PATH_SCRIPTS

# this is a in-house modified version of CPC with a parameter for the blast-database path, I only modified the run_predict.sh file
# (the database was too huge and I had to put it in a KDI space)
CPC_HOME=$( cd "$(dirname "$0")" ; pwd -P )/cpc
echo -e "CPC (coding potential calculator) scripts are in $CPC_HOME\n"
chmod 755 $CPC_HOME


###### check for errors (empty variables, not existing files). Remarks : one-line is used there to avoid to have codes on many lines, just for checks (sorry by advance) ##########

####programs/bins
#check samtools 
if [ "$samtools" == "" ];then samtools=$(which samtools);if [ "$samtools" == "" ];then echo "samtools not found !";exit ;else echo "samtools version is :"; mytest=$($samtools 2>&1|grep -i "^version") ;echo -e "$mytest\n"; fi ; else mytest=$($samtools 2>&1 |grep -i "^version");if [ "$mytest" == "" ];then echo "no tool detected for samtools !!"; exit; else echo -e "samtools version is $mytest\n";fi ;fi

#check bedtools bin
if [ "$BEDTOOLS_PATH" == "" ];then BEDTOOLS_PATH=$(which bedtools);if [ "$BEDTOOLS_PATH" == "" ];then echo "bedtools bin not found !"; exit ; else BEDTOOLS_PATH=$(dirname $BEDTOOLS_PATH);echo -e "bedtools bin is : ${BEDTOOLS_PATH}\n" ; fi; else echo -e "bedtools directory is : $BEDTOOLS_PATH\n" ;fi


#check blastall bin
if [ "$BLASTALL_PATH" == "" ];then BLASTALL_PATH=$(which blastall);if [ "$BLASTALL_PATH" == "" ];then echo "blastall bin not found !"; exit ; else BLASTALL_PATH=$(dirname $BLASTALL_PATH);echo -e "blastall bin is : ${BLASTALL_PATH}\n" ; fi; else echo -e "blastall directory is : $BLASTALL_PATH\n" ;fi


#adapt $PYTHONPATH first before to check htseq-count
if [ -z "$python_path" ];then echo -e "\npython_path (environment variable PYTHONPATH) is not modified for htseq-count, hope it's well installed on the cluster, otherwise, you will have to check \n";python_path="";fi


#we don't check htseq-count deeply like the others (code commented), because if the module numpy isn't available (it's often the case on the master node), we will have an error
#check htseq-count
#if [ "$HTSEQCOUNT" == "" ];then HTSEQCOUNT=$(which htseq-count);if [ "$HTSEQCOUNT" == "" ];then echo "htseq-count not found !";exit ;else echo "htseq-count version is : ";mytest=$($HTSEQCOUNT -help);mytest=$($HTSEQCOUNT -help|tail -n1);echo -e "$mytest\n" ; fi ; else mytest=$($HTSEQCOUNT -help);if [[ $(echo $?) -gt 0 ]];then echo "no tool detected for htseq-count, or maybe you should adapt $PYTHONPATH in the config file (variable python_path)!!" ; exit;else mytest=$($HTSEQCOUNT -help|tail -n1);echo -e "htseq-count version is : $mytest\n"; fi ;fi

if [ -z "$HTSEQCOUNT" ];then echo -e "\nHTSEQCOUNT is empty, check your config file !!!!\n";exit 1;fi

#check tophat2

if [ "$TOPHAT_BIN" == "" ];then TOPHAT_BIN=$(which tophat2);if [ "$TOPHAT_BIN" == "" ];then echo "tophat2 not found !";exit ;else echo "tophat2 version is :"; mytest=$($TOPHAT_BIN --version) ;echo -e "$mytest\n"; fi ; else mytest=$($TOPHAT_BIN --version);if [[ $(echo $?) -gt 0 ]];then echo "no tool detected for tophat2 !!"; exit; else echo -e "tophat2 version is $mytest\n";fi ;fi

#check bedtools bin
if [ "$BOWTIE2_PATH" == "" ];then BOWTIE2_PATH=$(which bowtie2);if [ "$BOWTIE2_PATH" == "" ];then echo "bowtie bin not found !"; exit ; else BOWTIE2_PATH=$(dirname $BOWTIE2_PATH);echo -e "bowtie2 bin is : ${BOWTIE2_PATH}\n" ; fi; else echo -e "bowtie2 directory is : $BOWTIE2_PATH\n" ;fi

#check Rscript
if [ "$r_script" == "" ];then r_script=$(which Rscript);if [ "$r_script" == "" ];then echo "Rscript not found !";exit ;else echo "Rscript version is :"; mytest=$($r_script --version 2>&1);echo -e "$mytest\n" ; fi ; else mytest=$($r_script --version 2>&1);if [[ $(echo $?) -gt 0 ]];then echo "no tool detected for Rscript !!" ; exit; else echo -e "Rscript version is : $mytest\n";fi ;fi



####supplied files

if [ "$WORK_DIR" == "" ];then echo -e "\nWORK_DIR is empty (output directory of HolDuP results), check your config file !!\n"; exit 1;else echo -e "\noutput dir. is $WORK_DIR\n";fi

if [ "$PROJECT_NAME" == "" ];then echo -e "\nPROJECT_NAME is empty, check your config file !!!!\n";exit 1;fi

if [ -z "$ANNOTATION_BED_FILE" ] || [ ! -f "$ANNOTATION_BED_FILE" ];then echo -e "\nANNOTATION_BED_FILE is empty or doesn't exit, check your config file !! !!\n";exit 1;fi

if [ -z "$ANNOTATION_GTF_FILE" ] || [ ! -f "$ANNOTATION_GTF_FILE" ];then echo -e "\nANNOTATION_GTF_FILE is empty or doesn't exit, check your config file !! !!\n";exit 1;fi

if [ -z "$EST_BED_FILE" ] || [ ! -f "$EST_BED_FILE" ];then echo -e "\nEST_BED_FILE is empty or doesn't exit, check your config file !! !!\n";exit 1;fi

if [ -z "$FASTA_FILE" ] || [ ! -f "$FASTA_FILE" ];then echo -e "\nFASTA_FILE is empty or doesn't exit, check your config file !! !!\n";exit 1;fi

if [ -z "$BOWTIE_INDEX" ];then echo -e "\nBOWTIE_INDEX is empty, check your config file !! !!\n";exit 1;fi

if [ -z "$CPC_DATA" ];then echo -e "\nCPC_DATA is empty, check your config file !! !!\n";exit 1;fi


#export the environment variables (we still have to export them in the subscripts that will be run by the qsub command though, but just like "export PATH=$PATH", no need to add all the directories already set)
export PATH="${PATH}:${BEDTOOLS_PATH}:${BLASTALL_PATH}:${BOWTIE2_PATH}"

export PYTHONPATH="${PYTHONPATH}:${python_path}"
echo -e "variable PYTHONPATH is $PYTHONPATH\n"


#CPC threshold
CPC_THRESHOLD=1.34365

# FASTQ_INDEX will be in job names and file names
# FASTQ_GROUP is used from the expression levels filtering
# FASTQ_FILE_R1 and FASTQ_FILE_R2 are paths to R1 and R2 fastq files
# these 4 arrays must have the same number of items
# FASTQ_INDEX[0]=prostate_cancer_tumor_1a
# FASTQ_GROUP[0]=prostate_cancer_tumor
# FASTQ_FILE_R1[0]=/bioinfo/guests/mdescrim/PROSTATE_CANCER/tmp_nfs/PROSTATE_CANCER/fastq/A65T1_R1.fastq.gz
# FASTQ_FILE_R2[0]=/bioinfo/guests/mdescrim/PROSTATE_CANCER/tmp_nfs/PROSTATE_CANCER/fastq/A65T1_R2.fastq.gz
# LIBRARY_TYPE[0]=fr-firststrand


FASTQ_INDEX=()
FASTQ_GROUP=()
FASTQ_FILE_R1=()
FASTQ_FILE_R2=()
LIBRARY_TYPE=()


#put the table in an bash array
#/data/tmp/mgabriel/design.txt
readarray my_config < $sample_design

i=1

while read line;do

 
		 #1st file is R1, 2nd file is R2, and so on, so we have to change the index only at R2
		 if [[ $(echo "$line" |awk '{print $3}'|grep -i "mate1"|wc -l) -gt 0 ]];then
		 
		   FASTQ_FILE_R1+=($(echo "$line" |awk '{print $1}'))
		 
		 
		 elif [[ $(echo "$line" |awk '{print $3}'|grep -i "mate2"|wc -l) -gt 0 ]];then
		 
		   FASTQ_FILE_R2+=($(echo "$line" |awk '{print $1}'))
		   
		   FASTQ_GROUP+=($(echo "$line" |awk '{print $2}'))

		   LIBRARY_TYPE+=($(echo "$line" |awk '{print $4}'))
		   
		   FASTQ_INDEX+=($(echo "$line" |awk -v my_index=$i '{print $2"_"my_index}'))
		   
		   echo -e "\nindex is $(echo "$line" |awk -v my_index=$i '{print $2"_"my_index}')"

		   i=$((i+1))
		 
		 fi
		 

 
          #process only non-empty line
done < <(echo "${my_config[*]}"|sed 's/^ //g'|LC_ALL=C grep -v "^$")

echo -e "\nfiles R1 are : \n$(echo -e ${FASTQ_FILE_R1[*]}|sed 's/ /\n/g')\n============\n"

echo -e "\nfiles R2 are : \n$(echo -e ${FASTQ_FILE_R2[*]}|sed 's/ /\n/g')\n============\n"


echo -e "\n====\nfastq group (${#FASTQ_GROUP[*]} ) : \n\n${FASTQ_GROUP[*]}
         \n====\nlibrary orientation (${#LIBRARY_TYPE[*]} ) : \n\n${LIBRARY_TYPE[*]}\n======\n"

#added by Marc G.
#exit

### IT STARTS HERE !!! ####
OUT_MAPPINGS=$WORK_DIR/$PROJECT_NAME/mappings
OUT_CONTIGS=$WORK_DIR/$PROJECT_NAME/contigs
OUT_COUNTS=$WORK_DIR/$PROJECT_NAME/counts
OUT_COUNTS_GENCODE=$WORK_DIR/$PROJECT_NAME/counts_gencode
OUT_CPC=$WORK_DIR/$PROJECT_NAME/CPC
BASH_SCRIPT_DIR=$WORK_DIR/$PROJECT_NAME/bash_script
TORQUE_LOG_DIR=$WORK_DIR/$PROJECT_NAME/torque_log

mkdir -p $BASH_SCRIPT_DIR
mkdir -p $TORQUE_LOG_DIR

N_SAMPLES=${#FASTQ_INDEX[@]}
IND_MAX=`expr $N_SAMPLES - 1`

# Mapping
for i in `seq 0 $IND_MAX`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_mapping_${FASTQ_INDEX[$i]}.$PROJECT_NAME.sh
  OUT_TOPHAT=$OUT_MAPPINGS/${FASTQ_INDEX[$i]}
  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N mapping_${FASTQ_INDEX[$i]}.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  if [ ! -f $OUT_TOPHAT/accepted_hits.bam ]
  then
    echo "#"PBS -l walltime=72:00:00 >> $bashFileName
    echo "#"PBS -l mem=60gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=12 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName
  
    echo -e "export PATH=$PATH\n" >>$bashFileName

 
    echo "mkdir -p $OUT_TOPHAT" >> $bashFileName
    echo "$TOPHAT_BIN -N 3 --read-edit-dist 3 -r 155 --mate-std-dev 80 -p 12 --library-type ${LIBRARY_TYPE[$i]} -g 1 -o $OUT_TOPHAT $BOWTIE_INDEX ${FASTQ_FILE_R1[$i]} ${FASTQ_FILE_R2[$i]}" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi
  
  JOB_ID_MAPPING[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_MAPPING[$i]} $bashFileName
done


# Mapping treatment
for i in `seq 0 $IND_MAX`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_traitement_mapping_${FASTQ_INDEX[$i]}.$PROJECT_NAME.sh
  OUT_TOPHAT=$OUT_MAPPINGS/${FASTQ_INDEX[$i]}
  
  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N traitement_mapping_${FASTQ_INDEX[$i]}.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_MAPPING[$i]} >> $bashFileName
  if [ ! -f $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.bam ]
  then
    echo "#"PBS -l walltime=12:00:00 >> $bashFileName
    echo "#"PBS -l mem=4gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName

    echo "$samtools view $OUT_TOPHAT/accepted_hits.bam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.noheader.sam" >> $bashFileName
    echo "$samtools view -H $OUT_TOPHAT/accepted_hits.bam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.header.sam" >> $bashFileName  
    echo "cat $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.header.sam $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.noheader.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.sam" >> $bashFileName  
    echo "$samtools view -bS $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.bam" >> $bashFileName
    echo "$samtools sort $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.bam $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique" >> $bashFileName
    echo "$samtools index $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.bam" >> $bashFileName
    echo "$samtools flagstat $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.bam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.flagstat" >> $bashFileName
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.header.sam" >> $bashFileName
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.noheader.sam" >> $bashFileName
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.sam" >> $bashFileName
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.unsorted.bam" >> $bashFileName
    echo "$samtools view -H $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.bam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.header.sam" >> $bashFileName
    echo "$samtools view $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.bam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.noheader.sam" >> $bashFileName
    echo "grep \"XS:A:+\" $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.noheader.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.noheader.sam" >> $bashFileName
    echo "grep \"XS:A:-\" $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.noheader.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.noheader.sam" >> $bashFileName
    echo "cat $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.header.sam $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.noheader.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.sam" >> $bashFileName
    echo "cat $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.header.sam $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.noheader.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.sam" >> $bashFileName
    echo "$samtools view -bS $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.bam" >> $bashFileName
    echo "$samtools view -bS $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.bam" >> $bashFileName
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.*.sam" >> $bashFileName
    echo "$samtools index $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.bam" >> $bashFileName
    echo "$samtools index $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.bam" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi
  
  JOB_ID_TRAITEMENT_MAPPING[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_TRAITEMENT_MAPPING[$i]} $bashFileName
done

# Recherche contigs
for i in `seq 0 $IND_MAX`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_recherche_contigs_${FASTQ_INDEX[$i]}.$PROJECT_NAME.sh
  OUT_CONTIG=$OUT_CONTIGS/${FASTQ_INDEX[$i]}
  
  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N recherche_contigs_${FASTQ_INDEX[$i]}.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_TRAITEMENT_MAPPING[$i]} >> $bashFileName
  if [ ! -f $OUT_CONTIG/contigs_total.bed ]
  then
    echo "#"PBS -l walltime=12:00:00 >> $bashFileName
    echo "#"PBS -l mem=4gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName

    echo "mkdir -p $OUT_CONTIG" >> $bashFileName
    echo "$BEDTOOLS_PATH/genomeCoverageBed -ibam $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.plus.unique.bam -d -g $PATH_SCRIPTS/genome.tab > $OUT_CONTIG/depth_plus.txt" >> $bashFileName
    echo "$BEDTOOLS_PATH/genomeCoverageBed -ibam $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.minus.unique.bam -d -g $PATH_SCRIPTS/genome.tab > $OUT_CONTIG/depth_minus.txt" >> $bashFileName
    echo "" >> $bashFileName
    echo "perl $PATH_SCRIPTS/covgenerate.pl $OUT_CONTIG/depth_plus.txt 101 + > $OUT_CONTIG/unite_plus.gff" >> $bashFileName
    echo "perl $PATH_SCRIPTS/gtfTobed.pl $OUT_CONTIG/unite_plus.gff > $OUT_CONTIG/unite_plus.bed" >> $bashFileName
    echo "" >> $bashFileName
    echo "perl $PATH_SCRIPTS/covgenerate.pl $OUT_CONTIG/depth_minus.txt 101 - > $OUT_CONTIG/unite_minus.gff" >> $bashFileName
    echo "perl $PATH_SCRIPTS/gtfTobed.pl $OUT_CONTIG/unite_minus.gff > $OUT_CONTIG/unite_minus.bed" >> $bashFileName
    echo "" >> $bashFileName

    echo "sort -k1,1 -k2,2n $OUT_CONTIG/unite_plus.bed > $OUT_CONTIG/unite_sorted_plus.bed" >> $bashFileName
    echo "$BEDTOOLS_PATH/mergeBed -s -i $OUT_CONTIG/unite_sorted_plus.bed -scores min -d 100 > $OUT_CONTIG/unite_sorted_merged_plus.bed" >> $bashFileName
    echo "awk '{print \$1\"\t\"\$2\"\t\"\$3\"\tcontig_plus_\"sprintf(\"%021d\",NR)\"_${FASTQ_INDEX[$i]}\t\"\$4\"\t\"\$5}' $OUT_CONTIG/unite_sorted_merged_plus.bed > $OUT_CONTIG/contigs_plus.bed" >> $bashFileName
    echo "sort -k1,1 -k2,2n $OUT_CONTIG/unite_minus.bed > $OUT_CONTIG/unite_sorted_minus.bed" >> $bashFileName
    echo "$BEDTOOLS_PATH/mergeBed -s -i $OUT_CONTIG/unite_sorted_minus.bed -scores min -d 100 > $OUT_CONTIG/unite_sorted_merged_minus.bed" >> $bashFileName
    echo "awk '{print \$1\"\t\"\$2\"\t\"\$3\"\tcontig_minus_\"sprintf(\"%021d\",NR)\"_${FASTQ_INDEX[$i]}\t\"\$4\"\t\"\$5}' $OUT_CONTIG/unite_sorted_merged_minus.bed > $OUT_CONTIG/contigs_minus.bed" >> $bashFileName
    echo "cat $OUT_CONTIG/contigs_plus.bed $OUT_CONTIG/contigs_minus.bed > $OUT_CONTIG/contigs_total.bed" >> $bashFileName
    echo "rm -f $OUT_CONTIG/depth_plus.txt $OUT_CONTIG/depth_minus.txt" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi
  
  JOB_ID_RECHERCHE_CONTIGS[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_RECHERCHE_CONTIGS[$i]} $bashFileName
done

# extrait les contigs qui ne sont pas dans l'annotation (novel contigs)
for i in `seq 0 $IND_MAX`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_extract_novel_contigs_${FASTQ_INDEX[$i]}.$PROJECT_NAME.sh
  OUT_CONTIG=$OUT_CONTIGS/${FASTQ_INDEX[$i]}
  
  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N extract_novel_contigs_${FASTQ_INDEX[$i]}.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_RECHERCHE_CONTIGS[$i]} >> $bashFileName
  if [ ! -f $OUT_CONTIG/contigs_sans_sens.length_transcripts.txt ]
  then
    echo "#"PBS -l walltime=12:00:00 >> $bashFileName
    echo "#"PBS -l mem=4gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName

    echo "$BEDTOOLS_PATH/intersectBed -v -s -a $OUT_CONTIG/contigs_total.bed -b $ANNOTATION_BED_FILE > $OUT_CONTIG/contigs_sans_sens.bed" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi
  
  JOB_ID_EXTRACT_NOVEL_CONTIGS[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_EXTRACT_NOVEL_CONTIGS[$i]} $bashFileName
done

# merge des novel contigs de toutes les banques
WAIT_QSUB="depend=afterok"
for i in `seq 0 $IND_MAX`
do
  WAIT_QSUB=$WAIT_QSUB:${JOB_ID_EXTRACT_NOVEL_CONTIGS[$i]}
done

bashFileName=$BASH_SCRIPT_DIR/bash_merge_novel_contigs.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N merge_novel_contigs.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W $WAIT_QSUB >> $bashFileName
if [ ! -f $OUT_CONTIGS/contigs_sans_sens.merged.bed ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=30gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "cat $OUT_CONTIGS/*/contigs_sans_sens.bed > $OUT_CONTIGS/contigs_sans_sens.cat.bed" >> $bashFileName
  echo "sort -k1,1 -k2,2n $OUT_CONTIGS/contigs_sans_sens.cat.bed > $OUT_CONTIGS/contigs_sans_sens.cat.sorted.bed" >> $bashFileName
  echo "$BEDTOOLS_PATH/mergeBed -s -i $OUT_CONTIGS/contigs_sans_sens.cat.sorted.bed -scores min -d 100 | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\tcontig_\"sprintf(\"%021d\",NR)\"_merged\\t\"\$4\"\\t\"\$5}' > $OUT_CONTIGS/contigs_sans_sens.merged.bed" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_MERGE_NOVEL_CONTIGS=$(qsub $bashFileName)
echo ${JOB_ID_MERGE_NOVEL_CONTIGS} $bashFileName

# extrait les merged novel contigs qui ne sont pas dans l'annotation
# et classe antisens/intergéniques
bashFileName=$BASH_SCRIPT_DIR/bash_class_as_int_contigs.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class_as_int_contigs.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_MERGE_NOVEL_CONTIGS} >> $bashFileName
if [ ! -f $OUT_CONTIGS/contigs_sans_sens.inter.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=12gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "$BEDTOOLS_PATH/intersectBed -v -s -a $OUT_CONTIGS/contigs_sans_sens.merged.bed -b $ANNOTATION_BED_FILE > $OUT_CONTIGS/contigs_sans_sens.novel.bed" >> $bashFileName
  echo "$BEDTOOLS_PATH/intersectBed -wo -S -a $ANNOTATION_BED_FILE -b $OUT_CONTIGS/contigs_sans_sens.novel.bed > $OUT_CONTIGS/contigs_sans_sens.anti.bed" >> $bashFileName
  echo "$BEDTOOLS_PATH/intersectBed -v -a $OUT_CONTIGS/contigs_sans_sens.novel.bed -b $ANNOTATION_BED_FILE > $OUT_CONTIGS/contigs_sans_sens.inter.bed" >> $bashFileName
  echo "awk '{print \$7\"\\t.\\ttranscript\\t\"\$8+1\"\\t\"\$9\"\\t.\\t\"\$12\"\\t.\\tgene_id \\\"\"\$10\"\\\";gene_type \\\"lncRNA\\\"\"}' $OUT_CONTIGS/contigs_sans_sens.anti.bed | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.anti.gtf" >> $bashFileName
  echo "awk '{print \$1\"\\t.\\ttranscript\\t\"\$2+1\"\\t\"\$3\"\\t.\\t\"\$6\"\\t.\\tgene_id \\\"\"\$4\"\\\";gene_type \\\"lncRNA\\\"\"}' $OUT_CONTIGS/contigs_sans_sens.inter.bed | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.inter.gtf" >> $bashFileName
  echo "awk '{print \$10\"\\t\"\$9-\$8}' $OUT_CONTIGS/contigs_sans_sens.anti.bed | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.anti.length_transcripts.txt" >> $bashFileName
  echo "awk '{print \$4\"\\t\"\$3-\$2}' $OUT_CONTIGS/contigs_sans_sens.inter.bed | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.inter.length_transcripts.txt" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS_AS_INT_CONTIGS=$(qsub $bashFileName)
echo ${JOB_ID_CLASS_AS_INT_CONTIGS} $bashFileName

# convertit le BED gencode en GTF
bashFileName=$BASH_SCRIPT_DIR/bash_convert_gencode_non_coding_bed2gtf.$PROJECT_NAME.sh
echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N convert_gencode_non_coding_bed2gtf.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
if [ ! -f $OUT_COUNTS_GENCODE/gencode.protein_coding.listID.txt ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=1gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "mkdir -p $OUT_COUNTS_GENCODE" >> $bashFileName
  echo "awk '{print \$1\"\\t.\\ttranscript\\t\"\$2+1\"\\t\"\$3\"\\t.\\t\"\$6\"\\t.\\tgene_id \\\"\"\$4\"\\\"\"}' $ANNOTATION_BED_FILE > $OUT_COUNTS_GENCODE/gencode.gtf" >> $bashFileName
  echo "awk '{print \$4\"\\t\"\$3-\$2}' $ANNOTATION_BED_FILE > $OUT_COUNTS_GENCODE/gencode.length_transcripts.txt" >> $bashFileName
  echo "perl $PATH_SCRIPTS/gtf_to_gene_id_gene_type.pl $ANNOTATION_GTF_FILE > $OUT_COUNTS_GENCODE/gencode.gene_id_gene_type.txt" >> $bashFileName
  echo "grep protein_coding $OUT_COUNTS_GENCODE/gencode.gene_id_gene_type.txt | awk '{print \$1}' > $OUT_COUNTS_GENCODE/gencode.protein_coding.listID.txt" >> $bashFileName  
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CONVERT_GENCODE=$(qsub $bashFileName)
echo ${JOB_ID_CONVERT_GENCODE} $bashFileName

# compte les reads pour chaque contigs anti inter et pour gencode
for i in `seq 0 $IND_MAX`
do
  bashFileName=${BASH_SCRIPT_DIR}/bash_read_count_${FASTQ_INDEX[$i]}.$PROJECT_NAME.sh
  OUT_COUNT=${OUT_COUNTS}/${FASTQ_INDEX[$i]}
  htseqOut_anti[$i]=${OUT_COUNT}/count.contigs_sans_sens.anti.out.txt
  htseqErr_anti[$i]=${OUT_COUNT}/count.contigs_sans_sens.anti.err.txt
  htseqOut_inter[$i]=${OUT_COUNT}/count.contigs_sans_sens.inter.out.txt
  htseqErr_inter[$i]=${OUT_COUNT}/count.contigs_sans_sens.inter.err.txt
  OUT_COUNT_GENCODE=${OUT_COUNTS_GENCODE}/${FASTQ_INDEX[$i]}
  htseqOut_gencode[$i]=${OUT_COUNT_GENCODE}/count.gencode.out.txt
  htseqErr_gencode[$i]=${OUT_COUNT_GENCODE}/count.gencode.err.txt
  
  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N read_count_${FASTQ_INDEX[$i]}.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_CLASS_AS_INT_CONTIGS}:${JOB_ID_CONVERT_GENCODE} >> $bashFileName
  if [ ! -f ${htseqErr_inter[$i]} ]
  then
    echo "#"PBS -l walltime=24:00:00 >> $bashFileName
    echo "#"PBS -l mem=4gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName

    #added by Marc G. (if not, htseq-count cannot load the module !)
    echo -e "export PYTHONPATH=$PYTHONPATH\n" >>$bashFileName
    
    echo "mkdir -p $OUT_COUNT" >> $bashFileName
    echo "mkdir -p $OUT_COUNT_GENCODE" >> $bashFileName
    echo "$samtools view $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.bam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.sam" >> $bashFileName
    echo "perl $PATH_SCRIPTS/prep_reads_for_counts.pl $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.sam > $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.groomed.sam" >> $bashFileName
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.sam" >> $bashFileName
    echo "$HTSEQCOUNT -t transcript -i gene_id -m intersection-nonempty $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.groomed.sam $OUT_CONTIGS/contigs_sans_sens.anti.gtf 1> ${htseqOut_anti[$i]} 2> ${htseqErr_anti[$i]}" >> $bashFileName 
    echo "$HTSEQCOUNT -t transcript -i gene_id -m intersection-nonempty $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.groomed.sam $OUT_CONTIGS/contigs_sans_sens.inter.gtf 1> ${htseqOut_inter[$i]} 2> ${htseqErr_inter[$i]}" >> $bashFileName 
    echo "$HTSEQCOUNT -t transcript -i gene_id -m intersection-nonempty $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.groomed.sam $OUT_COUNTS_GENCODE/gencode.gtf 1> ${htseqOut_gencode[$i]} 2> ${htseqErr_gencode[$i]}" >> $bashFileName 
    echo "rm $OUT_MAPPINGS/${FASTQ_INDEX[$i]}.unique.groomed.sam" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi
  
  JOB_ID_READ_COUNT[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_READ_COUNT[$i]} $bashFileName
done

# calcule les RPKM
WAIT_QSUB="depend=afterok"
for i in `seq 0 $IND_MAX`
do
  WAIT_QSUB=${WAIT_QSUB}:${JOB_ID_READ_COUNT[$i]}
  OUT_COUNT=${OUT_COUNTS}/${FASTQ_INDEX[$i]}
  OUT_COUNT_GENCODE=${OUT_COUNTS_GENCODE}/${FASTQ_INDEX[$i]}
  rpkmOut_anti[$i]=${OUT_COUNT}/rpkm.contigs_sans_sens.anti.out.txt
  rpkmOut_inter[$i]=${OUT_COUNT}/rpkm.contigs_sans_sens.inter.out.txt
  rpkmOut_gencode[$i]=${OUT_COUNT_GENCODE}/rpkm.gencode.out.txt
done

bashFileName=$BASH_SCRIPT_DIR/bash_calcul_rpkm.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N calcul_rpkm.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W $WAIT_QSUB >> $bashFileName

if [ ! -f $OUT_CONTIGS/contigs_sans_sens.anti.length_transcripts.txt ]
then
  echo "#"PBS -l walltime=48:00:00 >> $bashFileName
  echo "#"PBS -l mem=16gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName
  echo "$r_script $PATH_SCRIPTS/calcul_rpkm.R ${htseqOut_anti[*]} ${rpkmOut_anti[*]} $OUT_CONTIGS/contigs_sans_sens.anti.length_transcripts.txt" >> $bashFileName
  echo "$r_script $PATH_SCRIPTS/calcul_rpkm.R ${htseqOut_inter[*]} ${rpkmOut_inter[*]} $OUT_CONTIGS/contigs_sans_sens.inter.length_transcripts.txt" >> $bashFileName
  echo "$r_script $PATH_SCRIPTS/calcul_rpkm.R ${htseqOut_gencode[*]} ${rpkmOut_gencode[*]} $OUT_COUNTS_GENCODE/gencode.length_transcripts.txt" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CALCUL_RPKM=$(qsub $bashFileName)
echo ${JOB_ID_CALCUL_RPKM} $bashFileName

# cherche les transcrit du type voulu
for i in `seq 0 $IND_MAX`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_filtre_transcrit_type_voulu_${FASTQ_INDEX[$i]}.$PROJECT_NAME.sh
  OUT_COUNT_GENCODE=${OUT_COUNTS_GENCODE}/${FASTQ_INDEX[$i]}
  htseqOut_gencode_protein_coding[$i]=${OUT_COUNT_GENCODE}/count.gencode.protein_coding.out.txt
  rpkmOut_gencode_protein_coding[$i]=${OUT_COUNT_GENCODE}/rpkm.gencode.protein_coding.out.txt
  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N filtre_transcrit_type_voulu_${FASTQ_INDEX[$i]}.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_CALCUL_RPKM} >> $bashFileName
  if [ ! -f ${rpkmOut_gencode_protein_coding[$i]} ]
  then
    echo "#"PBS -l walltime=2:00:00 >> $bashFileName
    echo "#"PBS -l mem=1gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName
    
    echo "grep -Fwf $OUT_COUNTS_GENCODE/gencode.protein_coding.listID.txt ${htseqOut_gencode[$i]} > ${htseqOut_gencode_protein_coding[$i]}" >> $bashFileName
    echo "grep -Fwf $OUT_COUNTS_GENCODE/gencode.protein_coding.listID.txt ${rpkmOut_gencode[$i]} > ${rpkmOut_gencode_protein_coding[$i]}" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi
  JOB_ID_GENCODE_FILTRE_TRANSCRIPT_VOULU[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_GENCODE_FILTRE_TRANSCRIPT_VOULU[$i]} $bashFileName
done

WAIT_QSUB="depend=afterok:${JOB_ID_CALCUL_RPKM}"
for i in `seq 0 $IND_MAX`
do
  WAIT_QSUB=$WAIT_QSUB:${JOB_ID_GENCODE_FILTRE_TRANSCRIPT_VOULU[$i]}
done

#fait le filtre d'expression et produit les class3w
bashFileName=$BASH_SCRIPT_DIR/bash_filtre_rpkm_quantile.$PROJECT_NAME.sh
echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N filtre_rpkm_quantile.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W $WAIT_QSUB >> $bashFileName
if [ ! -f $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.0.1.proper.bed ]
then
  echo "#"PBS -l walltime=48:00:00 >> $bashFileName
  echo "#"PBS -l mem=16gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do

    #added by Marc G.
    echo -e "rpkmOut_gencode_protein_coding : \n${rpkmOut_gencode_protein_coding[*]}\n\nhtseqOut_anti : \n${htseqOut_anti[*]}\n\nFASTQ_GROUP : \n${FASTQ_GROUP[*]}\n\nrpkmOut_anti : \n${rpkmOut_anti[*]}\n\n==========================\n\n"
 

    echo "$r_script $PATH_SCRIPTS/filtre_expression_ajuste_par_quantile.R 10 $RPKM_QUANTILE_THRESHOLD ${rpkmOut_gencode_protein_coding[*]} ${htseqOut_anti[*]} ${FASTQ_GROUP[*]} ${rpkmOut_anti[*]} > $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.listID.txt" >> $bashFileName
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.listID.txt $OUT_CONTIGS/contigs_sans_sens.anti.gtf > $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName

    echo "$r_script $PATH_SCRIPTS/filtre_expression_ajuste_par_quantile.R 10 $RPKM_QUANTILE_THRESHOLD ${rpkmOut_gencode_protein_coding[*]} ${htseqOut_inter[*]} ${FASTQ_GROUP[*]} ${rpkmOut_inter[*]} > $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.listID.txt" >> $bashFileName
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.listID.txt $OUT_CONTIGS/contigs_sans_sens.inter.gtf > $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName

    echo "cp $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.gtf $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
    echo "cp $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.gtf $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
  echo "perl $PATH_SCRIPTS/gtf2bed.pl $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.0.1.gtf > $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.0.1.proper.bed" >> $bashFileName
  echo "perl $PATH_SCRIPTS/gtf2bed.pl $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.0.1.gtf > $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.0.1.proper.bed" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

#added by Marc G.
#exit

JOB_FILTRE_QUANTILE_RPKM=$(qsub $bashFileName)
echo ${JOB_FILTRE_QUANTILE_RPKM} $bashFileName

# compare to junctions obtained by the mappings
bashFileName=$BASH_SCRIPT_DIR/bash_compare_quantile_to_junctions.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N compare_quantile_to_junctions.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_FILTRE_QUANTILE_RPKM} >> $bashFileName
if [ ! -f $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_junctions.listID.0.1.txt ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "cat $OUT_MAPPINGS/*/junctions.bed |cut -f 1,2,3,6 |sort | uniq | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\tJUNC\"NR\"\\t0\\t\"\$4}' > $OUT_MAPPINGS/junctions_totales.bed" >> $bashFileName
  echo "${BEDTOOLS_PATH}/intersectBed -s -f 1 -wo -a $OUT_MAPPINGS/junctions_totales.bed -b $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.0.1.proper.bed | cut -f 10 | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_junctions.listID.0.1.txt" >> $bashFileName
  echo "${BEDTOOLS_PATH}/intersectBed -s -f 1 -wo -a $OUT_MAPPINGS/junctions_totales.bed -b $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.0.1.proper.bed | cut -f 10 | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_junctions.listID.0.1.txt" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_COMPARE_QUANTILE_TO_JUNCTIONS=$(qsub $bashFileName)
echo ${JOB_ID_COMPARE_QUANTILE_TO_JUNCTIONS} $bashFileName

#produit les classe2aw
bashFileName=$BASH_SCRIPT_DIR/bash_class2aw_quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class2aw_quantile.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_COMPARE_QUANTILE_TO_JUNCTIONS} >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2aw.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_junctions.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class2aw.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
  
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_junctions.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2aw.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS_2AW_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_CLASS_2AW_QUANTILE} $bashFileName

# compare to ESTs
bashFileName=$BASH_SCRIPT_DIR/bash_compare_quantile_to_EST.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N compare_quantile_to_EST.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_FILTRE_QUANTILE_RPKM} >> $bashFileName
if [ ! -f $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_EST.listID.0.1.txt ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "${BEDTOOLS_PATH}/intersectBed -s -f 1 -wo -a $EST_BED_FILE -b $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.0.1.proper.bed | awk '\$10>1' | cut -f 16 | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_EST.listID.0.1.txt" >> $bashFileName
  echo "${BEDTOOLS_PATH}/intersectBed -s -f 1 -wo -a $EST_BED_FILE -b $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.0.1.proper.bed | awk '\$10>1' | cut -f 16 | sort | uniq > $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_EST.listID.0.1.txt" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_COMPARE_QUANTILE_TO_EST=$(qsub $bashFileName)
echo ${JOB_ID_COMPARE_QUANTILE_TO_EST} $bashFileName

#produit les classe2bw
bashFileName=$BASH_SCRIPT_DIR/bash_class2bw_quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class2bw_quantile.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_COMPARE_QUANTILE_TO_EST} >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2bw.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class2bw.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
  
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2bw.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS_2BW_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_CLASS_2BW_QUANTILE} $bashFileName

# make class1w lncRNA files
bashFileName=$BASH_SCRIPT_DIR/bash_class1w_quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class1w_quantile.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_CLASS_2AW_QUANTILE}:${JOB_ID_CLASS_2BW_QUANTILE} >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class1w.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class2aw.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class1w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
  
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2aw.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class1w.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS1W_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_CLASS1W_QUANTILE} $bashFileName

# extract sequence for each transcript
bashFileName=$BASH_SCRIPT_DIR/bash_extract_sequence_contigs.class3w.quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N extract_sequence_contigs.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_FILTRE_QUANTILE_RPKM} >> $bashFileName
if [ ! -f $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "" >> $bashFileName
  echo "perl $PATH_SCRIPTS/gtf2bed.pl $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3w.quantile.0.1.gtf > $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.bed" >> $bashFileName
  echo "perl $PATH_SCRIPTS/gtf2bed.pl $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3w.quantile.0.1.gtf > $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.bed" >> $bashFileName

  echo "" >> $bashFileName
  echo "$BEDTOOLS_PATH/bedtools getfasta -s -name -fi $FASTA_FILE -bed $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.bed -fo $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa" >> $bashFileName
  echo "cat $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa | tr '\\n' ' ' | tr '>' '\\n' | tail -n +2 > $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp" >> $bashFileName
  echo "awk -F, '{print > \"$OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp\"NR%100}' $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp" >> $bashFileName
  echo "rm $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp" >> $bashFileName
  echo "for i in \`seq 0 99\`" >> $bashFileName
  echo "do" >> $bashFileName
  echo "  if [ -e $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp\$i ]" >> $bashFileName
  echo "  then" >> $bashFileName
  echo "    sed 's/ $//' $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp\$i | sed 's/^/>/' | tr ' ' '\\n' > $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.\$i" >> $bashFileName
  echo "    rm $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.tmp\$i" >> $bashFileName
  echo "  fi" >> $bashFileName
  echo "done" >> $bashFileName

  echo "" >> $bashFileName
  echo "$BEDTOOLS_PATH/bedtools getfasta -s -name -fi $FASTA_FILE -bed $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.bed -fo $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa" >> $bashFileName
  echo "cat $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa | tr '\\n' ' ' | tr '>' '\\n' | tail -n +2 > $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp" >> $bashFileName
  echo "awk -F, '{print > \"$OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp\"NR%100}' $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp" >> $bashFileName
  echo "rm $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp" >> $bashFileName
  echo "for i in \`seq 0 99\`" >> $bashFileName
  echo "do" >> $bashFileName
  echo "  if [ -e $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp\$i ]" >> $bashFileName
  echo "  then" >> $bashFileName
  echo "    sed 's/ $//' $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp\$i | sed 's/^/>/' | tr ' ' '\\n' > $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.\$i" >> $bashFileName
  echo "    rm $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.tmp\$i" >> $bashFileName
  echo "  fi" >> $bashFileName
  echo "done" >> $bashFileName
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_EXTRACT_SEQUENCE_CONTIGS_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_EXTRACT_SEQUENCE_CONTIGS_QUANTILE} $bashFileName

# coding potential calculator anti
for i in `seq 0 99`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_coding_potential_calculator_anti_quantile_$i.$PROJECT_NAME.sh

  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N coding_potential_calculator_anti_$i.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_EXTRACT_SEQUENCE_CONTIGS_QUANTILE} >> $bashFileName
  if [ ! -f $OUT_CPC/anti_quantile_01_$i/CPC_results.tab ]
  then
    echo "#"PBS -l walltime=12:00:00 >> $bashFileName
    echo "#"PBS -l mem=15gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName

    #added by Marc G.
    echo -e "export PATH=$PATH\n" >>$bashFileName

    echo "mkdir -p $OUT_CPC/anti_quantile_01_$i/tmp" >> $bashFileName
    echo "if [ -e $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.$i ]" >> $bashFileName
    echo "then" >> $bashFileName
    echo "  $CPC_HOME/bin/run_predict.sh $OUT_CONTIGS/contigs_sans_sens.anti.class3w.quantile.0.1.fa.$i $OUT_CPC/anti_quantile_01_$i/CPC_results.tab $OUT_CPC/anti_quantile_01_$i/tmp $OUT_CPC/anti_quantile_01_$i/result_evidence_plot $CPC_DATA" >> $bashFileName
    echo "fi" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi

  JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE[$i]} $bashFileName
done

# coding potential calculator inter
for i in `seq 0 99`
do
  bashFileName=$BASH_SCRIPT_DIR/bash_coding_potential_calculator_inter_quantile_$i.$PROJECT_NAME.sh

  echo "#"!/bin/bash > $bashFileName
  echo "#"PBS -N coding_potential_calculator_inter_$i.$PROJECT_NAME >> $bashFileName
  echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
  echo "#"PBS -W depend=afterok:${JOB_ID_EXTRACT_SEQUENCE_CONTIGS_QUANTILE} >> $bashFileName
  if [ ! -f $OUT_CPC/inter_quantile_01_$i/CPC_results.tab ]
  then
    echo "#"PBS -l walltime=12:00:00 >> $bashFileName
    echo "#"PBS -l mem=15gb >> $bashFileName
    echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
    echo "#"PBS -q batch >> $bashFileName
    echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
    echo "#"PBS -m ae >> $bashFileName
    echo "#"PBS -j oe >> $bashFileName

    #added by Marc G.
    echo -e "export PATH=$PATH\n" >>$bashFileName

    echo "mkdir -p $OUT_CPC/inter_quantile_01_$i/tmp" >> $bashFileName
    echo "if [ -e $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.$i ]" >> $bashFileName
    echo "then" >> $bashFileName
    echo "  $CPC_HOME/bin/run_predict.sh $OUT_CONTIGS/contigs_sans_sens.inter.class3w.quantile.0.1.fa.$i $OUT_CPC/inter_quantile_01_$i/CPC_results.tab $OUT_CPC/inter_quantile_01_$i/tmp $OUT_CPC/inter_quantile_01_$i/result_evidence_plot $CPC_DATA" >> $bashFileName
    echo "fi" >> $bashFileName
  else
    echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
  fi

  JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE[$i]=$(qsub $bashFileName)
  echo ${JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE[$i]} $bashFileName
done

# filtre les non-codants et produit class3 anti
WAIT_QSUB="depend=afterok"
for i in `seq 0 99`
do
  WAIT_QSUB=$WAIT_QSUB:${JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE[$i]}
done

bashFileName=$BASH_SCRIPT_DIR/bash_coding_potential_calculator_anti_quantile_filter.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N coding_potential_calculator_anti_filter.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W $WAIT_QSUB >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "cat $OUT_CPC/anti_quantile_01_*/CPC_results.tab > $OUT_CPC/anti_quantile_01_CPC_results.tab" >> $bashFileName
  echo "awk '\$4 <= $CPC_THRESHOLD {print \$1}' $OUT_CPC/anti_quantile_01_CPC_results.tab > $OUT_CPC/non_coding_contigs_anti.quantile.0.1.listID.txt" >> $bashFileName
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CPC/non_coding_contigs_anti.quantile.0.1.listID.txt $OUT_CONTIGS/contigs_sans_sens.anti.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE_FILTER=$(qsub $bashFileName)
echo ${JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE_FILTER} $bashFileName

# filtre les non-codants et produit class3 inter
WAIT_QSUB="depend=afterok"
for i in `seq 0 99`
do
  WAIT_QSUB=$WAIT_QSUB:${JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE[$i]}
done

bashFileName=$BASH_SCRIPT_DIR/bash_coding_potential_calculator_inter_quantile_filter.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N coding_potential_calculator_inter_filter.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W $WAIT_QSUB >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  echo "cat $OUT_CPC/inter_quantile_01_*/CPC_results.tab > $OUT_CPC/inter_quantile_01_CPC_results.tab" >> $bashFileName
  echo "awk '\$4 <= $CPC_THRESHOLD {print \$1}' $OUT_CPC/inter_quantile_01_CPC_results.tab > $OUT_CPC/non_coding_contigs_inter.quantile.0.1.listID.txt" >> $bashFileName
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CPC/non_coding_contigs_inter.quantile.0.1.listID.txt $OUT_CONTIGS/contigs_sans_sens.inter.expressed_quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE_FILTER=$(qsub $bashFileName)
echo ${JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE_FILTER} $bashFileName

# produit les class2a
bashFileName=$BASH_SCRIPT_DIR/bash_class2a_quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class2a.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE_FILTER}:${JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE_FILTER}:${JOB_ID_COMPARE_QUANTILE_TO_JUNCTIONS} >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2a.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_junctions.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class2a.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
  
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_junctions.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2a.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS_2A_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_CLASS_2A_QUANTILE} $bashFileName

# produit les class2b
bashFileName=$BASH_SCRIPT_DIR/bash_class2b_quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class2b.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_CODING_POTENTIAL_CALCULATOR_ANTI_QUANTILE_FILTER}:${JOB_ID_CODING_POTENTIAL_CALCULATOR_INTER_QUANTILE_FILTER}:${JOB_ID_COMPARE_QUANTILE_TO_EST} >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2b.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class3.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class2b.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
  
  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class3.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2b.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS_2B_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_CLASS_2B_QUANTILE} $bashFileName

# make class1 lncRNA files
bashFileName=$BASH_SCRIPT_DIR/bash_class1_quantile.$PROJECT_NAME.sh

echo "#"!/bin/bash > $bashFileName
echo "#"PBS -N class1.$PROJECT_NAME >> $bashFileName
echo "#"PBS -o $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -e $TORQUE_LOG_DIR >> $bashFileName
echo "#"PBS -W depend=afterok:${JOB_ID_CLASS_2A_QUANTILE}:${JOB_ID_CLASS_2B_QUANTILE} >> $bashFileName
if [ ! -f $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class1.quantile.0.1.gtf ]
then
  echo "#"PBS -l walltime=12:00:00 >> $bashFileName
  echo "#"PBS -l mem=4gb >> $bashFileName
  echo "#"PBS -l nodes=1:ppn=1 >> $bashFileName
  echo "#"PBS -q batch >> $bashFileName
  echo "#"PBS -M mdescrim@curie.fr >> $bashFileName
  echo "#"PBS -m ae >> $bashFileName
  echo "#"PBS -j oe >> $bashFileName

  for RPKM_QUANTILE_THRESHOLD in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
  do
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.anti.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class2a.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.anti.class1.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
    echo "grep -Fwf $OUT_CONTIGS/contigs_sans_sens.quantile.inter.with_EST.listID.0.1.txt $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class2a.quantile.$RPKM_QUANTILE_THRESHOLD.gtf > $WORK_DIR/$PROJECT_NAME/contigs_sans_sens.inter.class1.quantile.$RPKM_QUANTILE_THRESHOLD.gtf" >> $bashFileName
  done
else
  echo "sleep 180 > /dev/null 2>&1" >> $bashFileName
fi

JOB_ID_CLASS1_QUANTILE=$(qsub $bashFileName)
echo ${JOB_ID_CLASS1_QUANTILE} $bashFileName
