There are 8 parts in this README (from 1) to 8)), please read carefully at least until 7) included, to avoid unwanted behaviors, or to gain time !


###########################
### 1) Purpose :
############################

HoLDuP is a pipeline for ab initio transcripts assembly, from mapped paired-end reads, based on BEDTools facilities.

############################



###########################
### 2) Required dependancies :
############################

	- bash (version >= 4.3.46)

	- awk (version >= 4.1.3)

	- samtools (version 0.1.8)
	
	- tophat (version >= 2.1.1)
	
	- bowtie2 (version >= 2.2.5)

	- bedtools (version  >= 2.17)
	
	- htseq-count (version 0.6.1p1)
	
	- blastall (version 2.2.25)
	
	- cluster with Torque
	
	- uncompress cpc.tar.gz (coding potential calculator, http://cpc.cbi.pku.edu.cn) in the "HoLDuP_pipeline" directory

	
	

###############################



###########################
### 3) Required inputs :
############################
	
	- annotation in gtf format
	
	- annotation in bed6 format (convert the gtf above) ; the 4th column contains the gene id, one line can be for example : chr1	11868	14409	ENSG00000223972.5	.	+

	- ESTs file in uncompressed fasta format
	
		-> can be found here for hg38 http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/all_est.txt.gz
		-> can be found here for hg19 http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/all_est.txt.gz
		
    - uncompressed human genome in fasta format
		
	- bowtie index of the human genome (use a bowtie version compatible with the tophat2 you are using)
	
	- database of proteins from uniprot (ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz)
	

###############################


###############################
### 4) Processings by the user before usage :
###############################

You must create a database of the proteins (ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz), with the "formatdb" module of "blastall", in order to check the coding potential of the contigs :

uncompress the fasta file, then run :

/path/to/formatdb -p T -t prot_db-i uniref90.fasta

put the output database (files *phr, *pin, *psq, *pin, *pal) in a directory, and give this directory (full path) to the variable "CPC_DATA" in the config file "programs_config.txt" below in the Usage


###############################



###############################
### 4) Usage: 
###############################

on the master of the cluster (the pipeline will produce subscripts with PBS parameters as header, and these will be called with qsub commands by the pipeline, with dependancies between the jobs) :

./HoLDuP_pipeline.sh -f < sample design > -c < programs config >


 	Required arguments :

                  
                  -f < sample design  (you can adapt the file "sample_design.txt" supplied with the pipeline) >

                      		example of supplied file (1st column fastq file with full path, 2nd the condition, 3rd mate indication, 4th the library type):

                      
				/home/path/to/B67T11.R1.fastq.gz	prostate_tumor	mate1	fr-firststrand
                   
				/home/path/to/B67T11.R2.fastq.gz	prostate_tumor	mate2	fr-firststrand
                   
				/home/path/to/B67T14.R1.fastq.gz	prostate_normal	mate1	fr-firststrand
                   
				/home/path/to/B67T14.R2.fastq.gz	prostate_normal	mate2	fr-firststrand

                   
                                          	********* warnings *********
                   
                            		-> keep the order of mate1 & mate2 (don't put mate2 before mate1)
                             		-> between each column, we have tabulations (not spaces)
                   
                                          	***************************

                  -c < programs config (adapt the "programs_config.txt" supplied with the pipeline) >

                  	-> text file containing the location of the programs (or their bin)
                  	
                  	
    optional arguments :
    
                    -h < show Usage >                
                    
                    -v < show version >


	Results :

                  -  4 gtf files of intergenic & antisense contigs (class 1 & 3)
                  
                      -> contigs_sans_sens.anti.class3w.quantile.0.2.gtf
                      -> contigs_sans_sens.inter.class3w.quantile.0.2.gtf
                      -> contigs_sans_sens.anti.class1w.quantile.0.2.gtf
                      -> contigs_sans_sens.inter.class1w.quantile.0.2.gtf
                      
                      remarks : 
                      
                        - class 3 contigs : contigs with RPKM expression above the value of a quantile threshold (we use 0.2)
                      
                        - class 1 contigs : subset of class 3 ; contigs with RPKM expression above value of quantile threshold + presence of EST(s) + presence of junction(s)
                        

###############################



###############################
### 5)  useful directories
###############################

("WORK_DIR" is the variable you have to set in the config file "programs_config.txt" for the output directory of HolDuP results)
- $WORK_DIR

    ("PROJECT_NAME"  is the variable you have to set in the config file "programs_config.txt" for the name of the project ; it will be inside the output directory. So you could run many projects in this ouput directory)
	- $PROJECT_NAME 

				- bash_script/			-> location of the scripts ran by the qsub command in the pipeline (from these, you can find the log file name)
				
				- contigs/				-> directory where the contigs are processed (intermediate files)
				
				- counts			 	-> results from htseq-count on the inferred contigs - RPKM normalization
				
				- counts_gencode/		-> results from htseq-count on the official annotation - RPKM normalization
				
				- mappings/				-> directory of the BAM files for each sample
				
				- torque_log/			-> directory of the logs (we have log files for each submitted job by the qsub command in the pipeline). you can check here if something went wrong during the run.

###############################


                  
###############################
### 6) Description of the results (the interesting files are directly in ${WORK_DIR}/${PROJECT_NAME}/) : 
###############################

gtf files with prefix :

	- "contigs_sans_sens.anti"  -> antisense contigs in gtf format

	- "contigs_sans_sens.inter" -> intergenic contigs in gtf format
                  
gtf files with suffix :

    - "class3w"  -> contigs with RPKM expression above the value of a quantile threshold (the pipeline is called for 0.1 =< threshold <=0.9, but we use 0.2)
    
    - "class2aw" -> subset of 3w ; contigs with junction(s)
    
    - "class2bw" -> subset of 3w ; contigs with EST(s)
    
	- "class1w"  -> subset of 3w ; contigs with RPKM expression above the value of a quantile threshold + EST(s) + junctions
	
	
	This can be shematized like this :
	
	
							 ----
							| 3w | RPKM expression above the value of the quantile threshold (we use 0.2)
							 ----
							   |
							   v
			                  ----------------------------------
			                  |                                |
			                  v                                v
	                                ----                              ---- 
		                       |2bw | EST(s)                     |2aw | junction(s)
		                        ----                              ----
	                                  |                                 |
	                                   ----------------------------------
	                                                  |
	                                                  v
							 ----
							| 1w | RPKM expression above the value of the quantile threshold + EST(s) + junction(s)
							 ---- 
							 
	
Files without the "w" after the "class", mean their contigs with a value of coding potential above the threshold (1.34365 by default) have been filtered out ;
Usually we don't use them, because the part with the coding potential computing can be very long. So, the script can be stopped after the command "echo ${JOB_ID_CLASS1W_QUANTILE} $bashFileName"

We use the files with contigs above the value of 0.2 quantile of RPKM (4 files) :

contigs_sans_sens.anti.class3w.quantile.0.2.gtf
contigs_sans_sens.inter.class3w.quantile.0.2.gtf

contigs_sans_sens.anti.class1w.quantile.0.2.gtf
contigs_sans_sens.inter.class1w.quantile.0.2.gtf


###############################



###############################
### 7) Tips
###############################

 - The pipeline can be stopped after the command "echo ${JOB_ID_CLASS1W_QUANTILE} $bashFileName", because the rest will create the same files, excepted that the contigs with a coding potential above the threshold (variable "CPC_THRESHOLD" in the pipeline) will be filtered out, as said in 6) (files without "w" after the "class").

###############################
	


###############################

### 8) Improvements to do : 

###############################

- put the interesting 4 files in a dedicated directory (named "results" for example)

- update the tools :
           
           - replace htseq-count by featureCount (for speed) http://bioinf.wehi.edu.au/featureCounts/
           
           - replace blastall (this one is really old, not maintained anymore, and slow) by ncbi-blast ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
           
           - replace the bedtools version used (as the latest versions have some changes in the parameters compared to the 2.17 used there, some tests need to be done)
           
           - replace the samtools version used (as the latest versions have some changes in the parameters compared to the 0.1.8 used there, some tests need to be done)
           
           - replace the CPC version (coding potential calculator), by a faster one, CPC2 for example, as it's 1000 time faster as mentioned here https://www.ncbi.nlm.nih.gov/pubmed/28521017
           
           - create a version of the pipeline (or adapt this one) that could be run without a cluster
           
- simplify some inputs :

             - instead of the user, the annotation in bed6 format could be done by the pipeline, thanks to the supplied gtf annotation
             
- error logs :
            
             - add errors when some intermediate files are created, but are empty (or are not created at all)
         
###############################
