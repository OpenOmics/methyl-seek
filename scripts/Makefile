##################################################################################################################
#### Makefile to merge FASTQ.GZ files across different lanes
#### 1. Copy this code to your project directory 
#### 2. Specify the path to folder containing sequencing fastq.gz files: sampleName_S1_L00*_R?_001.fastq.gz
#### 3. Specify the path to output directory path to create merged files. Directory will be created if not absent
#### 
####
##################################################################################################################


FQFILES_DIR=/data/NHLBI_IDSS/rawdata/NHLBI-103/MKJ-9666
MERGEFILES_DIR=/data/NHLBI_IDSS/rawdata/NHLBI-103/mergedNHLBI-103/MKJ-9666

#mkdir -p $(MERGEFILES_DIR)
FILELIST=`ls -1 ${FQFILES_DIR}/ | grep '_R2_001.fastq.gz' | sed 's/_R2_001.fastq.gz//g'`
SAMPLELIST=`ls -1 ${FQFILES_DIR}/ | grep '_R2_001.fastq.gz' | sed 's/_R2_001.fastq.gz//g' | sed s'/_L00/\t/g' | cut -f 1 | sort | uniq`
LANES=`ls -1 ${FQFILES_DIR}/ | grep '_R2_001.fastq.gz' | sed 's/_R2_001.fastq.gz//g' | sed s'/_L00/\t/g' | cut -f 2 | sort | uniq`


all:
	@echo "Makefile to merge *.fastq.gz files across different sequencing lanes"
	@echo ""
	@echo "SAMPLE list as follows:    "
	@for sample in $(SAMPLELIST); do \
		echo "$${sample}"; \
		done ;
	@echo ""
	@echo ""
	@echo "Number of LANES:    "
	@echo "${LANES}"
	@echo ""
	@echo ""
	@echo ""
	@echo "make merge_fastqs"
	@echo "make sample_file"
 
 
merge_fastqs:
	@echo "Merge sample fastq.gz files across different lanes"
	@test -d $(MERGEFILES_DIR) || mkdir -p $(MERGEFILES_DIR)
	@echo "#!/bin/bash" > merge.sh
	@for sample in $(SAMPLELIST); do \
		echo "cat $(FQFILES_DIR)/$${sample}_L00*_R1_001.fastq.gz > $(MERGEFILES_DIR)/$${sample}.R1.fastq.gz" >> merge.sh;\
		echo "cat $(FQFILES_DIR)/$${sample}_L00*_R2_001.fastq.gz > $(MERGEFILES_DIR)/$${sample}.R2.fastq.gz" >> merge.sh;\
		done ;
	@echo "swarm --gb-per-process 10 --threads-per-process 2 --time 2-00:00:00 --bundle 2 merge.sh"


sample_file:
	@echo "Make sample.txt file for methyl-seek"
	@echo -e "samples\tgroup" > sample.txt
	@for sample in $(SAMPLELIST); do \
		echo "$${sample}" >> sample.txt;\
		done ;
	@echo "sample.txt file is created!"

