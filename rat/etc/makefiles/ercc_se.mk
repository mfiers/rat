-include ./config.mk
-include ./rat_core.mk

FASTQ_DIR ?= ../00.fastq
FASTQ_EXT ?= .fastq.gz
BOWTIE_MODULE ?= Bowtie2/2.2.6-foss-2014a
SAMTOOLS_MODULE ?= SAMtools/1.3-foss-2014a
MODULES ?= SAMtools/1.3-foss-2014a Bowtie2/2.2.6-foss-2014a

FASTQ_FILES = $(wildcard $(FASTQ_DIR)/*$(FASTQ_EXT))
BASENAMES = $(subst $(FASTQ_EXT),,$(notdir $(FASTQ_FILES)))

SHELL := /bin/bash

COUNT_FILES = $(addprefix out/, $(addsuffix .count, $(BASENAMES)))

ERCC_URL = https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
MAGPIE_GROUP = singlecell_ercc
.PHONY: ercc_db run prepare

all: prepare map ercc_summary

prepare: ercc_db
	@echo -e "\n\n# STAR RUN $(date)" >> mad.sh
	@echo "module load $(SAMTOOLS_MODULE)" >> mad.sh
	@echo -e "module load $(BOWTIE_MODULE)\n\n" >> mad.sh
	@mkdir -p out
	@echo "Found $(words $(FASTQ_FILES)) file(s) in $(FASTQ_DIR)"


ercc_db: db/ercc.1.bt2

ercc_summary: ercc.group.tsv

ercc.group.tsv: ../sampleids.tsv
	./summarize.py

db/ercc.1.bt2:
	module load $(BOWTIE_MODULE) \
	  && mkdir -p db \
		&& cd db \
		&& wget -O ercc.zip $(ERCC_URL) \
		&& unzip -o ercc.zip \
		&& bowtie2-build ERCC92.fa ercc

		echo -e '\n\n# $@' >> mad.sh
		echo 'mad ta add \
					--executable bowtie2-build \
					--cl "bowtie2-build ERCC92.fa ercc" \
					--input ERCC92.fa \
					--output ercc.1.bt2 \
					--output ercc.2.bt2 \
					--output ercc.3.bt2 \
					--output ercc.4.bt2 \
					--output ercc.rev.1.bt2 \
					--output ercc.rev.2.bt2' \
				>> mad.sh

map: $(COUNT_FILES)

$(COUNT_FILES): out/%.count: $(FASTQ_DIR)/%$(FASTQ_EXT)
	@echo "Running bowtie for : " $<
	@echo "Name               : " $*
	@echo "Output bam         : " $@
	@echo "Bowtie module      : " $(BOWTIE_MODULE)
	@echo "Magpie             : " $(MAGPIE) "--" $(MCMD)
	$(eval $@_CL := \
			bowtie2 -p 8 \
			  --sensitive-local \
				-x db/ercc \
				-U $< \
				| samtools view -S -q 10 - \
				| cut -f 3 | sort | uniq -c \
				| sed "s/ *\(.*\) \(.*\)/\2\t\1/" \
				> $@ )

	module load $(BOWTIE_MODULE) \
		&& module load $(SAMTOOLS_MODULE) \
		&& $($@_CL)

	echo -e "\n\n# $@" >> mad.sh
	echo "mad ta add --executable samtools \
		  --executable bowtie2 \
			--cl '$($@_CL)' \
			--input $< \
			--output $@" \
		>> mad.sh
