-include ./config.mk

FASTQ_DIR ?= ../00.fastq
FASTQ_EXT ?= .fastq.gz
FASTQC_MODULE ?= FastQC/0.11.4

FASTQ_FILES := $(wildcard $(FASTQ_DIR)/*$(FASTQ_EXT))
BASENAMES = $(subst $(FASTQ_EXT),,$(notdir $(FASTQ_FILES)))

SHELL:=/bin/bash

FASTQC_FILES = $(addprefix out/,$(addsuffix _fastqc.html, $(BASENAMES)))

fastqc: fastqc_prepare fastqc_run fastqc_summary

fastqc_prepare:
	@mkdir -p out

	@echo -e "\n\n# FASTQC RUN " >> mad.sh
	@date +'# %y.%m.%d %H:%M:%S' >> mad.sh
	@echo -e "module load $(FASTQC_MODULE)\n\n" >> mad.sh

	@echo "Found $(words $(FASTQ_FILES)) file(s) in $(FASTQ_DIR)"

fastqc_run: $(FASTQC_FILES)
	@echo "Finished running fastqc"

fastqc_summary: fastqc.summary.tsv
	@echo "Creating summary table"
	@mad table -o fastqc.mad.summary.tsv $(FASTQ_FILES)

fastqc.summary.tsv: ../sampleids.tsv
	./summary.py out

$(FASTQC_FILES): out/%_fastqc.html: $(FASTQ_DIR)/%$(FASTQ_EXT)
	@echo "Running fastqc for : " $<
	@echo "Name               : " $*
	@echo "Output html        : " $@
	@echo "        zip        : " out/$*_fastqc.zip
	$(eval $@_CL := fastqc -o out/ $<)
	module load $(FASTQC_MODULE) \
		&& ( $($@_CL) > out/$*_fastqc.log 2>&1 )

	echo "mad ta add \
			  --input $< \
				--cl '$($@_CL)' \
				--output html:$@ \
				--output zip:out/$*_fastqc.zip \
				--executable fastqc" \
			>> mad.sh

		echo "rat fastqc_zip out/$*_fastqc.zip \
		  	| mad apply_from_table -p fastqc_ - $<" \
			>> mad.sh
