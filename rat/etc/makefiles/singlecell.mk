-include ./config.mk
-include ./rat_core.mk

#? Input directory with parameters
FASTQ_DIR ?= 00.fastq

#? Extension of gipped fastq files
FASTQ_EXT ?= .fastq.gz

#? Regex to find sample ids
SAMPLE_REGEX ?= (.*).fastq.gz

#? R Module
R_MODULE ?= R/3.2.4-foss-2014a-noX

###end parameters

FASTQ_FILES := $(wildcard $(FASTQ_DIR)/*$(FASTQ_EXT))

.PHONY: report ercc star fastqc samplid check all

all: check sampleid fastqc star ercc report

check: $(addprefix check_,$(notdir $(FASTQ_FILES)))
	[[ -e $(FASTQ_DIR) ]]

check_%: $(FASTQ_DIR)/%
	@[[ ! -z $< ]]

sampleid: sampleids.tsv

sampleids.tsv: config.mk
	if [[ -e sampleids.tsv ]]; then \
		mv sampleids.tsv sampleids.backup.tsv; fi

	@echo "Assigning sampleids from filename"
	@echo "using regex: $(SAMPLE_REGEX)"
	@mad apply_using_filename_regex sample "$(SAMPLE_REGEX)" $(FASTQ_FILES) \
		>> sampleids.backup.tsv
	cat sampleids.backup.tsv | sort | uniq > sampleids.tsv
	-rm -rf sampleids.backup.tsv

clean_stat:
	-rm *.tsv
	-rm */*.tsv
	-rm *.html

fastqc:
	mkdir -p 05.fastqc
	cd 05.fastqc && rat mf install -f fastqc && $(MAKE)

star:
	mkdir -p 10.star
	cd 10.star && rat mf install -f starse && $(MAKE)

ercc:
	mkdir -p 20.ercc
	cd 20.ercc && rat mf install -f ercc_se && $(MAKE)

report: sc_report.ipynb
	k3 t ~/k3/brennecke_ercc.k3 20.ercc/ercc.group.tsv 10.star/counts.group.tsv
	jupyter nbconvert --ExecutePreprocessor.enabled=True --to=html $<
