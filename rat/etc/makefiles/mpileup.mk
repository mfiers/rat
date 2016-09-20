-include ./config.mk
-include ./rat_core.mk

### Input BAM directory
BAM_DIR ?= ../00.fastq
BAM_WILDCARD ?= *
BAM_EXT ?= .bam
DBSNP_VCF ?= ~/vib11/resources/human/hg19/dbSNP/All.vcf.gz
GENOME_FASTA ?=  /user/leuven/306/vsc30690/vib11/resources/human/hg19.fa
SNPEFF_GENOME ?= hg19
BCFTOOLS_MODULE ?= BCFtools/1.3-foss-2014a
SAMTOOLS_MODULE ?= SAMtools/1.3-foss-2014a
SNPEFF_MODULE ?= SnpEff/4.2
BAM_FILES = $(wildcard $(BAM_DIR)/$(BAM_WILDCARD)$(BAM_EXT))

SHELL:=/bin/bash

.PHONY: run prepare mpileup

all: prepare mpileup.snpeff.vcf

prepare:
	@echo "bam input: $(BAM_DIR)/$(BAM_WILDCARD)$(BAM_EXT)"
	@echo "Discovered $(words $(BAM_FILES)) bam files"


mpileup.snpeff.vcf: mpileup.vcf
	module load $(SNPEFF_MODULE) \
		&& snpEff ann $(SNPEFF_GENOME) \
			 -csvStats mpileup.sneff.stats.vcf  \
			 > $@ < $<

mpileup.vcf: mpileup.raw
	module load $(SNPEFF_MODULE) \
		&& snpSift annotate -tabix $(DBSNP_VCF) - \
			< $< > $@

mpileup.raw.2: mpileup.raw
	module load $(BCFTOOLS_MODULE) \
		&& module load $(SAMTOOLS_MODULE) \
		&& bcftools filter -s LowQual -e '%QUAL<20 || DP>100' $< \
			> $@

mpileup.raw: $(BAM_FILES)
	module load $(BCFTOOLS_MODULE) && \
		&& module load $(SAMTOOLS_MODULE) \
		&& samtools mpileup -q 5 -u -t 'DP,AD,INFO/AD,ADF,ADR,SP' \
			--max-depth 1000 \
			-f $(GENOME_FASTA) \
			$^ \
			| bcftools call -mv \
			> mpileup.raw
