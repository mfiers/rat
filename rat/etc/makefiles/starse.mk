-include ./config.mk
-include ./rat_core.mk

FASTQ_DIR ?= ../00.fastq
FASTQ_EXT ?= .fastq.gz
STAR_DB ?= ~/vib11/resources/mouse/ensembl-84/STAR/
STAR_MODULE ?= STAR/2.5.1b-foss-2014a
SAMTOOLS_MODULE ?= SAMtools/1.3-foss-2014a
SUBREAD_MODULE ?= Subread/1.4.6-p2-foss-2014a
GTF_FILE ?= ~/vib11/resources/mouse/ensembl-84/Mus_musculus.GRCm38.84.chr.gtf

FASTQ_FILES = $(wildcard $(FASTQ_DIR)/*$(FASTQ_EXT))
BASENAMES = $(subst $(FASTQ_EXT),,$(notdir $(FASTQ_FILES)))

SHELL:=/bin/bash

BAM_FILES = $(addprefix out/, $(addsuffix _Aligned.sortedByCoord.out.bam, $(BASENAMES)))
BAI_FILES = $(subst .out.bam,.out.bam.bai,$(BAM_FILES))
STATS_FILES = $(subst .out.bam,.out.bam.stats,$(BAM_FILES))
IDX_FILES = $(subst .out.bam,.out.bam.idx,$(BAM_FILES))

all: star_prepare star_run samtools_index samtools_stats samtools_idx process_counts

star_prepare:
	@echo -e "\n\n# STAR RUN $(date)" >> mad.sh
	@echo "module load $(SAMTOOLS_MODULE)" >> mad.sh
	@echo -e "module load $(STAR_MODULE)\n\n" >> mad.sh
	@mkdir -p out
	@echo "Found $(words $(FASTQ_FILES)) file(s) in $(FASTQ_DIR)"


process_counts: counts.raw.tsv counts.group.tsv


counts.group.tsv counts.raw.tsv: counts/output.tsv
	./summarize.py
	k3 

counts/output.tsv: $(BAM_FILES)
	@mkdir -p counts
	module load $(SUBREAD_MODULE) \
	  && featureCounts \
			-a $(GTF_FILE) \
			-o counts/output.tsv \
			-Q 10 \
			-T 8 \
			-g gene_name \
			$^


samtools_idx: $(IDX_FILES)
samtools_stats: $(STATS_FILES)
samtools_index: $(BAI_FILES)

$(IDX_FILES): out/%.bam.idx: out/%.bam
	@echo "Run samtools idxstats for : " $<
	@echo "generating                : " $@
	$(eval $@_CL := samtools idxstats $< > $@ )

	@if [[ -s "$<" ]]; then \
		module load $(SAMTOOLS_MODULE) \
			&& $($@_CL); \
		echo -e "\n\n# $@" >> mad.sh; \
		echo "mad ta add --executable samtools \
				--input $< --output $@ \
				--cl '$($@_CL)'" >> mad.sh; \
		echo 'rat samtools_idxstats $@ \
			| mad apply_from_table -p samtools_ - $<' \
			>> mad.sh; \
	fi


$(STATS_FILES): out/%.bam.stats: out/%.bam
	@echo "Run samtools stats for : " $<
	@echo "generating             : " $@

	$(eval $@_CL := samtools stats $< > $@)

	module load $(SAMTOOLS_MODULE) && $($@_CL)

	@if [[ -s "$<" ]]; then \
		echo "mad ta add --executable samtools --input $< --output $@ --cl '$($@_CL)'" \
			>> mad.sh; \
		echo "rat samtools_stats $@ | mad apply_from_table -p samtools_ - $<" \
			>> mad.sh; \
	fi



filled = $(shell [[ -s "$(1)" ]] && echo $(1))


$(BAI_FILES): out/%.bam.bai: out/%.bam
	@echo "Run samtools index for : " $<
	@echo "Generating             : " $@
	@echo "input file contain data: " $(call filled,$<)
	$(eval $@_CL := samtools index $< $@)

	if [[ -s "$<" ]]; \
	then \
		module load $(SAMTOOLS_MODULE) && $($@_CL); \
		echo -e "\n\n# $@" >> mad.sh; \
		echo "mad ta add --executable samtools --input $< --output $@ \
			--cl '$($@_CL)'" >> mad.sh; \
	fi

star_run: $(BAM_FILES)

$(BAM_FILES): out/%_Aligned.sortedByCoord.out.bam: $(FASTQ_DIR)/%$(FASTQ_EXT)
	@echo "Running star for   : " $<
	@echo "Name               : " $*
	@echo "Output bam         : " $@
	@echo "STAR module        : " $(STAR_MODULE)
	$(eval $@_CL := \
	module load $(STAR_MODULE) \
	  && module load $(SAMTOOLS_MODULE) \
		&& $($@_CL)

	echo -e "\n\n# $@" >> mad.sh
	echo "mad ta add --executable STAR \
			--cl '$($@_CL)' \
			--db $(STAR_DB)/SA \
			--db $(STAR_DB)/SAindex \
			--db $(STAR_DB)/Genome \
			--db $(STAR_DB)/exonInfo.tab \
			--db $(STAR_DB)/geneInfo.tab \
			--db $(STAR_DB)/exonGeTrInfo.tab \
			--db $(STAR_DB)/sjdbInfo.txt \
			--db $(STAR_DB)/transcriptInfo.tab \
			--db $(STAR_DB)/sjdbList.fromGTF.out.tab \
			--db $(STAR_DB)/sjdbList.out.tab \
			--input $<  \
			--output bam:$@ \
			--output junctions:out/$*_SJ.out.tab \
			--output chimeric:out/$*_Chimeric.out.junction \
			--output chimeric:out/$*_Chimeric.out.sam" \
				>> mad.sh
