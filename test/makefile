SHELL := /bin/bash
R = test_ref
O = test_out# out bin
D = test_diffs
READS_FA = $R/reads.fa
GENOME = $R/genome.fa
ANNOTATION = $R/annotation.gtf
ANNOTATION_INCOMPLETE = $R/annotation.incomplete.gtf
FLAIR = ../flair.py
FLAIR_BIN = ../bin/

UNCORRECTED_READS_BED = $R/test.align.bed
READS_BED = $R/test.corrected.bed

ISOFORMS_FA = $R/test.collapse.isoforms.fa
ISOFORMS_BED = $R/test.collapse.isoforms.bed
READS_MANIFEST = $R/reads_manifest.txt
COUNTS_MATRIX = $R/counts_matrix.2.tsv

SALMON_PATH = ~/bl/bin/salmon-0.12.0_linux_x86_64/bin/salmon
FIRSTPASS_Q_SAM = hello/tmp5_ls3qcw.firstpass.q.sam
FIRSTPASS_BED = $R/test.collapse.firstpass.bed

MMLOC = /cluster/home/jeltje/.conda/envs/mappy_env/bin/minimap2
SAMTOOLSLOC = ~/.conda/envs/mappy_env/bin/samtools

test-align-base :
	python $(FLAIR) align -r $(READS_FA) -g $(GENOME) -t 4 -o $O/test.align --psl #> $O/test.align.stdout 2> $O/test.align.stderr
	diff -q <(sort $R/test.align.bed) <(sort $O/test.align.bed) > $D/test.align.bed.diff
	diff -q <(sort $R/test.align.psl) <(sort $O/test.align.psl) > $D/test.align.psl.diff

test-align-params:
	python $(FLAIR) align -r $(READS_FA) --minimap2 $(MMLOC) --samtools $(SAMTOOLSLOC) --mm_index $R/genome.mmi --chromsizes $R/chrom.sizes --nvrna --quality 2 -N 1 -t 4 -o $O/test.align2
	diff <(sort $R/test.align2.bed) <(sort $O/test.align2.bed) > $D/test.align2.diff
# TODO: add a --junction_bed file (don't know how to make for this set yet)


test-correct-short-reads :
	python $(FLAIR) correct -q $(UNCORRECTED_READS_BED) -j $R/shortread_junctions.tab -g $(GENOME) -o $O/test.correct.shortread
	diff -q <(sort $R/test.correct.shortread_all_inconsistent.bed ) <(sort $O/test.correct.shortread_all_inconsistent.bed) > $D/test.correct_shortread.inconsistent.diff
	diff -q <(sort $R/test.correct.shortread_all_corrected.bed ) <(sort $O/test.correct.shortread_all_corrected.bed) > $D/test.correct_shortread.corrected.diff

test-correct-annotation :
	python $(FLAIR) correct -q $(UNCORRECTED_READS_BED) -f $(ANNOTATION_INCOMPLETE) -g $(GENOME) -o $O/test.correct.annotation
	diff -q <(sort $R/test.correct.annotation_all_inconsistent.bed) <(sort $O/test.correct.annotation_all_inconsistent.bed) > $D/test.correct_annotation.inconsistent.diff
	diff -q <(sort $R/test.correct.annotation_all_corrected.bed) <(sort $O/test.correct.annotation_all_corrected.bed) > $D/test.correct_annotation.corrected.diff

test-correct-short-reads-annotation :
	python $(FLAIR) correct -q $(UNCORRECTED_READS_BED) -j $R/shortread_junctions.tab  -f $(ANNOTATION_INCOMPLETE) -g $(GENOME) -o $O/test.correct.shortread.annotation
	diff -q <(sort $R/test.correct.shortread.annotation_all_inconsistent.bed) <(sort $O/test.correct.shortread.annotation_all_inconsistent.bed) > $D/test.correct_shortread_annotation.inconsistent.diff
	diff -q <(sort $R/test.correct.shortread.annotation_all_corrected.bed) <(sort $O/test.correct.shortread.annotation_all_corrected.bed) > $D/test.correct_shortread_annotation.corrected.diff

test-correct-params :
	python $(FLAIR) correct --query $(UNCORRECTED_READS_BED) --shortread $R/shortread_junctions.tab --genome $(GENOME) --output $O/test.correct.params --chromsizes $R/chrom.sizes --nvrna --ss_window 8
	diff -q <(sort $R/test.correct.params_all_inconsistent.bed ) <(sort $O/test.correct.params_all_inconsistent.bed) > $D/test.correct.params.diff

test-correct-window-25 :
	python $(FLAIR) correct -q $(UNCORRECTED_READS_BED) -j $R/shortread_junctions.tab  -f $(ANNOTATION) -g $(GENOME) -w 25 -o $O/test.correct.w25
	diff -q <(sort $R/test.correct.w25_all_corrected.bed) <(sort $O/test.correct.w25_all_corrected.bed) > $D/test.correct_w25.diff

# COUNT SAM

test-count-sam-transcripts :
	python $(FLAIR_BIN)count_sam_transcripts.py -s $(FIRSTPASS_Q_SAM) -t 4 --generate_map $O/test.count_sam_transcripts.map.txt -o $O/test.count_sam_transcripts 
	diff -q <(sort $R/test.count_sam_transcripts) <(sort $O/test.count_sam_transcripts) > $D/test.count_sam_transcripts.diff
	diff -q <(sort $R/test.count_sam_transcripts.map.txt) <(sort $O/test.count_sam_transcripts.map.txt) > $D/test.count_sam_transcripts.map.diff

test-count-sam-transcripts-trust-ends :
	python $(FLAIR_BIN)count_sam_transcripts.py -s $(FIRSTPASS_Q_SAM) -t 4 --generate_map $O/test.count_sam_transcripts.trust_ends.map.txt -o $O/test.count_sam_transcripts.trust_ends --trust_ends
	diff -q <(sort $R/test.count_sam_transcripts.trust_ends.map.txt) <(sort $O/test.count_sam_transcripts.trust_ends.map.txt) > $D/test.count_sam_transcripts.trust_ends.diff

test-count-sam-transcripts-stringent :
	python $(FLAIR_BIN)count_sam_transcripts.py -s $(FIRSTPASS_Q_SAM) -t 4 --generate_map $O/test.count_sam_transcripts.stringent.map.txt -o $O/test.count_sam_transcripts.stringent --stringent -i $(FIRSTPASS_BED)
	diff -q <(sort $R/test.count_sam_transcripts.stringent.map.txt) <(sort $O/test.count_sam_transcripts.stringent.map.txt) > $D/test.count_sam_transcripts.stringent.diff

# COLLAPSE

test-collapse-isoforms-precise-longest :
	python $(FLAIR_BIN)collapse_isoforms_precise.py -q $(READS_BED) -w 40 -n longest -o $O/test.collapse_isoforms_precise.longest.bed
	diff -q <(sort $R/test.collapse_isoforms_precise.longest.bed) <(sort $O/test.collapse_isoforms_precise.longest.bed) > $D/test.collapse_isoforms_precise.longest.diff

test-collapse-isoforms-precise-best-only :
	python $(FLAIR_BIN)collapse_isoforms_precise.py -q $(READS_BED) -w 40 -n best_only -o $O/test.collapse_isoforms_precise.best_only.bed
	diff -q <(sort $R/test.collapse_isoforms_precise.best_only.bed) <(sort $O/test.collapse_isoforms_precise.best_only.bed) > $D/test.collapse_isoforms_precise.best_only.diff

test-collapse-isoforms-precise-m2 :
	python $(FLAIR_BIN)collapse_isoforms_precise.py -q $(READS_BED) -w 40 -m 2 -s 0.05 -o $O/test.collapse_isoforms_precise.max2.bed
	diff -q <(sort $R/test.collapse_isoforms_precise.max2.bed) <(sort $O/test.collapse_isoforms_precise.max2.bed) > $D/test.collapse_isoforms_precise.max2.diff

test-collapse-sensitive :
	python $(FLAIR) collapse -r $(READS_FA) -q $(READS_BED) -g $(GENOME) -t 4 --generate_map --temp_dir $O/temp_collapse --keep_intermediate -f $(ANNOTATION) -o $O/test.collapse --quiet
	diff -q <(sort $R/test.collapse.isoforms.bed) <(sort $O/test.collapse.isoforms.bed) > $D/test.collapse.sensitive.diff
	diff -q <(sort $R/test.collapse.firstpass.bed) <(sort $O/test.collapse.firstpass.bed) > $D/test.collapse.firstpass.bed.diff
	diff -q <(sort $R/test.collapse.firstpass.q.counts) <(sort $O/test.collapse.firstpass.q.counts) > $D/test.collapse.firstpass.q.counts.diff
	diff -q <(sort $R/test.collapse.firstpass.unfiltered.bed) <(sort $O/test.collapse.firstpass.unfiltered.bed) > $D/test.collapse.firstpass.unfiltered.bed.diff
	diff -q <(sort $R/test.collapse.isoform.read.map.txt) <(sort $O/test.collapse.isoform.read.map.txt) > $D/test.collapse.isoform.read.map.txt.diff
	diff -q <(sort $R/test.collapse.isoforms.bed) <(sort $O/test.collapse.isoforms.bed) > $D/test.collapse.isoforms.bed.diff
	diff -q <(sort $R/test.collapse.isoforms.gtf) <(sort $O/test.collapse.isoforms.gtf) > $D/test.collapse.isoforms.gtf.diff

test-collapse-stringent :
	python $(FLAIR) collapse -r $(READS_FA) -q $(READS_BED) -g $(GENOME) -t 4 --generate_map --temp_dir $O/temp_collapse --keep_intermediate -f $(ANNOTATION) --stringent -o $O/test.collapse.stringent --quiet
	diff -q <(sort $R/test.collapse.stringent.isoforms.bed) <(sort $O/test.collapse.stringent.isoforms.bed) > $D/test.collapse.stringent.diff

test-collapse-remove :
	python $(FLAIR) collapse -r $(READS_FA) -q $(READS_BED) -g $(GENOME) -t 4 --generate_map --temp_dir $O/temp_collapse_to_remove -f $(ANNOTATION) -o $O/test.collapse.remove --quiet
	diff -q <(sort $R/test.collapse.isoforms.bed) <(sort $O/test.collapse.remove.isoforms.bed) > $D/test.collapse.remove.diff

test-collapse-annotation-reliant-fasta :
	python $(FLAIR) collapse -r $(READS_FA) -q $(READS_BED) -g $(GENOME) -t 4 --generate_map --annotation_reliant $R/annotation.fa --temp_dir $O/temp_collapse --keep_intermediate -f $(ANNOTATION) --stringent -o $O/test.collapse.annotation.reliant.fasta --quiet
	diff -q <(sort $R/test.collapse.annotation.reliant.fasta.isoforms.bed) <(sort $O/test.collapse.annotation.reliant.fasta.isoforms.bed) > $D/test.collapse.annotation.reliant.fasta.diff

test-collapse-annotation-reliant-generate :
	python $(FLAIR) collapse -r $(READS_FA) -q $(READS_BED) -g $(GENOME) -t 4 --generate_map --annotation_reliant generate --temp_dir $O/temp_collapse --keep_intermediate -f $(ANNOTATION) --stringent -o $O/test.collapse.annotation.reliant.generate --quiet
	diff -q <(sort $R/test.collapse.annotation.reliant.generate.annotated_transcripts.bed) <(sort $O/test.collapse.annotation.reliant.generate.annotated_transcripts.bed) > $D/test.collapse.annotation.reliant.generate.annotated_transcripts.diff
	diff -q <(sort $R/test.collapse.annotation.reliant.generate.annotated_transcripts.supported.bed) <(sort $O/test.collapse.annotation.reliant.generate.annotated_transcripts.supported.bed) > $D/test.collapse.annotation.reliant.generate.annotated_transcripts.supported.diff
	diff -q <(sort $R/test.collapse.annotation.reliant.generate.unassigned.bed) <(sort $O/test.collapse.annotation.reliant.generate.unassigned.bed) > $D/test.collapse.annotation.reliant.generate.unassigned.diff
	diff -q <(sort $R/test.collapse.annotation.reliant.generate.annotated_transcripts.alignment.counts) <(sort $O/test.collapse.annotation.reliant.generate.annotated_transcripts.alignment.counts) > $D/test.collapse.annotation.reliant.generate.annotated_transcripts.alignment.counts.diff
	diff -q <(sort $R/test.collapse.annotation.reliant.generate.annotated_transcripts.isoform.read.map.txt) <(sort $O/test.collapse.annotation.reliant.generate.annotated_transcripts.isoform.read.map.txt) > $D/test.collapse.annotation.reliant.generate.annotated_transcripts.isoform.read.map.txt.diff
	diff -q <(sort $R/test.collapse.annotation.reliant.generate.combined.isoform.read.map.txt) <(sort $O/test.collapse.annotation.reliant.generate.combined.isoform.read.map.txt) > $D/test.collapse.annotation.reliant.generate.combined.isoform.read.map.txt.diff

# test promoter region collapse
test-collapse-promoter :
	python $(FLAIR) collapse -r $(READS_FA) -q $(READS_BED) -g $(GENOME) -p $R/promoter_regions.bed -t 4 --temp_dir $O/temp_collapse --keep_intermediate -f $(ANNOTATION) -o $O/test.collapse.promoter
	diff -q <(sort $R/test.collapse.promoter.promoter_supported.bed) <(sort $O/test.collapse.promoter.promoter_supported.bed) > $D/test.collapse.promoter.promoter_supported.diff

### END OF MAIN TESTS ####
# this is not called (and it fails)
test-collapse-range :
	python $(FLAIR) 3.5 -r $R/test.align.bam -q $(READS_BED) -g $(GENOME) -t 4 --temp_dir $O/temp_collapse_range --keep_intermediate -f $(ANNOTATION) -o $O/test.collapse.range --quiet
	diff -q <(sort $R/test.collapse.isoforms.bed) <(sort $O/test.collapse.range.isoforms.bed) > $D/test.collapse.range.diff

test-quantify-mapq :
	python $(FLAIR) quantify -r $(READS_MANIFEST) -i $(ISOFORMS_FA) --temp_dir $O/test.collapse.temp --tpm -o $O/test.quantify.mapq.counts_matrix.tsv
	diff -q <(sort $R/test.quantify.salmon.counts_matrix.tsv) <(sort $O/test.quantify.salmon.counts_matrix.tsv) > $D/test.quantify.salmon.diff

test-quantify-mapq-trust-ends :
	python $(FLAIR) quantify -r $(READS_MANIFEST) -i $(ISOFORMS_FA) --temp_dir $O/test.collapse.temp --tpm -o $O/test.quantify.mapq.trust_ends.counts_matrix.tsv
	diff -q <(sort $R/test.quantify.mapq.trust_eds.counts_matrix.tsv) <(sort $O/test.quantify.mapq.trust_ends.counts_matrix.tsv) > $D/test.quantify.mapq.trust_ends.diff

test-quantify-salmon :
	python $(FLAIR) quantify -r $(READS_MANIFEST) -i $(ISOFORMS_FA) --temp_dir $O/test.collapse.temp --salmon $(SALMON_PATH) -o $O/test.quantify.salmon.counts_matrix.tsv
	diff -q <(sort $R/test.quantify.salmon.counts_matrix.tsv) <(sort $O/test.quantify.salmon.counts_matrix.tsv) > $D/test.quantify.salmon.diff

test-diffexp :
	python $(FLAIR) diffExp -q $(COUNTS_MATRIX) -o $O/test.diffexp -e 1 -of
	diff -q <(sort $R/test.collapse.isoforms.bed) <(sort $O/test.123.isoforms.bed) > $D/test.diffexp.diff

test-diffexp-e300 :
	python $(FLAIR) diffExp -q $(COUNTS_MATRIX) -o $O/test.diffexp -e 300 -of
# 	diff -q <(sort $R/bed) <(sort $O/bed) > $D/test.diffexp.diff

test-diffsplice :
	python $(FLAIR) diffSplice -i $(ISOFORMS_BED) -q $(COUNTS_MATRIX) --test -o $O/test.diffsplice
	diff -q <(sort $R/test.diffsplice.alt3.A_v_B_drimseq2_results.tsv) <(sort $O/test.diffsplice.alt3.A_v_B_drimseq2_results.tsv) > $D/test.diffsplice.alt3.diff

test-predict-productivity : 
	python $(FLAIR_BIN)predictProductivity.py -i $(ISOFORMS_BED) -g $(ANNOTATION) -f $(GENOME) > $D/test.predict_productivity.bed --longestORF
	#diff -q <(sort $R/bed) <(sort $O/bed) > $D/test..diff

test-diff-iso-usage : 
	python $(FLAIR_BIN)diff_iso_usage.py $(COUNTS_MATRIX) 1 4 $O/test.diff_iso_usage
	#diff -q <(sort $R/bed) <(sort $O/bed) > $D/test..diff

test-plot-isoform-usage : 
	python $(FLAIR_BIN)plot_isoform_usage.py $(ISOFORMS_BED) $(COUNTS_MATRIX) ENSG00000133703.12 $O/test.plot_isoform_usage.kras
	#diff -q <(sort $R/bed) <(sort $O/bed) > $D/test..diff

test-123 :
	python $(FLAIR) 123 -o $O/test.123 -r $(READS_FA) -g $(GENOME) -f $(ANNOTATION) --keep_intermediate
	diff -q $R/test.collapse.isoforms.bed $O/test.123.isoforms.bed > $D/test.123.diff

test-12346 :
	python $(FLAIR) 12346  -o $O/test.12346 --temp_dir $O/temp_collapse -r $(READS_FA) -g $(GENOME) -f $(ANNOTATION) --reads_manifest $(READS_MANIFEST) --quiet
	diff -q $R/test.collapse.isoforms.bed $O/test.12346.isoforms.bed > $D/test.12346.diff

test-123.54 :
	python $(FLAIR) 123.54 -o $O/test.123.54 -r $(READS_FA) -g $(GENOME) -f $(ANNOTATION) --reads_manifest $(READS_MANIFEST) 
	diff -q $R/test.collapse.isoforms.bed $O/test.123.54.isoforms.bed > $D/test.123.54.diff

test-234 :
	python $(FLAIR) 234 -o $O/test.234 -r $(READS_FA) -g $(GENOME) -f $(ANNOTATION) --reads_manifest $(READS_MANIFEST) 
	diff -q <(sort $R/test.collapse.isoforms.bed) <(sort $O/test.234.isoforms.bed) > $D/test.234.diff

test-align : test-align-base test-align-params
test-correct : test-correct-short-reads test-correct-annotation test-correct-short-reads-annotation test-correct-window-25 test-correct-params
test-cip : test-collapse-isoforms-precise-longest test-collapse-isoforms-precise-best-only test-collapse-isoforms-precise-m2
test-cst : test-count-sam-transcripts test-count-sam-transcripts-stringent test-count-sam-transcripts-trust-ends
test-collapse : test-cip test-cst test-collapse-sensitive test-collapse-stringent test-collapse-remove test-collapse-promoter
test-quantify : test-quantify-mapq test-quantify-mapq-trust-ends test-quantify-salmon
test-standalone : test-predict-productivity test-diff-iso-usage test-plot-isoform-usage
test-consecutive : test-123 test-12346 test-123.54 test-234

test : test-align test-correct test-collapse #test-quantify test-diffexp test-diffsplice test-standalone

clean :
	rm -rf $O/* $D/*
