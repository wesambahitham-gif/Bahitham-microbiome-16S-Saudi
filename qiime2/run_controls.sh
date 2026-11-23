#!/usr/bin/env bash
# =============================================================================
# run_controls.sh — process the sequencing-run controls with the same QIIME 2
# pipeline and parameters as the study samples (qiime2_pipeline.sh).
#
# Controls on the MiSeq run: 4 no-template PCR controls (NEG1–NEG4) and
# 4 positive controls (POS1–POS4; ZymoBIOMICS Microbial Community Standard, D6300).
# Usage: put the control FASTQs (*_R1_001.fastq.gz / *_R2_001.fastq.gz) and the
# classifier below in one folder, then run this script from that folder.
# Classifier: https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-human-stool-weighted-classifier.qza
# Outputs: controls_denoising_stats.tsv, controls_genus_counts.tsv
# =============================================================================
set -euo pipefail
THREADS=8
O=out; mkdir -p "$O"

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
for r1 in *_R1_001.fastq.gz; do s=${r1%%_*}; echo -e "$s\t$PWD/$r1\t$PWD/${r1/_R1_/_R2_}" >> manifest.tsv; done

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv --input-format PairedEndFastqManifestPhred33V2 --output-path $O/demux.qza
qiime cutadapt trim-paired --i-demultiplexed-sequences $O/demux.qza \
  --p-front-f ^GTGYCAGCMGCCGCGGTAA --p-front-r ^GGACTACNVGGGTWTCTAAT \
  --p-error-rate 0.1 --p-overlap 3 --p-discard-untrimmed --p-cores $THREADS \
  --o-trimmed-sequences $O/trimmed.qza
qiime dada2 denoise-paired --i-demultiplexed-seqs $O/trimmed.qza \
  --p-trunc-len-f 230 --p-trunc-len-r 200 --p-max-ee-f 2 --p-max-ee-r 2 --p-trunc-q 2 \
  --p-min-overlap 12 --p-pooling-method independent --p-chimera-method consensus --p-n-threads $THREADS \
  --o-table $O/table.qza --o-representative-sequences $O/rep.qza --o-denoising-stats $O/stats.qza
qiime tools export --input-path $O/rep.qza --output-path $O/seqs
awk '/^>/{id=substr($0,2);next}{len[id]+=length($0)}END{print "feature-id";for(i in len)if(len[i]>=250&&len[i]<=256)print i}' \
  $O/seqs/dna-sequences.fasta > $O/keep.tsv
qiime feature-table filter-features --i-table $O/table.qza --m-metadata-file $O/keep.tsv --o-filtered-table $O/table-len.qza
qiime feature-table filter-seqs --i-data $O/rep.qza --m-metadata-file $O/keep.tsv --o-filtered-data $O/rep-len.qza
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-nb-human-stool-weighted-classifier.qza \
  --i-reads $O/rep-len.qza --p-confidence 0.7 --p-n-jobs $THREADS --o-classification $O/taxonomy.qza
qiime taxa filter-table --i-table $O/table-len.qza --i-taxonomy $O/taxonomy.qza \
  --p-include p__ --p-exclude mitochondria,chloroplast,eukaryota --p-mode contains --o-filtered-table $O/table-filt.qza
qiime taxa collapse --i-table $O/table-filt.qza --i-taxonomy $O/taxonomy.qza --p-level 6 --o-collapsed-table $O/l6.qza

qiime tools export --input-path $O/stats.qza --output-path $O/stats
qiime tools export --input-path $O/l6.qza --output-path $O/l6
biom convert -i $O/l6/feature-table.biom -o controls_genus_counts.tsv --to-tsv
cp $O/stats/stats.tsv controls_denoising_stats.tsv
# Results (Supplementary Methods): positive controls 360–55,968 reads, 100% assigned to the
# eight D6300 bacterial taxa; no-template controls 1,337–2,135 reads (fecal-like carry-over).
