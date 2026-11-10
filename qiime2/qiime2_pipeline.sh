#!/usr/bin/env bash
# =============================================================================
# qiime2_pipeline.sh — 16S rRNA (V4, 515F/806R) processing, QIIME 2 2024.10.1
#
# Parameters are taken from the provenance recorded in the output artifacts
# (qiime tools peek / action.yaml) and from the run logs.
# Input: demultiplexed paired-end FASTQ (2 x 300 bp, Illumina MiSeq), 205 samples.
# Raw reads: NCBI SRA BioProject PRJNA1459468.
# =============================================================================
set -euo pipefail

THREADS=8
OUT=reprocess_v2
mkdir -p "$OUT"

# 0. Import (paired-end FASTQ manifest; Phred33) -------------------------------
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path demux.qza

# 1. Primer removal: anchored 5' primers; pairs lacking either primer discarded
#    (98.8-99.2% of pairs retained per sample)
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f ^GTGYCAGCMGCCGCGGTAA \
  --p-front-r ^GGACTACNVGGGTWTCTAAT \
  --p-error-rate 0.1 \
  --p-overlap 3 \
  --p-discard-untrimmed \
  --p-cores "$THREADS" \
  --o-trimmed-sequences "$OUT/demux-trimmed.qza" \
  --verbose > "$OUT/cutadapt_log.txt" 2>&1

# 2. DADA2 (DADA2 1.30.0): truncation 230/200, maxEE 2, consensus chimera removal
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$OUT/demux-trimmed.qza" \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 200 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-min-overlap 12 \
  --p-pooling-method independent \
  --p-chimera-method consensus \
  --p-n-threads 0 \
  --o-table "$OUT/table.qza" \
  --o-representative-sequences "$OUT/rep-seqs.qza" \
  --o-denoising-stats "$OUT/denoising-stats.qza" \
  --verbose > "$OUT/dada2_log.txt" 2>&1

# 3. Length filter: keep ASVs of the expected V4 amplicon length (250-256 bp)
qiime tools export --input-path "$OUT/rep-seqs.qza" --output-path "$OUT/export_seqs"
awk '/^>/{id=substr($0,2); next} {len[id]+=length($0)}
     END{print "feature-id"; for (i in len) if (len[i]>=250 && len[i]<=256) print i}' \
  "$OUT/export_seqs/dna-sequences.fasta" > "$OUT/keep_len_ids.tsv"

qiime feature-table filter-seqs \
  --i-data "$OUT/rep-seqs.qza" \
  --m-metadata-file "$OUT/keep_len_ids.tsv" \
  --o-filtered-data "$OUT/rep-seqs-len.qza"

qiime feature-table filter-features \
  --i-table "$OUT/table.qza" \
  --m-metadata-file "$OUT/keep_len_ids.tsv" \
  --o-filtered-table "$OUT/table-len.qza"

# 4. Taxonomy: naive Bayes, SILVA 138 99% OTUs, human-stool-weighted
#    (pre-trained weighted classifier, q2-clawback / readytowear collection)
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-human-stool-weighted-classifier.qza \
  --i-reads "$OUT/rep-seqs-len.qza" \
  --p-confidence 0.7 \
  --o-classification "$OUT/taxonomy.qza"

# 5. Remove ASVs without a phylum assignment and mitochondria/chloroplast/eukaryota
#    (archaeal ASVs, e.g. Methanobrevibacter, are retained)
qiime taxa filter-table \
  --i-table "$OUT/table-len.qza" \
  --i-taxonomy "$OUT/taxonomy.qza" \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --p-mode contains \
  --o-filtered-table "$OUT/table-filt.qza"

qiime feature-table filter-seqs \
  --i-data "$OUT/rep-seqs-len.qza" \
  --i-table "$OUT/table-filt.qza" \
  --o-filtered-data "$OUT/rep-seqs-filt.qza"
#    -> 3,868 ASVs; 6,666,482 reads (median 28,925; range 13,241-174,297)

# 6. Phylogeny (MAFFT + FastTree, midpoint-rooted)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$OUT/rep-seqs-filt.qza" \
  --p-mask-max-gap-frequency 1.0 \
  --p-mask-min-conservation 0.4 \
  --o-alignment "$OUT/aligned-rep-seqs.qza" \
  --o-masked-alignment "$OUT/masked-aligned-rep-seqs.qza" \
  --o-tree "$OUT/unrooted-tree.qza" \
  --o-rooted-tree "$OUT/rooted-tree.qza"

# 7. Core diversity metrics, rarefied to 13,000 reads (retains all 205 samples)
#    metadata.tsv: SampleID + Status only (clinical metadata not distributed)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$OUT/rooted-tree.qza" \
  --i-table "$OUT/table-filt.qza" \
  --p-sampling-depth 13000 \
  --m-metadata-file metadata.tsv \
  --output-dir "$OUT/core_metrics"

qiime diversity alpha \
  --i-table "$OUT/core_metrics/rarefied_table.qza" \
  --p-metric chao1 \
  --o-alpha-diversity "$OUT/core_metrics/chao1_vector.qza"

# 8. Genus- and family-level tables (unrarefied) for differential abundance
qiime taxa collapse --i-table "$OUT/table-filt.qza" --i-taxonomy "$OUT/taxonomy.qza" \
  --p-level 6 --o-collapsed-table "$OUT/table-l6.qza"
qiime taxa collapse --i-table "$OUT/table-filt.qza" --i-taxonomy "$OUT/taxonomy.qza" \
  --p-level 5 --o-collapsed-table "$OUT/table-l5.qza"

# Downstream statistics (R) are in ../Reviewers/.
