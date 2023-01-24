module load bedtools
rbc=/oak/stanford/groups/smontgom/amarder/t21_multiome/output/finemap/finemappedtraits_hg19/hg38/rbc.PP001.hg38.bed
bedtools intersect -a $rbc -b ~/tmp/peaks.bed -wa -wb > ~/tmp/rbc_peaks_overlap.bed