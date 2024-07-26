### install ABC-model as this tutorial: https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/usage/getting_started.html

### prepare the ATAC-seq data
awk '$1~/chr/' $i.tagAlign|  bedtools sort -faidx reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv -i -   |bgzip -c  > $i.tagAlign.gz 

###download average HiC data

##prepare the config files: config.hsc.HSCall.yaml config_biosamples_hsc.HSCsall.tsv


snakemake --configfile config.hsc.HSCall.yaml

