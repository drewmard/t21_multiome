# t21_multiome

RNA + ATAC in T21 and healthy fetal liver samples

## Notes

- I made a directory "scripts" with subdirectories "scripts/andrew" and "scripts/jon"
- I also made a directory "output" that is in .gitignore so any output here will not be pushed to github (useful for generating large data & saving)

## Paths

#### Healthy datasets
- **ATAC (Scanpy):** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.ATAC_only2.h.h5ad
- **RNA (Scanpy):** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_only2.h.h5ad
- **chromvar:** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/h.ChromVAR.txt

#### T21 datasets:
- **ATAC (Scanpy):** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.ATAC_only2.ds.h5ad
- **RNA (Scanpy):** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_only2.ds.h5ad
- **chromVAR:** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/ds.ChromVAR.txt

#### Other data: 
- **TF motif ID to TF name dataframe:** /oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/motif_gene.ChromVAR.txt