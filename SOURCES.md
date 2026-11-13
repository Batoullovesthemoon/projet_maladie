# Sources

## The dataset

- **Accession:** GSE66494 (NCBI GEO)
- **Series page:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66494
- **Platform:** GPL6480 (Agilent-014850 Whole Human Genome 4x44K)
- **What it contains:** gene-expression arrays from kidney biopsies, 61 samples in total
  (53 CKD biopsies and 8 controls)

## Files we use

| File | URL | Downloaded / verified | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| GSE66494_RAW.tar | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66494/suppl/GSE66494_RAW.tar | 2026-11-11 | 137799680 | 7a9462e1da5fc1608ff6d85ce896a8dc9345232d4468e8d26ab7de39af8749b6 |
| GSE66494_series_matrix.txt.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66494/matrix/GSE66494_series_matrix.txt.gz | 2026-11-11 | 14014463 | 6772d256790ffd008755eac4d39fd2db9e47266c9b586f93af9b999a513b95d9 |
| GPL6480.annot.gz | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz | 2026-11-11 | 7380104 | aacec81b2d95f183f0827b43e769388a92917811b81920c13ecc8fe148c5a3ab |

Together the three files weigh about 160 MB.

## Good to know

- The data files are **not** stored in the repository. `run.sh` downloads them from the
  URLs above and checks them with `sha256sum -c checksums.sha256`.
- At the start of the project we copied the tar (on 2026-11-07) and the series matrix
  by hand. On 2026-11-11 we downloaded both again from the URLs above and confirmed that
  they have exactly the same SHA-256 as our copies. The pipeline itself always fetches
  them with `run.sh`.
- GEO data are public. The usage conditions are described at
  https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html

## Publication linked to the dataset

Nakagawa S, Nishihara K, Miyata H, Shinke H, Tomita E, Kajiwara M, et al.
**Molecular Markers of Tubulointerstitial Fibrosis and Tubular Cell Damage in Patients
with Chronic Kidney Disease.** *PLoS ONE* 10(8): e0136994 (2015).
https://doi.org/10.1371/journal.pone.0136994
