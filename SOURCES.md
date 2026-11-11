# Sources

## Dataset

- **Accession**: GSE66494 (NCBI GEO)
- **Series page**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66494
- **Platform**: GPL6480 (Agilent-014850 Whole Human Genome 4x44K)
- **Content**: kidney biopsy gene-expression arrays, 61 samples (CKD and control)

## Files used

| File | URL | Downloaded / verified | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| GSE66494_RAW.tar | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66494/suppl/GSE66494_RAW.tar | 2026-11-11 | 137799680 | 7a9462e1da5fc1608ff6d85ce896a8dc9345232d4468e8d26ab7de39af8749b6 |
| GSE66494_series_matrix.txt.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66494/matrix/GSE66494_series_matrix.txt.gz | 2026-11-11 | 14014463 | 6772d256790ffd008755eac4d39fd2db9e47266c9b586f93af9b999a513b95d9 |
| GPL6480.annot.gz | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz | 2026-11-11 | 7380104 | aacec81b2d95f183f0827b43e769388a92917811b81920c13ecc8fe148c5a3ab |

## Notes

- The data files are not stored in the repository. `run.sh` downloads them from the URLs above and verifies them with `sha256sum -c checksums.sha256`.
- The tar (copied into the project on 2026-11-07) and the series matrix were first obtained by hand, then checked on 2026-11-11 to be byte-identical to a fresh download from the URLs above (same SHA-256). The pipeline itself downloads them with `run.sh`.
- GEO data are public. See https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html

## Associated publication

To be completed from the "Citation" field of the GEO series page (see the series URL above).
