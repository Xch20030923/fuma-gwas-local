# FUMA local project contract

## Scope

The fixed production profile is FUMA v2.1.6-compatible, GRCh37/hg19, 1000 Genomes Phase 3 EUR, 503 samples, Ensembl v102, dbSNP146-labelled input/reference policy, CADD v1.4, RegulomeDB v1.1, Roadmap 127-state ChromHMM, and GWAS Catalog e110. The official vendor commit recorded by the project is `9ccff60570ea06a43bd7fa77aeb62920ad271df4`.

The local result is deliberately classified as a practical candidate reproduction. The FUMA private dbSNP146/RsMerge146 processed snapshot, private EUR precomputed LD/frequency archive, and private processed annotation snapshots are not all publicly recoverable. The production workflow therefore preserves the strongest reproducible result and records these limits in provenance.

## Screenshot/output contract

The standard FUMA output directory is expected to contain the following files. The local worker keeps these names unchanged.

| File | Role | Fidelity status for the audited AF-anxious same-input comparison |
|---|---|---|
| `IndSigSNPs.txt` | independent significant SNPs | exact rows, fields, and key set |
| `leadSNPs.txt` | lead SNP groups | exact rows, fields, and key set |
| `GenomicRiskLoci.txt` | merged genomic risk loci | exact rows, fields, and key set |
| `genes.txt` | V2G gene table | exact file in the audited comparison |
| `snps.txt` | candidate SNP annotation | key set exact; snapshot-dependent field differences in 977 rows |
| `ld.txt` | candidate LD rows | key set exact; 3 r2 differences at one multiallelic coordinate due private LD-row selection |
| `annov.txt` | ANNOVAR consequence table | candidate membership and annotation; 12 current-only rows and 48 field differences versus historical snapshot |
| `annot.txt` | functional annotation matrix | key set exact; 327 field differences, mainly CADD formatting/value and RegulomeDB snapshot effects |
| `annov.stats.txt` | ANNOVAR summary | produced when candidate background finalization completes |
| `EUR.annov.count` | EUR ANNOVAR background counts | produced when candidate background finalization completes |
| `gwascatalog.txt` | GWAS Catalog overlap | produced when the catalog resource is enabled |
| `params.config`, `README` | run parameters and readme | preserved from the local run |

The exact core result for the audit input was 19 independent significant SNPs, 13 lead groups, and 9 loci. The candidate set was 1,068 SNPs, with 20 genes in `genes.txt`. These numbers are the key practical result for that input, not a promise for every future GWAS.

## Input contract

The detector accepts a single file or a directory. For each candidate it records the source hash and selected member if the source is a zip archive. It recognizes:

- chromosome: `chr`, `CHR`, `chrom`, `chromosome`
- position: `pos`, `BP`, `bp`, `position`, `base_pair_location`
- variant: `SNP`, `rsID`, `MarkerName`, `variant_id`, or coordinate IDs such as `2:179464954`
- p value: `P`, `pval`, `pvalue`, `p.placo`, `p_placo`, `PLACO_P`, and common meta-analysis names

The canonical file written for the production launcher is tab-delimited with exactly `chr`, `pos`, `SNP`, and `p.placo`. Chromosomes `chr1`-`chr22`, `23`, and `X` are normalized. `P=0` is clamped to `1e-300` and recorded. Invalid rows are dropped with reason counts. Duplicate variant rows retain the smallest p value and the number removed is recorded.

The required biological interpretation is GRCh37/hg19. If a source is GRCh38, liftover is not silently attempted; the job must be stopped for explicit conversion or marked as incompatible.

## Runtime and parallelism

Runtime depends mainly on input size, the number of significant candidates, and whether the shared reference/background checkpoints are warm.

- Detection and normalization: seconds to minutes; an 8.9-million-row compressed input is expected to take minutes rather than hours.
- Warm single job: use approximately 10-30 minutes as an initial planning range for the local production workflow after references and ANNOVAR background checkpoints are ready. The exact job timestamps are authoritative.
- First-time reference/background preparation: hours to days. It is a one-time server preparation and should not be confused with each GWAS's analytical runtime.
- Recommended parallelism: two jobs. Three or four are supported by the queue but can increase memory pressure and PLINK contention. The queue never starts more than four.

The queue is detached with a scheduler PID and log. User or Codex suspension only stops observation; it does not intentionally stop the detached worker or its child production processes. Later status is read from `/media/desk16/iy19619/iyun8003/post-GWAS分析/FUMA_local_runs/jobs/*.json`.

## Downstream post-GWAS interface

For each completed run, use the generated interface directory as the `path` argument for the existing R function `process_fuma_coloc_annov_all` in:

`/media/desk16/iy19619/iyun8003/post-GWAS分析/post-GWAS代码文件/【241117】【重要代码保存计划】/【完美】FUMA结果一键出所有的结果.R`

The directory contains symlinked raw files and these generated files:

- `GenomicRiskLoci_with_NearestGene_and_Coloc.csv`
- `GenomicRiskLoci_with_NearestGene_and_Coloc.xlsx`
- `<run_id>_Final.xlsx` with `Nearest_V2G_Exonic` and `GenomicRiskLoci+coloc` sheets
- `FUMA_V2G_Nearest_Exonic.csv`
- `FUMA_postGWAS_manifest.json`

The generated locus table includes `LeadSNP_NearestGene`, matching the existing R code. `cytoBand` is initially `NA` with a provenance note; the R workflow can fill it using `QTLMR::data_info_cytoBand(build=37)`. A coloc file is optional and is joined by `risk_locus` or `GenomicLocus`.

## Current reference/download decision

All reference classes needed by the production profile have already been downloaded or promoted to the local audited candidate/production bundle. No new large download is required before submitting ordinary future GWAS jobs. The remaining non-exact classes are an availability/reproduction limitation, not a reason to block key locus results.
