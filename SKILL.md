---
name: fuma-gwas-local
description: Run the fixed local FUMA-compatible SNP2GENE workflow for GWAS, PLACO, or other summary-statistics files. Use when a user asks to recognize GWAS columns, prepare FUMA input, run FUMA SNP2GENE, process one or more GWAS pairs, compare key FUMA outputs, or connect FUMA results to the local post-GWAS R workflows. The skill uses the audited GRCh37/EUR/v2.1.6-compatible project, launches detached jobs, supports parallel queues, and preserves provenance for later sessions.
---

# FUMA GWAS Local

Use this skill as the default entry point for local FUMA/SNP2GENE work. It is intentionally tied to the server's audited project and post-GWAS directory so a later session can recover jobs and context without reconstructing paths.

## Fixed locations

- FUMA project: `/media/desk16/iy19619/FUMA-compatible-SNP2GENE`
- Production profile: `/media/desk16/iy19619/FUMA-compatible-SNP2GENE/config/profiles/fuma_v2.1.6_eur_hg19_ensembl_v102.yaml`
- Post-GWAS code: `/media/desk16/iy19619/iyun8003/post-GWAS分析/post-GWAS代码文件`
- Detached job root: `/media/desk16/iy19619/iyun8003/post-GWAS分析/FUMA_local_runs`
- Main script: `/media/desk16/iy19619/.codex/skills/fuma-gwas-local/scripts/fuma_gwas_local.py`

Read `references/project_contract.md` when explaining fidelity, runtime, required inputs, or downstream filenames. Do not search for a second FUMA installation unless the fixed project is unavailable.

## Operating rules

1. For a new GWAS file, run `detect` first when the file is unfamiliar; otherwise submit it directly. The detector accepts tab, comma, whitespace, gzip, and zip-wrapped text and recognizes common `chr`, `pos/bp`, `SNP/rsID`, `P/P_VALUE`, and PLACO column names.
2. Submit jobs with `submit`. It writes a canonical `chr`, `pos`, `SNP`, `p.placo` input, records source and normalized SHA256 values, and starts a detached queue worker. Do not keep a Codex turn open merely to watch a long-running job.
3. Use `workers=2` by default. Up to four jobs are supported, but production FUMA work is memory and PLINK-I/O intensive; raise parallelism only when the server has headroom.
4. The worker calls the existing production launcher with `FUMA_PROFILE`, `FUMA_INPUT`, and `FUMA_RUN_ID`. It never silently substitutes a different genome build, panel, database release, or historical snapshot.
5. After completion, the worker creates a downstream interface directory containing symlinked raw FUMA files plus `GenomicRiskLoci_with_NearestGene_and_Coloc.csv`, XLSX workbooks, `FUMA_V2G_Nearest_Exonic.csv`, and `FUMA_postGWAS_manifest.json`.
6. Report the result as `practical_candidate` unless a genuine same-input and same-reference golden comparison proves strict identity. Core locus membership is the priority result; snapshot-dependent annotation fields must remain flagged rather than presented as exact.
7. Later sessions should use `status`, `report`, and the job JSON under the fixed job root. A detached worker writes independent logs and survives the end of the current Codex session.

## Commands

```bash
python3 /media/desk16/iy19619/.codex/skills/fuma-gwas-local/scripts/fuma_gwas_local.py detect /path/to/gwas.tsv --output /path/to/detection.json
python3 /media/desk16/iy19619/.codex/skills/fuma-gwas-local/scripts/fuma_gwas_local.py submit /path/to/gwas.tsv --workers 2
python3 /media/desk16/iy19619/.codex/skills/fuma-gwas-local/scripts/fuma_gwas_local.py status
python3 /media/desk16/iy19619/.codex/skills/fuma-gwas-local/scripts/fuma_gwas_local.py status --run-id RUN_ID
python3 /media/desk16/iy19619/.codex/skills/fuma-gwas-local/scripts/fuma_gwas_local.py report RUN_ID
```

Pass a coloc summary to `submit --coloc FILE` or later to `interface RUN_ID --coloc FILE`. Coloc rows are matched by `risk_locus` or `GenomicLocus` and are written with a `coloc_` prefix when needed.

## Completion checks

A complete job has state `completed`, a non-empty `reproducibility_report.json`, all key files `IndSigSNPs.txt`, `leadSNPs.txt`, `GenomicRiskLoci.txt`, `snps.txt`, `ld.txt`, `genes.txt`, `annov.txt`, and `annot.txt`, and a passing `FUMA_postGWAS_manifest.json`. A failed production process or failed interface step must remain visible as `failed` or `postprocess_failed`.

Do not claim exact screenshot-level reproduction solely from a passing local run. Use the comparison matrix in the reference contract and the run's own provenance.
