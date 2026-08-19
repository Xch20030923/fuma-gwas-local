# fuma-gwas-local

Codex skill for the fixed local FUMA-compatible SNP2GENE workflow.

This repository intentionally contains only the skill package:

- `SKILL.md`: trigger, operational, and completion contract.
- `agents/openai.yaml`: Codex skill metadata.
- `scripts/fuma_gwas_local.py`: GWAS detector, detached queue runner, and post-GWAS interface generator.
- `references/project_contract.md`: audited reference and fidelity contract.

It intentionally excludes the local FUMA project, all reference databases, GWAS inputs, historical runs, logs, and post-GWAS code. Those resources must be transferred and installed separately on the target Linux server.

## Required target-side resources

The skill's current production contract expects these paths to exist on the target server, or for the skill/package to be adapted before use:

- FUMA project: `/media/desk16/iy19619/FUMA-compatible-SNP2GENE`
- Production profile: `config/profiles/fuma_v2.1.6_eur_hg19_ensembl_v102.yaml` inside that project
- Post-GWAS root: `/media/desk16/iy19619/iyun8003/post-GWAS分析`
- Post-GWAS R code: `/media/desk16/iy19619/iyun8003/post-GWAS分析/post-GWAS代码文件`

The exact reference bundle and its known non-exact snapshot limitations are documented in `references/project_contract.md`.
