# What to commit on the C branch

Add **only `code_c/` at the team repository root**. `code_c_upload.zip` contains that folder, including code, tests, bilingual documentation and small recorded tables/figures. Do not commit the ZIP itself.

The merged reviewed A/B project must provide `R/code_b_config.R`, `R/forecast_mfbm.R`, `R/rolling_window.R` and `data/processed/dj30_logvol.csv`. These remain colleagues' files. Do not upload the entire local RS workspace: it contains older versions, branch snapshots, the article, temporary files and historical outputs.

## Steps

1. Switch to your C branch, based on the team's agreed current base. C can be committed before A/B are merged, but execution requires their reviewed interfaces.
2. Extract the archive so the path is `repository/code_c/`, without an extra `RS/` parent. Preserve the same structure for a GitHub web upload.
3. From the repository root, test, run the small workflow, stage C only, inspect the staged list and commit:

```sh
Rscript code_c/scripts/run_tests.R
Rscript code_c/scripts/run_experiments.R smoke all

git add -- code_c
git diff --cached --name-only
git diff --cached --check
git commit -m "Add Code C mfOU forecasting and research experiments"
```

Every staged path for this commit should start with `code_c/`. If unrelated changes were already staged, remove them from this commit's staging first. Branch names and push commands depend on your repository and are not invented here.

Adding this folder does not undo earlier commits on an old C branch. For minimal conflict risk, use a fresh branch from the agreed current team base and add only this overlay. If continuing an existing C branch, inspect its entire diff against the base: earlier shared-file edits still participate in the merge.

## Included files

- `R/`, `scripts/`, `tests/`: implementation and checks.
- `README*.md`, `EXPLAINED_*.md`, `RESULTS_*.md`, `COMMIT_GUIDE_*.md`, `docs/`: bilingual usage, explanations, actual outcomes and methods.
- `examples/`: small executed-research artifacts with provenance and limits; useful for review, not required for a fresh run.
- `.gitignore`, `FILE_MANIFEST.csv`: local exclusions and SHA-256 inventory of all other payload files; the manifest excludes itself.

Exclude local A/B snapshot folders, legacy root-level C files, `.Rlib`, `tmp/`, bulk `results/`, checkpoints and duplicate raw data. A/B already supply processed data. Generated C outputs go to `results/code_c_v2/`; staging only `code_c/` avoids committing them without editing the team's shared ignore file.

## Integration limits

C adds new paths without changing A/B, their configuration, root README or common entry point. Compatibility is verified against the supplied reviewed versions, including A/B located at the repository root. Snapshot fingerprints establish that colleagues' originals remain unchanged.

Unseen remote changes cannot be guaranteed conflict-free. Another branch adding `code_c/`, or an incompatible B API change, would need coordination. The local evidence covers the supplied versions rather than future commits.

RS has no `.git`, and no remote URL or branch name was supplied. No actual commit or push was made. The upload payload is ready for your branch; no artificial repository history was created.
