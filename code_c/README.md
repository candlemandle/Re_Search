# Code C — additional forecasting experiments

Start with [the simple explanation in Russian](EXPLAINED_RU.md) or [English](EXPLAINED_EN.md).

- [Русский README: запуск, данные, архитектура, результаты](README_RU.md)
- [English README: usage, data, architecture, outcomes](README_EN.md)
- [Выполненные результаты / executed results](RESULTS_RU.md) · [English](RESULTS_EN.md)
- [Что загрузить в свою ветку](COMMIT_GUIDE_RU.md) · [English](COMMIT_GUIDE_EN.md)
- [Формулы и основания экспериментов](docs/METHODS_RU.md) · [English](docs/METHODS_EN.md)

This folder is an additive extension for reviewed A/B v2. A/B provide the data, parameter utilities and authoritative trading-calendar mfBm baseline. Code C adds stationary mfOU forecasting, window sensitivity, controlled simulations and asset-selection experiments. It does not modify A/B files.

Run from the repository root after merging the reviewed A/B branches:

```sh
Rscript code_c/scripts/run_tests.R
Rscript code_c/scripts/run_experiments.R smoke all
```

Tests require `testthat`; experiments use base R and A/B's existing functions. The supplied recorded study is exploratory and does not establish a general winning model.
