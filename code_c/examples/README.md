# Recorded research artifacts / сохранённые результаты

These are small artifacts from actual executions on 9 October 2026. They are included for review and presentation; the runtime does not read them. Reproduction commands and input units are in [README_EN](../README_EN.md) and [README_RU](../README_RU.md). Interpretation: [English](../RESULTS_EN.md), [русский](../RESULTS_RU.md).

| Folder/file | Contents / содержание |
|---|---|
| `dense_study/` | Daily 2019–2021, one/two assets, 500 rows; paired losses and 20/40/60-day bootstrap sensitivity / ежедневное сравнение |
| `window_study/` | Twelve origins, all five nested sets, 250/500/1000 rows; descriptive comparisons / чувствительность к окну и размерности |
| `simulation_study/` | Three generators, thirty replications each; losses, MCSE, parameter bias/RMSE and status counts / контролируемая симуляция |
| `selection_study/` | Three targets, twelve origins, fixed/correlation/H-gap rules and fBm/HAR / выбор акций |
| `verification/` | Passed C/A/B tests, baseline audit and snapshot identity evidence / проверки совместимости |
| `executed_runs.csv` | Actual scope, counts, dates and raw run signatures / реестр выполненных запусков |
| `full_workload_preflight.csv` | Planned full workload, not executed full results / оценка объёма полного запуска |

Forecast RDS files, checkpoints, libraries and full per-origin forecast CSVs are intentionally outside the upload payload. Running the documented commands regenerates them under `results/code_c_v2/`. `configuration.txt` records each saved study's exact settings. In the preflight table, fields ending `_MB` were computed with a 1024² denominator and express MiB.

The original forecasting `source_identity.csv` and current `analysis_source_identity.csv` are distinct. Daily forecasts were aligned to B's shared past-only eligibility calendar, preserving the original raw RDS locally; sixteen one-asset rows from one origin were removed. No retained forecast was refitted during that alignment. The window study's alignment removed zero rows. Independent precision audits subsequently refitted representative real windows and matched the saved mfOU predictions exactly.

Raw signatures therefore retain their original computation provenance; they should not be interpreted as fingerprints of a fresh full refit with every final reporting change. Pointwise intervals are exploratory, sparse studies have no intervals, and no full-paper all-window daily result is claimed.

Это небольшие файлы для ревью и защиты. Программа их не использует при новом запуске. Полные прогнозы и кэш остаются локально и создаются заново документированными командами. `period=full` в таблицах означает всю выполненную выборку конкретного запуска. Полный период статьи со всеми окнами ежедневно здесь не рассчитывался.

Исходный код прогнозирования и код последующего анализа имеют отдельные отпечатки. После исправления общего календаря исходные RDS сохранены локально; оставшиеся прогнозы не изменены. Последующий независимый пересчёт на реальных окнах подтвердил их совпадение. Интервалы исследовательские и точечные, без поправки за множество сравнений. Малое число симуляций отражено в MCSE.
