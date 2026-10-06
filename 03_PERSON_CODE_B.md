# Участник 3 — код B: forecasts и результаты статьи

## Ответственность

Полностью закрыть прогнозную и эмпирическую часть статьи:

- optimal forecast simulations Sections 4.1–4.2;
- fBm, bfBm, mfBm3, mfBm4, mfBm5;
- rolling-window evaluation;
- HAR и все vector/common-factor/log-HAR варианты статьи;
- DJ30 empirical forecasting;
- Magnificent 7 и другие robustness experiments;
- итоговые forecast tables и figures.

## Дедлайны

- **3 октября, 22:00** — все forecasting scripts запущены в smoke/reduced режиме.
- **4 октября, 20:00** — рабочий pipeline создаёт fBm/mfBm/HAR outputs.
- **5 октября, 12:00** — закончены самопроверка, тесты, robustness и сравнение с оригиналом.
- **5 октября, 15:00 — ваш основной дедлайн** — review-ready PR передан участнику 5.
- **6 октября, 10:00** — получить первый review.
- **6 октября, 15:00 — дедлайн исправлений** — исправлены blocker/major замечания.
- **6 октября, 19:00** — получить финальный статус reviewer.
- **6 октября, 21:00** — объяснить докладчику принятые forecast graphs и metrics.

После 5 октября, 15:00 — только исправления и закрытие замечаний reviewer.

## Зависимость

Использовать данные и функции estimators участника 2. Если они ещё не готовы, начать с mock/synthetic input и заменить его без переписывания forecasting pipeline.

## Куда загружать

```text
R/forecast_mfbm.R
R/forecast_har.R
R/rolling_window.R
R/metrics.R
scripts/05_run_forecast_simulations.R
scripts/06_run_empirical_forecasts.R
scripts/07_run_robustness.R
results/tables/code_b_*.csv
results/figures/code_b_*.pdf
results/logs/code_b.log
tests/testthat/test-metrics.R
tests/testthat/test-rolling-window.R
```

## Обязательные outputs

- theoretical forecast simulation results;
- fBm/bfBm/mfBm3–5 comparison;
- HAR-family comparison;
- RMSFE tables;
- empirical forecast figures;
- robustness results;
- краткая таблица `original vs ours` по своему блоку.

## Как проверить самостоятельно

- Все модели оцениваются на одинаковых test dates.
- Одинаковы target, forecast horizon и metric.
- `max(training_date) < forecast_date` для каждого прогноза.
- RMSFE проверен вручную на маленьком примере.
- Forecast и actual имеют одинаковую длину.
- Нет неожиданных `NA/Inf`.
- Все модели запускаются одной функцией/единым pipeline.
- Все результаты сохранены в `results/`.
- Все свои таблицы и фигуры отметить в coverage matrix.

## Как сдавать reviewer

В PR указать:

```text
Команды запуска:
Модели и периоды:
Созданные outputs:
Runtime:
Что совпало с оригиналом:
Что отличается и почему:
Проверка look-ahead:
```

Для передачи докладчику подготовить:

- объяснение RMSFE и остальных метрик;
- расшифровку каждого финального forecast graph;
- список моделей от простейшей к наиболее сложной;
- два предложения о том, почему mfBm выигрывает или не выигрывает.
