# Часть C: выполненное сравнение

Режим: study . Фактические даты целевых значений: 2019-01-03 to 2021-11-05 . Горизонт считается в торговых днях.
Запланировано дат прогноза: 12 ; строк прогнозов: 1728 ; ошибок mfOU: 0 .

## Основное окно 500 дней

| dimension | h | metric | n | gain_pct | uncertainty_status |
| --- | --- | --- | --- | --- | --- |
| 1 |  1 | MSFE | 12 |  5.503 | unavailable_sparse_or_insufficient_blocks |
| 1 |  1 | QLIKE | 12 |  4.264 | unavailable_sparse_or_insufficient_blocks |
| 1 | 10 | MSFE | 12 | 36.855 | unavailable_sparse_or_insufficient_blocks |
| 1 | 10 | QLIKE | 12 | 17.936 | unavailable_sparse_or_insufficient_blocks |
| 1 | 20 | MSFE | 12 | 15.488 | unavailable_sparse_or_insufficient_blocks |
| 1 | 20 | QLIKE | 12 | 14.544 | unavailable_sparse_or_insufficient_blocks |
| 2 |  1 | MSFE | 12 |  5.519 | unavailable_sparse_or_insufficient_blocks |
| 2 |  1 | QLIKE | 12 |  4.077 | unavailable_sparse_or_insufficient_blocks |
| 2 | 10 | MSFE | 12 | 37.140 | unavailable_sparse_or_insufficient_blocks |
| 2 | 10 | QLIKE | 12 | 18.163 | unavailable_sparse_or_insufficient_blocks |
| 2 | 20 | MSFE | 12 | 14.402 | unavailable_sparse_or_insufficient_blocks |
| 2 | 20 | QLIKE | 12 | 13.733 | unavailable_sparse_or_insufficient_blocks |
| 3 |  1 | MSFE | 12 |  6.134 | unavailable_sparse_or_insufficient_blocks |
| 3 |  1 | QLIKE | 12 |  4.448 | unavailable_sparse_or_insufficient_blocks |
| 3 | 10 | MSFE | 12 | 38.010 | unavailable_sparse_or_insufficient_blocks |
| 3 | 10 | QLIKE | 12 | 18.200 | unavailable_sparse_or_insufficient_blocks |
| 3 | 20 | MSFE | 12 | 13.396 | unavailable_sparse_or_insufficient_blocks |
| 3 | 20 | QLIKE | 12 | 12.848 | unavailable_sparse_or_insufficient_blocks |
| 4 |  1 | MSFE | 12 |  5.722 | unavailable_sparse_or_insufficient_blocks |
| 4 |  1 | QLIKE | 12 |  4.118 | unavailable_sparse_or_insufficient_blocks |
| 4 | 10 | MSFE | 12 | 38.849 | unavailable_sparse_or_insufficient_blocks |
| 4 | 10 | QLIKE | 12 | 18.932 | unavailable_sparse_or_insufficient_blocks |
| 4 | 20 | MSFE | 12 | 13.209 | unavailable_sparse_or_insufficient_blocks |
| 4 | 20 | QLIKE | 12 | 12.867 | unavailable_sparse_or_insufficient_blocks |
| 5 |  1 | MSFE | 12 |  7.403 | unavailable_sparse_or_insufficient_blocks |
| 5 |  1 | QLIKE | 12 |  5.236 | unavailable_sparse_or_insufficient_blocks |
| 5 | 10 | MSFE | 12 | 41.587 | unavailable_sparse_or_insufficient_blocks |
| 5 | 10 | QLIKE | 12 | 20.416 | unavailable_sparse_or_insufficient_blocks |
| 5 | 20 | MSFE | 12 | 13.817 | unavailable_sparse_or_insufficient_blocks |
| 5 | 20 | QLIKE | 12 | 13.731 | unavailable_sparse_or_insufficient_blocks |

Положительный gain_pct означает меньшую ошибку mfOU. Отрицательный результат сохранён. n — число общих дат для данного горизонта, а не число независимых наблюдений.
Smoke проверяет запуск. Квартальная выборка — описательное исследование. Ежедневная выборка допускает парный блочный bootstrap при достаточном объёме и покрытии календаря. Пропущенные пары сохраняют своё положение в календаре; их нельзя сжать в соседние дни. Интервалы точечные и исследовательские, без поправки за множество сравнений.

## Что читать

metrics.csv — MSFE/RMSFE/QLIKE по горизонтам и периодам; comparisons.csv — выигрыш и статус неопределённости; window_sensitivity.csv — сравнение окон на общих датах; parameters.csv — параметры и слабая идентификация; estimation_failures.csv — сбои. Отсутствующий интервал не означает отсутствие статистической значимости: её здесь не оценивали.
Внутримодельная условная дисперсия не включает всю ошибку оценивания параметров. Данный отчёт не заменяет независимое ревью и не подтверждает улучшение статьи.
