# Cell-level comparison of saved Monte Carlo tables with the paper (no simulation).
# Usage: Rscript scripts/04b_compare_mc_with_paper.R [smoke|full]
# Reads <results>/tables/code_a_mc_*.csv written by scripts/04_run_monte_carlo.R and writes
#   code_a_mc_agreement.csv          one row per paper cell: ours, paper, MC standard errors, z
#   code_a_mc_agreement_summary.csv  counts per table
# Criterion and caveats: see mc_agreement() in R/monte_carlo.R.

source("R/code_a_config.R")
source_code_a()
mode <- code_a_mode()
code_a_set_mode(mode)
code_a_log("04b_compare_mc_with_paper: mode = ", mode)
a <- mc_agreement(code_a_results("tables"))
a$source_mode <- mode
s <- mc_agreement_summary(a)
s$source_mode <- mode
code_a_write(a, "code_a_mc_agreement.csv")
code_a_write(s, "code_a_mc_agreement_summary.csv")
code_a_log(sprintf("  MC cells: %d of %d with |z| <= 2 (about 95%% expected if both sides estimate the same quantity)",
                   sum(a$agree & a$metric != "asym_se"), sum(a$metric != "asym_se")))
code_a_log(sprintf("  asymptotic SE cells equal after rounding: %d of %d",
                   sum(a$agree & a$metric == "asym_se"), sum(a$metric == "asym_se")))
code_a_log_run("04b_compare_mc_with_paper", mode)
