c2_write_csv <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(x,path,row.names=FALSE)
}

c2_write_empirical <- function(x,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  saveRDS(x,file.path(out,"empirical_run.rds"))
  tables <- list(forecasts=x$forecasts,parameters=x$parameters,origin_audit=x$audit,estimation_failures=x$failures,
    metrics=c2_metrics(x$forecasts),comparisons=c2_compare_scopes(x$forecasts,x$cfg),
    window_sensitivity=c2_window_comparison(x$forecasts),cumulative_loss=c2_loss_paths(x$forecasts))
  if(!ncol(tables$estimation_failures)) tables$estimation_failures<-data.frame(dimension=integer(),window=integer(),
    origin_date=as.Date(character()),origin_index=integer(),model=character(),reason=character(),elapsed_seconds=numeric())
  for(name in names(tables)) c2_write_csv(tables[[name]],file.path(out,paste0(name,".csv")))
  if(nrow(x$audit)) {
    z<-x$audit;z$year<-format(z$origin_date,"%Y")
    # Horizon audit rows and distinct origins are deliberately separate counts.
    counts<-aggregate(list(n_audit_rows=rep(1,nrow(z))),z[,c("dimension","window","year","reason")],sum)
    origins<-unique(z[,c("dimension","window","year","reason","origin_date")])
    origins<-aggregate(list(n_distinct_origins=rep(1,nrow(origins))),origins[,c("dimension","window","year","reason")],sum)
    c2_write_csv(merge(counts,origins),file.path(out,"origin_coverage.csv"))
  }
  c2_plot_empirical(x,tables,out)
  c2_make_report(x,tables,out)
  identity<-c2_source_identity()
  c2_write_csv(data.frame(path=names(identity),md5=unname(identity)),file.path(out,"analysis_source_identity.csv"))
  invisible(tables)
}

c2_plot_empirical <- function(x,tables,out) {
  figures<-file.path(out,"figures");dir.create(figures,recursive=TRUE,showWarnings=FALSE)
  # Reusing an output folder must not show a previous run's unavailable figure.
  unlink(file.path(figures,c("cumulative_loss.png","rolling_parameters.png")))
  paths<-tables$cumulative_loss
  if(nrow(paths)) {
    png(file.path(figures,"cumulative_loss.png"),width=1500,height=750,res=140)
    par(mfrow=c(1,2),mar=c(4,4,3,1))
    for(metric in c("cumulative_se_difference","cumulative_qlike_difference")) {
      z<-paths[paths$dimension==max(paths$dimension) & paths$window==x$cfg$window & paths$h %in% c(1L,20L),]
      if(!nrow(z)) {plot.new();next}
      limits<-range(c(0,z[[metric]]));dates<-range(z$origin_date)
      plot(dates,limits,type="n",xlab="Forecast origin",ylab="Cumulative mfBm loss minus mfOU loss",
        main=if(metric=="cumulative_se_difference") "Squared error (positive favors mfOU)" else "QLIKE (positive favors mfOU)")
      abline(h=0,col="grey60")
      for(h in unique(z$h)) {y<-z[z$h==h,];lines(y$origin_date,y[[metric]],col=if(h==1) "#245a91" else "#bd541d",lwd=2)}
      legend("topleft",paste0("h = ",unique(z$h)),col=ifelse(unique(z$h)==1,"#245a91","#bd541d"),lty=1,bty="o",bg="white")
    }
    dev.off()
  }
  params<-x$parameters
  z<-if(nrow(params)) params[params$dimension==max(params$dimension) & params$window==x$cfg$window,,drop=FALSE] else data.frame()
  if(nrow(z)) {
    png(file.path(figures,"rolling_parameters.png"),width=1500,height=750,res=140)
    par(mfrow=c(1,2),mar=c(4,4,3,1));colors<-hcl.colors(length(unique(z$asset)),"Dark 3")
    for(metric in c("H","kappa")) {
      plot(range(z$origin_date),range(z[[metric]]),type="n",xlab="Forecast origin",ylab=if(metric=="H") "H" else "kappa per year",
        main=if(metric=="H") "Rolling roughness estimates" else "Mean-reversion estimates (bounds disclosed)")
      for(i in seq_along(unique(z$asset))) {
        y<-z[z$asset==unique(z$asset)[i],];lines(y$origin_date,y[[metric]],col=colors[i],lwd=2)
        points(y$origin_date,y[[metric]],col=colors[i],pch=16,cex=.5)
        points(y$origin_date[y$boundary],y[[metric]][y$boundary],col=colors[i],pch=4)
      }
      legend("topright",unique(z$asset),col=colors,lty=1,bty="o",bg="white",cex=.8)
    }
    dev.off()
  }
}

c2_make_report <- function(x,tables,out) {
  metrics<-tables$metrics;comparison<-tables$comparisons
  csv_table <- function(df) {
    if(!nrow(df)) return("No estimable rows.")
    labels<-names(df);rows<-apply(df,1,function(r) paste0("| ",paste(r,collapse=" | ")," |"))
    paste(c(paste0("| ",paste(labels,collapse=" | ")," |"),paste0("| ",paste(rep("---",length(labels)),collapse=" | ")," |"),rows),collapse="\n")
  }
  paired<-if(nrow(comparison)) comparison[comparison$window==x$cfg$window & comparison$period=="full" & comparison$h %in% c(1L,10L,20L),
    c("dimension","h","metric","n","gain_pct","uncertainty_status"),drop=FALSE]
    else data.frame()
  if(nrow(paired)) paired$gain_pct<-round(paired$gain_pct,3)
  limits<-if(nrow(x$forecasts)) paste(min(x$forecasts$target_date),max(x$forecasts$target_date),sep=" to ") else "none"
  study_scope<-if(x$cfg$mode=="study") paste0("This study uses prespecified origins every ",x$cfg$study_step,
    " trading day(s), within ",x$cfg$study_start," to ",x$cfg$study_end,"; it is not the full daily paper period.") else
    if(x$cfg$mode=="smoke") "The full empirical paper period is not established by a smoke run." else
    "Full mode requests the configured period; actual dates, coverage and counts describe the executed scope."
  writeLines(c("# Code C v2: executed empirical report", "",
    paste("Mode:",x$cfg$mode,". Actual target dates:",limits,". Calendar: trading."),
    paste("Scheduled origins:",length(x$origins),"; generated forecast rows:",nrow(x$forecasts),"; mfOU failures:",nrow(x$failures),"."),
    "","## Scope", "",
    paste("Smoke proves execution only.",study_scope,"Inference is pointwise exploratory, with paired calendar blocks and no multiple-testing correction."),
    "Sparse or short samples receive no confidence intervals. Numerical gains below are descriptive; negative gains are retained. Baseline and mfOU losses use their exact shared keys; raw baseline forecasts and all failures remain available.",
    "","## Main 500-day comparison", "",csv_table(paired),"",
    "## Estimation diagnostics", "",
    paste("Boundary parameter rows:",if(nrow(x$parameters)) sum(x$parameters$boundary) else 0,"of",nrow(x$parameters),"."),
    "A boundary or a poorly conditioned local moment Jacobian signals weak identification; it does not prove mean reversion. Conditional log variance is plug-in and excludes parameter uncertainty.",
    "","## Files", "",
    "metrics.csv includes MSFE, RMSFE and QLIKE by model, horizon and prespecified paper subperiod. window_sensitivity.csv intersects keys before comparing training lengths. cumulative_loss.csv and figures retain the direction of every outcome. parameters.csv and estimation_failures.csv expose numerical/estimation limits.",
    "",if(file.exists(file.path(out,"figures/cumulative_loss.png"))) "![Cumulative loss](figures/cumulative_loss.png)" else "No paired loss plot is available.",
    "",if(file.exists(file.path(out,"figures/rolling_parameters.png"))) "![Rolling parameters](figures/rolling_parameters.png)" else "No rolling parameter plot is available.",
    "","Scientific definition and references: code_c/docs/METHODS_EN.md and METHODS_RU.md. No independent acceptance is asserted."),
    file.path(out,"REPORT_EN.md"))
  writeLines(c("# Часть C: выполненное сравнение", "",
    paste("Режим:",x$cfg$mode,". Фактические даты целевых значений:",limits,". Горизонт считается в торговых днях."),
    paste("Запланировано дат прогноза:",length(x$origins),"; строк прогнозов:",nrow(x$forecasts),"; ошибок mfOU:",nrow(x$failures),"."),
    "","## Основное окно 500 дней", "",csv_table(paired),"",
    "Положительный gain_pct означает меньшую ошибку mfOU. Отрицательный результат сохранён. n — число общих дат для данного горизонта, а не число независимых наблюдений.",
    "Smoke проверяет запуск. Квартальная выборка — описательное исследование. Ежедневная выборка допускает парный блочный bootstrap при достаточном объёме и покрытии календаря. Пропущенные пары сохраняют своё положение в календаре; их нельзя сжать в соседние дни. Интервалы точечные и исследовательские, без поправки за множество сравнений.",
    "","## Что читать", "",
    "metrics.csv — MSFE/RMSFE/QLIKE по горизонтам и периодам; comparisons.csv — выигрыш и статус неопределённости; window_sensitivity.csv — сравнение окон на общих датах; parameters.csv — параметры и слабая идентификация; estimation_failures.csv — сбои. Отсутствующий интервал не означает отсутствие статистической значимости: её здесь не оценивали.",
    "Внутримодельная условная дисперсия не включает всю ошибку оценивания параметров. Данный отчёт не заменяет независимое ревью и не подтверждает улучшение статьи."),file.path(out,"REPORT_RU.md"))
}

c2_workload_plan <- function(panel,cfg) {
  origins<-c2_origins(panel,cfg)
  c2_bind(lapply(c2_scenarios(cfg),function(sc) {
    el<-c2_joint_eligibility(panel,cfg,sc$window,origins)
    size<-sc$window*length(sc$assets)
    data.frame(dimension=length(sc$assets),window=sc$window,n_scheduled=length(origins),
      n_training_eligible=sum(el$eligible),covariance_matrix_MB=8*size^2/1024^2,
      conservative_working_memory_MB=6*8*size^2/1024^2,
      scope="training_preflight_not_completed_forecasts")
  }))
}
