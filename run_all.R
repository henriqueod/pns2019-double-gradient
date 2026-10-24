# run_all.R
# Runs the full analysis pipeline in order, from the repository root:
#   Rscript run_all.R
# or, inside R/RStudio with the working directory set to the repository root:
#   source("run_all.R")
#
# Outputs go to results/ (CSV tables) under `root`, which defaults to the working
# directory; set the environment variable PNS_ROOT to use another folder.
# The PNS 2019 microdata are downloaded by PNSIBGE::get_pns() into data_pnsibge/
# on the first run (about 1 GB of disk space; internet access required).


required <- c("PNSIBGE", "survey", "dplyr", "ggplot2", "ranger")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the missing packages first: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))", call. = FALSE)
}

root <- Sys.getenv("PNS_ROOT", unset = getwd())
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

pipeline <- c(
  "R/01_prepare_cohort.R",
  "R/02_main_analysis.R",
  "R/03_income_bands.R",
  "R/04_depression_cascade.R",
  "R/05_detection_gap.R",
  "R/06_depression_model_checks.R",
  "R/07_random_forest.R",
  "R/08_referral_failure_blockD.R"
)

timing <- data.frame(script = pipeline, minutes = NA_real_)
for (i in seq_along(pipeline)) {
  f <- pipeline[i]
  message("\n==== ", f, " ====")
  t0 <- Sys.time()
  source(f, local = new.env(parent = globalenv()), encoding = "UTF-8")
  timing$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
}

# Record the computational environment next to the results.
dir.create(file.path(root, "session_info"), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "session_info", "sessionInfo.txt"))
write.csv(timing, file.path(root, "session_info", "pipeline_timing.csv"), row.names = FALSE)

message("\nPipeline finished. Tables are in results/. Figures: run\n",
        "  python python/make_figures.py\n  python python/graphical_abstract.py")
print(timing)
