# ensure required packages are installed and loaded
required_pkgs <- c("tidyverse", "haven", "survey", "srvyr")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

# -------- helper: cycle catalog --------
cycles <- tribble(
  ~cycle,        ~years,       ~suffix, ~base_url,
  "1999-2000",   "1999-2000",  "_A",    "https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/",
  "2001-2002",   "2001-2002",  "_B",    "https://wwwn.cdc.gov/Nchs/Nhanes/2001-2002/",
  "2003-2004",   "2003-2004",  "_C",    "https://wwwn.cdc.gov/Nchs/Nhanes/2003-2004/",
  "2005-2006",   "2005-2006",  "_D",    "https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/",
  "2007-2008",   "2007-2008",  "_E",    "https://wwwn.cdc.gov/Nchs/Nhanes/2007-2008/",
  "2009-2010",   "2009-2010",  "_F",    "https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/",
  "2011-2012",   "2011-2012",  "_G",    "https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/",
  "2013-2014",   "2013-2014",  "_H",    "https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/",
  "2015-2016",   "2015-2016",  "_I",    "https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/",
  "2017-2018",   "2017-2018",  "_J",    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/",
  "2021-2022",   "2021-2022",  "_M",    "https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/"
)

# File stems by component (DEMO, SMQ changed naming slightly over time but suffix pattern works)
file_stem <- list(
  demo = "DEMO",   # e.g., DEMO_J.XPT
  smq  = "SMQ"     # e.g., SMQ_J.XPT  (pre-2017 sometimes SMQ is split; this stem covers the main file)
)

# -------- helper: read one cycle (from web URL or local path) --------
read_cycle <- function(cyc, src="web"){
  suff <- cycles %>% filter(cycle==cyc) %>% pull(suffix)
  base <- cycles %>% filter(cycle==cyc) %>% pull(base_url)
  
  demo_url <- paste0(base, file_stem$demo, suff, ".XPT")
  smq_url  <- paste0(base, file_stem$smq,  suff, ".XPT")
  
  demo <- tryCatch(read_xpt(if (src=="web") demo_url else file.path(src, paste0(file_stem$demo, suff, ".XPT"))),
                   error=function(e) NULL)
  smq  <- tryCatch(read_xpt(if (src=="web") smq_url  else file.path(src,  paste0(file_stem$smq,  suff, ".XPT"))),
                   error=function(e) NULL)
  
  if (is.null(demo) || is.null(smq)) return(NULL)
  
  dat <- demo %>%
    select(SEQN, RIDAGEYR, SDMVSTRA, SDMVPSU, WTINT2YR) %>%
    left_join(
      smq %>% select(SEQN, SMQ020, SMQ040),
      by="SEQN"
    ) %>%
    mutate(
      adult = RIDAGEYR >= 20,
      current_smoker = case_when(
        SMQ020 == 1 & SMQ040 %in% c(1,2) ~ 1,
        SMQ020 %in% c(2) ~ 0,
        SMQ020 == 1 & SMQ040 == 3 ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(adult)
  
  dat
}

# -------- analyze one cycle --------
estimate_cycle <- function(cyc, src="web"){
  dat <- read_cycle(cyc, src=src)
  if (is.null(dat)) return(tibble(
    cycle = cyc, n_adults = NA_integer_, prev = NA_real_, lcl = NA_real_, ucl = NA_real_
  ))
  
  des <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTINT2YR,
    nest = TRUE, data = dat
  )
  est <- svymean(~I(current_smoker==1), design = des, na.rm = TRUE)
  ci  <- confint(est)
  tibble(
    cycle   = cyc,
    n_adults = sum(!is.na(dat$current_smoker)),
    prev = as.numeric(est)[1]*100,
    lcl  = ci[1]*100,
    ucl  = ci[2]*100
  )
}

# -------- run all standard cycles --------
results <- map_dfr(cycles$cycle, estimate_cycle)

# -------- OPTIONAL: 2017–March 2020 prepandemic combined file --------
# NCHS created combined 2017–2018 + 2019–Mar2020 files with special weights.
# We read the combined DEMO and SMQ directly (file names end with "_P" on the NHANES site).
prepand_base <- "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/"
prepand_demo <- paste0(prepand_base, "DEMO_P.XPT")
prepand_smq  <- paste0(prepand_base, "SMQ_P.XPT")

try({
  demoP <- read_xpt(prepand_demo)
  smqP  <- read_xpt(prepand_smq)
  
  datP <- demoP %>%
    select(SEQN, RIDAGEYR, SDMVSTRA, SDMVPSU, WTINT2YR) %>%  # WTINT2YR in _P is the combined interview weight
    left_join(smqP %>% select(SEQN, SMQ020, SMQ040), by="SEQN") %>%
    mutate(
      adult = RIDAGEYR >= 20,
      current_smoker = case_when(
        SMQ020 == 1 & SMQ040 %in% c(1,2) ~ 1,
        SMQ020 %in% c(2) ~ 0,
        SMQ020 == 1 & SMQ040 == 3 ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(adult)
  
  desP <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=datP)
  estP <- svymean(~I(current_smoker==1), design=desP, na.rm=TRUE)
  ciP  <- confint(estP)
  
  results <- bind_rows(
    results,
    tibble(
      cycle = "2017–Mar 2020 (prepandemic, combined)",
      n_adults = sum(!is.na(datP$current_smoker)),
      prev = as.numeric(estP)[1]*100,
      lcl  = ciP[1]*100,
      ucl  = ciP[2]*100
    )
  )
}, silent = TRUE)

# -------- save & plot --------
readr::write_csv(results %>% arrange(cycle), "nhanes_current_smoking_by_cycle.csv")

png("nhanes_current_smoking_by_cycle.png", width=1100, height=600)
plot(
  x = seq_len(nrow(results)),
  y = results$prev,
  type = "b",
  xaxt = "n",
  xlab = "NHANES cycle",
  ylab = "Current smoking prevalence (%)",
  main = "NHANES adults (≥20y): Current cigarette smoking by 2-year cycle"
)
axis(1, at = seq_len(nrow(results)), labels = results$cycle, las = 2, cex.axis = 0.8)
arrows(x0=seq_len(nrow(results)), y0=results$lcl, x1=seq_len(nrow(results)), y1=results$ucl, angle=90, code=3, length=0.05)
dev.off()

results %>% arrange(cycle)
