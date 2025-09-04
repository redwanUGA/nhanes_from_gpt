# NHANES: % current smokers by 2-year cycle (adults >=20y), with survey design
# Python + rpy2 (to use R's 'survey' package correctly)
# Outputs: nhanes_current_smoking_by_cycle.csv and optional plot nhanes_current_smoking_by_cycle.png

import io
import os
import sys
import requests
import pandas as pd
import numpy as np
import pyreadstat
import matplotlib.pyplot as plt

# ---- rpy2 / R setup ----
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
pandas2ri.activate()
# Import R packages
survey = importr('survey')

# ---- cycles catalog ----
cycles = pd.DataFrame([
    ("1999-2000", "1999-2000", "_A", "https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/"),
    ("2001-2002", "2001-2002", "_B", "https://wwwn.cdc.gov/Nchs/Nhanes/2001-2002/"),
    ("2003-2004", "2003-2004", "_C", "https://wwwn.cdc.gov/Nchs/Nhanes/2003-2004/"),
    ("2005-2006", "2005-2006", "_D", "https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/"),
    ("2007-2008", "2007-2008", "_E", "https://wwwn.cdc.gov/Nchs/Nhanes/2007-2008/"),
    ("2009-2010", "2009-2010", "_F", "https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/"),
    ("2011-2012", "2011-2012", "_G", "https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/"),
    ("2013-2014", "2013-2014", "_H", "https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/"),
    ("2015-2016", "2015-2016", "_I", "https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/"),
    ("2017-2018", "2017-2018", "_J", "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/"),
    ("2021-2022", "2021-2022", "_M", "https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/"),
], columns=["cycle","years","suffix","base_url"])

# Optional: include the NCHS combined prepandemic file (2017–Mar 2020)
PREPANDEMIC_COMBINED = True
prepandemic_entry = {
    "cycle": "2017–Mar 2020 (prepandemic, combined)",
    "demo_url": "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_P.XPT",
    "smq_url":  "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQ_P.XPT",
    "combined": True
}

# ---- helpers ----
def read_xpt_from_url(url: str) -> pd.DataFrame:
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    df, meta = pyreadstat.read_xport(io.BytesIO(r.content))
    return df

def load_cycle(cyc_row: pd.Series) -> pd.DataFrame | None:
    suff = cyc_row["suffix"]
    base = cyc_row["base_url"]
    demo_url = f"{base}DEMO{suff}.XPT"
    smq_url  = f"{base}SMQ{suff}.XPT"
    try:
        demo = read_xpt_from_url(demo_url)
        smq  = read_xpt_from_url(smq_url)
    except Exception as e:
        print(f"[warn] Could not load {cyc_row['cycle']}: {e}")
        return None

    # Keep needed variables and merge
    demo_keep = demo[["SEQN","RIDAGEYR","SDMVSTRA","SDMVPSU","WTINT2YR"]].copy()
    smq_keep  = smq[["SEQN","SMQ020","SMQ040"]].copy()
    df = demo_keep.merge(smq_keep, on="SEQN", how="left")

    # Define adult and current_smoker
    df["adult"] = df["RIDAGEYR"] >= 20
    # current smoker: SMQ020==1 (≥100 cigs lifetime) AND SMQ040 in {1=every day, 2=some days}
    df["current_smoker"] = np.where(
        (df["SMQ020"]==1) & (df["SMQ040"].isin([1,2])),
        1,
        np.where(
            # non-smokers: never or former (SMQ020==2 OR SMQ020==1 & SMQ040==3)
            (df["SMQ020"]==2) | ((df["SMQ020"]==1) & (df["SMQ040"]==3)),
            0,
            np.nan
        )
    )
    # Adults only
    df = df.loc[df["adult"]==True].copy()
    return df

def survey_weighted_prev(df: pd.DataFrame) -> tuple[float,float,float,int]:
    """
    Use R 'survey' via rpy2 to compute mean(current_smoker==1) with
    WTINT2YR weights and SDMVSTRA/SDMVPSU design, Taylor linearization.
    Returns: (prev_pct, lcl_pct, ucl_pct, n_nonmissing)
    """
    # Drop rows with missing outcome
    sub = df[~df["current_smoker"].isna()].copy()
    n = len(sub)
    if n == 0:
        return (np.nan, np.nan, np.nan, 0)

    # Send to R
    rdf = pandas2ri.py2rpy(sub)

    # Design: ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR
    ro.globalenv["dat"] = rdf
    ro.r('des <- survey::svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=dat)')
    # Indicator in R
    ro.r('dat$current1 <- as.numeric(dat$current_smoker==1)')
    ro.r('des <- survey::svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTINT2YR, nest=TRUE, data=dat)')
    # Estimate mean and CI
    est = ro.r('survey::svymean(~current1, design=des, na.rm=TRUE)')
    ci  = ro.r('confint(svymean(~current1, design=des, na.rm=TRUE))')

    prev = float(est[0]) * 100.0
    lcl  = float(ci[0]) * 100.0
    ucl  = float(ci[1]) * 100.0
    return (prev, lcl, ucl, n)

# ---- run all standard cycles ----
rows = []
for _, cyc in cycles.iterrows():
    print(f"[info] Processing {cyc['cycle']} ...")
    df = load_cycle(cyc)
    if df is None:
        rows.append({"cycle": cyc["cycle"], "n_adults": np.nan, "prev": np.nan, "lcl": np.nan, "ucl": np.nan})
        continue
    prev, lcl, ucl, n = survey_weighted_prev(df)
    rows.append({"cycle": cyc["cycle"], "n_adults": n, "prev": prev, "lcl": lcl, "ucl": ucl})

# ---- optional: add prepandemic combined estimate ----
if PREPANDEMIC_COMBINED:
    try:
        print("[info] Processing 2017–Mar 2020 (prepandemic, combined) ...")
        demoP = read_xpt_from_url(prepandemic_entry["demo_url"])
        smqP  = read_xpt_from_url(prepandemic_entry["smq_url"])
        demo_keep = demoP[["SEQN","RIDAGEYR","SDMVSTRA","SDMVPSU","WTINT2YR"]].copy()  # combined interview weight
        smq_keep  = smqP[["SEQN","SMQ020","SMQ040"]].copy()
        dfP = demo_keep.merge(smq_keep, on="SEQN", how="left")
        dfP["adult"] = dfP["RIDAGEYR"] >= 20
        dfP["current_smoker"] = np.where(
            (dfP["SMQ020"]==1) & (dfP["SMQ040"].isin([1,2])), 1,
            np.where((dfP["SMQ020"]==2) | ((dfP["SMQ020"]==1) & (dfP["SMQ040"]==3)), 0, np.nan)
        )
        dfP = dfP.loc[dfP["adult"]==True].copy()
        prev, lcl, ucl, n = survey_weighted_prev(dfP)
        rows.append({
            "cycle": "2017–Mar 2020 (prepandemic, combined)",
            "n_adults": n, "prev": prev, "lcl": lcl, "ucl": ucl
        })
    except Exception as e:
        print(f"[warn] Prepandemic combined file failed: {e}")

results = pd.DataFrame(rows)

# ---- save results to CSV ----
out_csv = "nhanes_current_smoking_by_cycle.csv"
results.sort_values("cycle").to_csv(out_csv, index=False)
print(f"[info] Saved results to {out_csv}")

# ---- optional: make a simple plot with error bars ----
try:
    ordered = results.reset_index(drop=True)
    x = np.arange(len(ordered))
    plt.figure(figsize=(12,6))
    plt.errorbar(x, ordered["prev"], yerr=[ordered["prev"]-ordered["lcl"], ordered["ucl"]-ordered["prev"]], fmt='o-')
    plt.xticks(x, ordered["cycle"], rotation=60, ha='right')
    plt.xlabel("NHANES cycle")
    plt.ylabel("Current smoking prevalence (%)")
    plt.title("NHANES adults (≥20y): Current cigarette smoking by 2-year cycle")
    plt.tight_layout()
    plt.savefig("nhanes_current_smoking_by_cycle.png", dpi=150)
    print("[info] Saved plot to nhanes_current_smoking_by_cycle.png")
except Exception as e:
    print(f"[warn] Plotting failed: {e}")
