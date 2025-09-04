# NHANES: % current smokers by 2-year cycle (adults >=20y), with survey design
# Pure Python implementation (Windows-friendly) using weighted means and
# Kish effective sample size to approximate 95% CI (no rpy2/R dependency).
# Outputs: nhanes_current_smoking_by_cycle.csv and optional plot nhanes_current_smoking_by_cycle.png

import io
import os
import sys
import requests
import pandas as pd
import numpy as np
import pyreadstat
import matplotlib.pyplot as plt

# ---- Survey computation: pure Python (no rpy2) ----
# We avoid rpy2/R on Windows for compatibility. We will compute a
# weighted prevalence and an approximate CI using Kish's effective
# sample size. This is a reasonable approximation when full Taylor
# linearization (as in R's 'survey') is not available.
_APPROX_WARNING_SHOWN = False

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
import re
from urllib.parse import urljoin, urlparse, parse_qs, unquote
def read_xpt_from_url(url: str) -> pd.DataFrame:
    # Use a session to preserve cookies the server may set on the doc page
    session = requests.Session()
    # Determine referer: for direct XPT, replace .XPT with .htm; for DownloadXpt.aspx, build from query params
    if 'downloadxpt.aspx' in url.lower():
        try:
            pr = urlparse(url)
            qs = parse_qs(pr.query)
            file_base = (qs.get('FileName',[None])[0] or '').strip()
            pth = unquote((qs.get('Path',[None])[0] or ''))
            ftpname = (qs.get('ftpname',[None])[0] or '').strip()
            if pth and file_base:
                ref = f"https://{pr.netloc}{pth}{file_base}.htm"
            elif ftpname and file_base:
                ref = f"https://{pr.netloc}/Nchs/Nhanes/{ftpname}/{file_base}.htm"
            else:
                ref = url
        except Exception:
            ref = url
    else:
        ref = re.sub(r'\.XPT$', '.htm', url, flags=re.IGNORECASE)
    base_headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36',
        'Accept': '*/*',
        'Connection': 'keep-alive',
    }
    # Prime cookies by visiting referer (documentation page), if looks like a CDC doc URL
    try:
        if ('wwwn.cdc.gov' in ref.lower()) and ref.lower().endswith('.htm'):
            session.get(ref, headers={**base_headers, 'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'}, timeout=60, allow_redirects=True)
    except Exception:
        pass
    # Now download the XPT with Referer header and session cookies
    xpt_headers = {**base_headers, 'Referer': ref}
    r = session.get(url, headers=xpt_headers, timeout=60, allow_redirects=True)
    ctype = r.headers.get('Content-Type', '')
    data = r.content
    # If the server returns a smart 404 with 200 status, detect HTML body
    if 'html' in ctype.lower() or (data[:15].lstrip().lower().startswith(b'<!doctype html') or data[:6].lstrip().lower().startswith(b'<html>')):
        raise RuntimeError(f"Server returned HTML instead of XPT for {url} (Content-Type={ctype}).")
    r.raise_for_status()
    # Try pyreadstat first (read from in-memory buffer)
    try:
        df, meta = pyreadstat.read_xport(io.BytesIO(data))
        return df
    except Exception as e1:
        # Fallback to pandas
        try:
            with io.BytesIO(data) as bio:
                df = pd.read_sas(bio, format='xport')
            return df
        except Exception as e2:
            # Attach content-type info for debugging
            snippet = ''
            try:
                text = data[:300].decode('utf-8', errors='ignore')
                snippet = text.replace('\n', ' ')[:200]
            except Exception:
                pass
            raise RuntimeError(f"Failed to read XPT from {url}. Content-Type={ctype}. pyreadstat: {e1}. pandas: {e2}. First bytes as text: {snippet}")

def list_ftp_xpt_urls(years: str) -> list[str]:
    """List .XPT file URLs in the NHANES ftp folder for a given years label."""
    bases = [
        f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/nhanes/{years}/",
        f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/NHANES/{years}/",
    ]
    urls: list[str] = []
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36'
    }
    for base in bases:
        try:
            r = requests.get(base, headers=headers, timeout=60)
            if r.status_code != 200:
                continue
            text = r.text
            # Try HTML hrefs first
            found = False
            for href in re.findall(r'href=[\"\']([^\"\']+)', text, flags=re.IGNORECASE):
                if href.lower().endswith('.xpt'):
                    urls.append(urljoin(base, href))
                    found = True
            if not found:
                # Fallback: scan plain text listing for .XPT tokens
                for token in re.findall(r'[^\s]+\.xpt', text, flags=re.IGNORECASE):
                    urls.append(urljoin(base, token))
        except Exception:
            continue
    # Deduplicate while preserving order
    seen = set()
    out = []
    for u in urls:
        if u not in seen:
            out.append(u)
            seen.add(u)
    return out

def discover_xpt_url(years: str, token: str, suffix: str) -> str | None:
    """Find the best .XPT URL in ftp listing for given token (e.g., DEMO, SMQ) and suffix (e.g., _A)."""
    try:
        all_urls = list_ftp_xpt_urls(years)
    except Exception:
        all_urls = []
    if not all_urls:
        return None
    token_up = token.upper()
    suffix_up = suffix.upper()
    # Priorities: exact token+suffix, then token_*, then contains token
    # 1) Exact match
    for u in all_urls:
        name = u.split('/')[-1].upper()
        if name == f"{token_up}{suffix_up}.XPT":
            return u
    # 2) token with any suffix like token_?.XPT or token?.XPT
    for u in all_urls:
        name = u.split('/')[-1].upper()
        if name.startswith(token_up) and name.endswith('.XPT'):
            return u
    # 3) contains token
    for u in all_urls:
        if token_up in u.upper() and u.upper().endswith('.XPT'):
            return u
    return None

def discover_xpt_from_doc(years: str, token: str, suffix: str) -> str | None:
    """
    Try multiple component documentation URL variants and parse for an .XPT or DownloadXpt link.
    """
    token_up = token.upper()
    suffix_up = suffix.upper()
    candidate_docs = [
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{years}/{token_up}{suffix_up}.htm",
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{years}/{token_up}.htm",
        f"https://wwwn.cdc.gov/nchs/nhanes/{years}/{token_up}{suffix_up}.htm",
        f"https://wwwn.cdc.gov/nchs/nhanes/{years}/{token_up}.htm",
    ]
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Referer': f"https://wwwn.cdc.gov/Nchs/Nhanes/continuousnhanes/default.aspx",
    }
    for doc_url in candidate_docs:
        try:
            r = requests.get(doc_url, headers=headers, timeout=60)
            if r.status_code != 200:
                continue
            html = r.text
            # Look for direct .XPT links
            matches = re.findall(r'href=[\"\']([^\"\']+\.xpt)\b', html, flags=re.IGNORECASE)
            for href in matches:
                name = href.split('/')[-1]
                if re.match(fr'^{re.escape(token_up)}{re.escape(suffix_up)}\.xpt$', name, flags=re.IGNORECASE) or name.upper().startswith(token_up):
                    return urljoin(doc_url, href)
            # Look for DownloadXpt.aspx links
            matches2 = re.findall(r'href=[\"\']([^\"\']*DownloadXpt\.aspx\?[^\"\']*)', html, flags=re.IGNORECASE)
            for href in matches2:
                if re.search(fr'FileName=({re.escape(token_up)}{re.escape(suffix_up)}|{re.escape(token_up)})', href, flags=re.IGNORECASE) and (
                    re.search(fr'Path=[^&]*/Nchs/Nhanes/{re.escape(years)}/', href, flags=re.IGNORECASE) or
                    re.search(fr'ftpname={re.escape(years)}', href, flags=re.IGNORECASE)
                ):
                    return urljoin(doc_url, href)
            # Fallback: any .XPT or DownloadXpt link on the page
            if matches:
                return urljoin(doc_url, matches[0])
            if matches2:
                return urljoin(doc_url, matches2[0])
        except Exception:
            continue
    return None

def build_downloadxpt_url(years: str, token: str, suffix: str) -> str:
    """Heuristically construct the CDC DownloadXpt.aspx URL for a component.
    Example (Path style): https://wwwn.cdc.gov/Nchs/Nhanes/DownloadXpt.aspx?FileName=DEMO_I&FileType=XPT&Path=/Nchs/Nhanes/2015-2016/
    """
    token_up = token.upper()
    suffix_up = suffix.upper()
    file_base = f"{token_up}{suffix_up}".strip()
    path = f"/Nchs/Nhanes/{years}/"
    return (
        "https://wwwn.cdc.gov/Nchs/Nhanes/DownloadXpt.aspx?" +
        f"FileName={file_base}&FileType=XPT&Path={requests.utils.quote(path, safe='/') }"
    )

def build_downloadxpt_url_ftpname(years: str, token: str, suffix: str) -> str:
    """Alternate DownloadXpt.aspx URL using ftpname parameter instead of Path."""
    token_up = token.upper()
    suffix_up = suffix.upper()
    file_base = f"{token_up}{suffix_up}".strip()
    return (
        "https://wwwn.cdc.gov/Nchs/Nhanes/DownloadXpt.aspx?" +
        f"FileName={file_base}&FileType=XPT&ftpname={years}"
    )

def load_cycle(cyc_row: pd.Series) -> pd.DataFrame | None:
    suff = cyc_row["suffix"]
    base = cyc_row["base_url"]
    years = cyc_row["years"]
    # Candidate URLs in priority order
    primary_with_suffix_demo = f"{base}DEMO{suff}.XPT"
    primary_with_suffix_smq  = f"{base}SMQ{suff}.XPT"
    primary_no_suffix_demo   = f"{base}DEMO.XPT"
    primary_no_suffix_smq    = f"{base}SMQ.XPT"
    ftp_base = f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/nhanes/{years}/"
    ftp_base_caps = f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/NHANES/{years}/"
    ftp_with_suffix_demo = f"{ftp_base}DEMO{suff}.XPT"
    ftp_with_suffix_smq  = f"{ftp_base}SMQ{suff}.XPT"
    ftp_no_suffix_demo   = f"{ftp_base}DEMO.XPT"
    ftp_no_suffix_smq    = f"{ftp_base}SMQ.XPT"
    ftp_caps_with_suffix_demo = f"{ftp_base_caps}DEMO{suff}.XPT"
    ftp_caps_with_suffix_smq  = f"{ftp_base_caps}SMQ{suff}.XPT"
    ftp_caps_no_suffix_demo   = f"{ftp_base_caps}DEMO.XPT"
    ftp_caps_no_suffix_smq    = f"{ftp_base_caps}SMQ.XPT"

    tried = []
    demo = smq = None

    # Explicit overrides for cycles with known DataFiles links
    if years == "2005-2006":
        override_demo = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DEMO_D.xpt"
        override_smq  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/SMQ_D.xpt"
        try:
            demo = read_xpt_from_url(override_demo)
            smq  = read_xpt_from_url(override_smq)
        except Exception as e:
            tried.append((override_demo, override_smq, str(e)))
            demo = smq = None

    if demo is None or smq is None:
        for d_url, s_url in [
            (primary_with_suffix_demo, primary_with_suffix_smq),
            (primary_no_suffix_demo, primary_no_suffix_smq),
            (ftp_with_suffix_demo, ftp_with_suffix_smq),
            (ftp_no_suffix_demo, ftp_no_suffix_smq),
            (ftp_caps_with_suffix_demo, ftp_caps_with_suffix_smq),
            (ftp_caps_no_suffix_demo, ftp_caps_no_suffix_smq),
        ]:
            try:
                demo = read_xpt_from_url(d_url)
                smq  = read_xpt_from_url(s_url)
                break
            except Exception as e:
                tried.append((d_url, s_url, str(e)))
                demo = smq = None
                continue

    if demo is None or smq is None:
        # Try the official DownloadXpt.aspx endpoint (Path style first, then ftpname style)
        try:
            d_dl = build_downloadxpt_url(years, "DEMO", suff)
            s_dl = build_downloadxpt_url(years, "SMQ", suff)
            demo = read_xpt_from_url(d_dl)
            smq  = read_xpt_from_url(s_dl)
        except Exception as e_dl:
            tried.append((d_dl if 'd_dl' in locals() else 'n/a', s_dl if 's_dl' in locals() else 'n/a', str(e_dl)))
            demo = smq = None
        if demo is None or smq is None:
            try:
                d_dl2 = build_downloadxpt_url_ftpname(years, "DEMO", suff)
                s_dl2 = build_downloadxpt_url_ftpname(years, "SMQ", suff)
                demo = read_xpt_from_url(d_dl2)
                smq  = read_xpt_from_url(s_dl2)
            except Exception as e_dl2:
                tried.append((d_dl2 if 'd_dl2' in locals() else 'n/a', s_dl2 if 's_dl2' in locals() else 'n/a', str(e_dl2)))
                demo = smq = None
        # Try discovery from documentation pages next
        if demo is None or smq is None:
            d_doc = discover_xpt_from_doc(years, "DEMO", suff)
            s_doc = discover_xpt_from_doc(years, "SMQ", suff)
            if d_doc and s_doc:
                try:
                    demo = read_xpt_from_url(d_doc)
                    smq  = read_xpt_from_url(s_doc)
                except Exception as e3:
                    tried.append((d_doc, s_doc, str(e3)))
                    demo = smq = None
        # If still not found, try discovery from ftp directory listing
        if demo is None or smq is None:
            d_disc = discover_xpt_url(years, "DEMO", suff)
            s_disc = discover_xpt_url(years, "SMQ", suff)
            if d_disc and s_disc:
                try:
                    demo = read_xpt_from_url(d_disc)
                    smq  = read_xpt_from_url(s_disc)
                except Exception as e4:
                    tried.append((d_disc, s_disc, str(e4)))
                    demo = smq = None
        # Last resort: use the static DataFiles path by first year of cycle
        if demo is None or smq is None:
            try:
                first_year = str(years).split('-')[0]
                letter = str(suff).strip('_')
                df_demo = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{first_year}/DataFiles/DEMO_{letter}.xpt"
                df_smq  = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{first_year}/DataFiles/SMQ_{letter}.xpt"
                demo = read_xpt_from_url(df_demo)
                smq  = read_xpt_from_url(df_smq)
            except Exception as e5:
                tried.append((df_demo if 'df_demo' in locals() else 'n/a', df_smq if 'df_smq' in locals() else 'n/a', str(e5)))
                demo = smq = None
        if demo is None or smq is None:
            print(f"[warn] Could not load {cyc_row['cycle']} after trying multiple sources:")
            for (d_url, s_url, err) in tried:
                print(f"       - DEMO: {d_url} | SMQ: {s_url} | err: {err}")
            if 'd_dl' in locals() or 's_dl' in locals():
                print(f"       - DownloadXpt: DEMO: {d_dl} | SMQ: {s_dl}")
            if 'd_doc' in locals() or 's_doc' in locals():
                print(f"       - Doc page discovery: DEMO: {d_doc} | SMQ: {s_doc}")
            if 'd_disc' in locals() or 's_disc' in locals():
                print(f"       - FTP discovery: DEMO: {d_disc} | SMQ: {s_disc}")
            if 'df_demo' in locals() or 'df_smq' in locals():
                print(f"       - DataFiles: DEMO: {df_demo} | SMQ: {df_smq}")
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
    Compute weighted prevalence of current_smoker (1/0) using WTINT2YR.
    95% CI approximated via Kish effective sample size (neff).
    Returns: (prev_pct, lcl_pct, ucl_pct, n_nonmissing)
    Note: This approximates R's survey Taylor SE; may be slightly narrower.
    """
    global _APPROX_WARNING_SHOWN
    if not _APPROX_WARNING_SHOWN:
        print("[note] Using approximate CI via Kish effective sample size (no rpy2/R 'survey').")
        _APPROX_WARNING_SHOWN = True
    # Drop rows with missing outcome
    sub = df[~df["current_smoker"].isna()].copy()
    n = len(sub)
    if n == 0:
        return (np.nan, np.nan, np.nan, 0)

    y = sub["current_smoker"].astype(float).to_numpy()
    w = sub["WTINT2YR"].astype(float).to_numpy()
    # Guard against nonpositive or NaN weights
    mask = np.isfinite(w) & (w > 0) & np.isfinite(y)
    y = y[mask]
    w = w[mask]
    if w.size == 0 or w.sum() == 0:
        return (np.nan, np.nan, np.nan, int(mask.sum()))

    # Weighted mean
    p = float(np.sum(w * y) / np.sum(w))

    # Kish effective sample size
    sumw = np.sum(w)
    sumw2 = np.sum(w ** 2)
    neff = (sumw ** 2) / sumw2 if sumw2 > 0 else np.nan

    # Standard error and 95% CI (Wald)
    if np.isnan(neff) or neff <= 0:
        lcl = ucl = p
    else:
        se = np.sqrt(max(p * (1 - p), 0.0) / neff)
        z = 1.96
        lcl = max(0.0, p - z * se)
        ucl = min(1.0, p + z * se)

    return (p * 100.0, lcl * 100.0, ucl * 100.0, int(len(sub)))

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
    fig, ax = plt.subplots(figsize=(12,6))
    # Bar chart with error bars (95% CI) and value labels on top
    heights = ordered["prev"].to_numpy()
    # Compute asymmetric error bars only if lcl/ucl present
    if {"lcl", "ucl"}.issubset(ordered.columns):
        yerr = [ordered["prev"] - ordered["lcl"], ordered["ucl"] - ordered["prev"]]
    else:
        yerr = None
    bars = ax.bar(x, heights, yerr=yerr, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(ordered["cycle"], rotation=60, ha='right')
    ax.set_xlabel("NHANES cycle")
    ax.set_ylabel("Current smoking prevalence (%)")
    ax.set_title("NHANES adults (≥20y): Current cigarette smoking by 2-year cycle")
    # Ensure some headroom for labels
    try:
        max_ucl = np.nanmax(ordered["ucl"]) if "ucl" in ordered.columns else np.nanmax(heights)
        ax.set_ylim(0, max_ucl + 5)
    except Exception:
        pass
    # Add value labels on top of each bar
    for bar, val in zip(bars, heights):
        if not np.isnan(val):
            ax.annotate(f"{val:.1f}",
                        xy=(bar.get_x() + bar.get_width() / 2, val),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    fig.tight_layout()
    fig.savefig("nhanes_current_smoking_by_cycle.png", dpi=150)
    print("[info] Saved plot to nhanes_current_smoking_by_cycle.png")
except Exception as e:
    print(f"[warn] Plotting failed: {e}")
