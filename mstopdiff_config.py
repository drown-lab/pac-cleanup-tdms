"""
Shared dataset registry for the MSTopDiff analysis scripts.

Auto-discovers the confident-mass MSTopDiff exports in flashdeconv/mstopdiff/
(<stem>_conf_mstopdiff.csv), labels and colours each via experimental_design.csv,
and orders them MCW -> Cytiva -> MagReSyn, by resuspension strength.

Run MSTopDiff on flashdeconv/filtered/<stem>_conf.tsv and drop the resulting
<stem>_conf_mstopdiff.csv into flashdeconv/mstopdiff/ to add a dataset; every
script picks it up automatically.
"""
import os
import sys
import glob
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MSTOPDIFF_DIR = os.path.join(HERE, "flashdeconv", "mstopdiff")
DESIGN = os.path.join(HERE, "experimental_design.csv")

# Short family names and a fixed resuspension display order.
FAMILY_SHORT = {"MCW": "MCW", "Cytiva Carboxyl": "Cyt", "MagReSyn Hydroxyl": "Hyd"}
RESUSP_ORDER = ["5% ACN, 0.1% FA", "0.5% TFA", "2% TFA", "10% FA", "20% FA"]

# Colours: MCW green; Cytiva blues; MagReSyn reds. 0.5% TFA matches the
# publication palette (MCW green / Cyt blue / Hyd red).
COLOR = {
    ("MCW", "5% ACN, 0.1% FA"): "#5AAF46",
    ("Cytiva Carboxyl", "0.5% TFA"): "#6482CD",
    ("Cytiva Carboxyl", "2% TFA"): "#3F5DA0",
    ("Cytiva Carboxyl", "10% FA"): "#27408B",
    ("Cytiva Carboxyl", "20% FA"): "#9DB0DE",
    ("MagReSyn Hydroxyl", "0.5% TFA"): "#ED645A",
    ("MagReSyn Hydroxyl", "2% TFA"): "#C0392B",
    ("MagReSyn Hydroxyl", "10% FA"): "#922B21",
    ("MagReSyn Hydroxyl", "20% FA"): "#F1948A",
}

# The three datasets intended for publication.
PUBLICATION = [("MCW", "5% ACN, 0.1% FA"),
               ("Cytiva Carboxyl", "0.5% TFA"),
               ("MagReSyn Hydroxyl", "0.5% TFA")]


def _design():
    d = pd.read_csv(DESIGN)
    d["stem"] = d["Name"].str.replace(r"\.raw$", "", regex=True)
    return d.set_index("stem", verify_integrity=True)


def _included(v):
    """Interpret the experimental_design Include flag (TRUE/FALSE strings)."""
    return str(v).strip().upper() not in ("FALSE", "0", "NO", "N", "NAN", "")


def discover():
    """Return ordered list of dataset dicts for every included CSV present.

    CSVs whose stem is absent from experimental_design.csv, or whose Include
    flag is FALSE, are skipped with a warning to stderr (so a misnamed/renamed
    file is never silently dropped).
    """
    design = _design()
    fam_rank = {"MCW": 0, "Cytiva Carboxyl": 1, "MagReSyn Hydroxyl": 2}
    out, skipped = [], []
    for path in sorted(glob.glob(os.path.join(MSTOPDIFF_DIR,
                                              "*_conf_mstopdiff.csv"))):
        stem = os.path.basename(path)[:-len("_conf_mstopdiff.csv")]
        if stem not in design.index:
            skipped.append((stem, "no matching row in experimental_design.csv"))
            continue
        if not _included(design.loc[stem, "Include"]):
            skipped.append((stem, "Include=FALSE"))
            continue
        cleanup = design.loc[stem, "Cleanup"]
        resusp = design.loc[stem, "Resuspension"]
        short = FAMILY_SHORT.get(cleanup, cleanup)
        label = short if cleanup == "MCW" else f"{short} {resusp}"
        out.append({
            "stem": stem, "csv": path, "cleanup": cleanup, "resusp": resusp,
            "label": label, "color": COLOR.get((cleanup, resusp), "#777777"),
            "is_pub": (cleanup, resusp) in PUBLICATION,
        })
    for stem, why in skipped:
        print(f"[mstopdiff_config] skipped {stem}_conf_mstopdiff.csv: {why}",
              file=sys.stderr)
    rr = {r: i for i, r in enumerate(RESUSP_ORDER)}
    out.sort(key=lambda d: (fam_rank.get(d["cleanup"], 9),
                            rr.get(d["resusp"], 9)))
    return out
