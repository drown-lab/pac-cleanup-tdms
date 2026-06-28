# Run the Python analysis steps. The external steps (MSConvert, FLASHDeconv,
# MSTopDiff, ProSightPD) must be done first and their outputs placed under data/
# (see README "Reproducing the analysis" and data/README.md).
#
# Override the interpreter if conda is not on PATH or you use a different env:
#   make all PY="conda run -n tdms python"
#   make all PY=python            # if the tdms env is already active

PY ?= conda run -n tdms python

.PHONY: all deconv-summary filter mods identification clean

all: deconv-summary filter mods identification

# --- deconvolution pipeline (needs data/flashdeconv/ from FLASHDeconv) ---
deconv-summary:
	$(PY) src/deconvolution/flashdeconv_summary.py

filter:
	$(PY) src/deconvolution/filter_features.py

# --- modification pipeline (needs data/flashdeconv/mstopdiff/ from MSTopDiff) ---
mods:
	$(PY) src/modifications/mstopdiff_compare.py
	$(PY) src/modifications/mstopdiff_unannotated.py

# --- identification pipeline (needs the .tdReport in data/) ---
identification:
	$(PY) src/identification/shared_proteoforms.py
	$(PY) src/identification/proteoform_physiochemical_props.py
	$(PY) src/identification/physiochemical_stats.py

# Remove generated outputs (figures + tables). Inputs under data/ are untouched.
clean:
	rm -f results/figures/* results/tables/*
