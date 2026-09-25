# MENR: Mechano-Epigenetic Nano-Rewriter (computational framework)

Purely computational modeling of a theoretical mechano-epigenetic nanoparticle
therapy, applied to 16,911 TCGA patient-gene pairs across 4 cancer types.

## Structure
- `model-a-trigger/` — Model A, stiffness-gated activation (sigmoid)
- `model-b-kinetics/` — Model B, two-wave DNMT3A/TET1 methylation kinetics
- `digital-twin/` — full pipeline applied to real patient data
- `validation/` — sensitivity analysis, statistical validation, dashboard
- `final-outputs/` — summary figures and results

## Data
Patient-level TCGA data (mutation, copy-number, methylation) was sourced via
cBioPortal (Cerami et al. 2012; Gao et al. 2013) and is not redistributed
here due to size. Data acquisition scripts expect this data at
`data-acquisition/`.

## Dashboard
Interactive Streamlit dashboard: see `validation/menr_dashboard/`.

This project is purely computational — no wet-lab work was conducted.
