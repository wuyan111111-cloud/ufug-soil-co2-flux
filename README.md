# Short-term soil CO₂ flux responses in an engineered urban green corridor

This repository contains the maintained analysis and figure-generation code
for the manuscript:

> *Short-term soil CO₂ flux responses to soil amendments in highly
> heterogeneous engineered urban green spaces: A stratified paired-effect
> analysis along a linear corridor*

## Study design represented by the code

- 11 existing configuration-site units (CSUs)
- 3 randomized spatial blocks per CSU
- 3 soil treatments per block: CK, adaptive DCD, and biochar
- 99 permanent soil collars
- 2 short monitoring windows: July and October 2023
- 198 Collar-level flux observations
- 132 within-stratum treatment/CK log-ratio effects

The code does not rank plant configurations. CSU-specific estimates are
descriptive because plant composition and site position were not spatially
replicated as independent treatments.

## Configuration-site unit backgrounds

Latin names follow accepted GBIF Backbone Taxonomy matches. Botanical names
are italicized here; the archived field workbook retains the original Chinese
labels for provenance.

| CSU | Existing plant background |
|---|---|
| Lawn | Managed lawn; turf species were not recorded |
| PB | *Pinus bungeana* |
| PG | *Phyllostachys glauca* |
| SJ | *Styphnolobium japonicum* |
| NZE | *Ligustrum lucidum* × *Iris tectorum* |
| GB×CH | *Ginkgo biloba* × *Cotoneaster horizontalis* |
| SJ×JN | *Styphnolobium japonicum* × *Jasminum nudiflorum* |
| OF×RC | *Osmanthus fragrans* × *Rosa chinensis* |
| AR×ND | *Acer palmatum* × *Nandina domestica* |
| GYL×EF | *Magnolia grandiflora* × *Euonymus fortunei* |
| GYL×ND | *Magnolia grandiflora* × *Nandina domestica* |

## Data

The frozen data and code snapshot used by the manuscript will be archived in
Mendeley Data. Its DOI and version will be added here when the dataset record
is published.

The repository includes the public, machine-readable
`data/collar_flux_198_public.csv`. The archived Chinese-header Excel workbook
will be preserved in Mendeley Data as the provenance copy.

The analysis accepts either file format.

## Reproduce the statistical outputs

```bash
conda env create -f environment.yml
conda activate ufug-co2-flux
python src/reanalysis.py data/collar_flux_198_public.csv --outdir results
```

The command writes:

- `collar_flux_198_public.csv`
- `paired_effects_132.csv`
- `descriptive_summary.csv`
- `primary_csu_clustered_effects.csv`
- `sensitivity_block_clustered_effects.csv`
- `window_interaction_tests.csv`
- `csu_window_effects.csv`
- `exploratory_csu_effects.csv`
- `data_integrity_checks.csv`
- `run_summary.json`

## Regenerate figures

```bash
python src/make_figures.py --results results --outdir figures
```

The script exports PNG, PDF, SVG, and 600-dpi TIFF files for Figure 1,
Figure 2, and Supplementary Figure S1.

## Statistical specification

The response is the log ratio
`ln(F_treatment / F_CK)` within each CSU × block × monitoring-window stratum.
The no-intercept mean model estimates four pre-specified treatment × window
means. The primary covariance is cluster-robust by CSU (11 clusters, t
inference with df = 10); the sensitivity covariance is clustered by randomized
block (33 clusters, df = 32). Small-sample covariance correction is enabled.

## Evidence boundary

The outputs describe post-treatment soil-surface CO₂ flux during two short
field windows. They do not estimate annual emissions, treatment mechanisms,
plant-configuration effects, or net ecosystem carbon sequestration.

## Software

The manuscript analysis was validated with Python 3.12.10 and the exact package
versions in `requirements.txt` and `environment.yml`.

## Citation and licence

The formal software citation and licence will be added after the author list,
institutional licence choice, and Mendeley Data DOI are confirmed. Until a
licence is added, standard copyright restrictions apply.
