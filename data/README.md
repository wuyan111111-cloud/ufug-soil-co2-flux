# Local data directory

The machine-readable `collar_flux_198_public.csv` is included here to make the
analysis immediately reproducible. The archived Chinese-header Excel workbook
will be preserved in Mendeley Data as the provenance copy.

The analysis accepts either of these files:

- `original_data.xlsx` — archived original workbook with Chinese headers; or
- `collar_flux_198_public.csv` — machine-readable public CSV with English
  headers.

Run:

```bash
python src/reanalysis.py data/collar_flux_198_public.csv --outdir results
```

The raw workbook must remain unchanged. Manuscript display-label
standardization is documented separately in the Mendeley Data metadata.
Generated public tables use the English `Plant_background` field and the Latin
names listed in `metadata/csu_plant_backgrounds.csv`.
