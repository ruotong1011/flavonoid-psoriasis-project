# RSC / Food & Function code and data checklist

- [x] Include a Data Availability Statement in the manuscript.
- [x] Provide public access to analysis code through GitHub.
- [x] Explain why UK Biobank raw data cannot be publicly shared.
- [x] Document all required input files in `data/README_data.md`.
- [x] Remove local absolute paths such as `C:/Users/...`.
- [x] Do not upload participant-level data, participant IDs, or restricted UK Biobank files.
- [x] Add a README explaining how to reproduce the analysis.
- [x] Use `renv` or equivalent to document R package dependencies.
- [x] Use `targets` to make the workflow easier to rerun.
- [ ] After acceptance or final submission, create a GitHub release and archive it on Zenodo to obtain a DOI.
