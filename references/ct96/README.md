# ct96 shortened FR reference set

This folder contains a shortened full-resolution reference case for `ct96-24-13`.

Files:
- `12664_ctd_shortened.txt`: 3 days of HR raw text data, covering 20 reference profiles.
- `ct96-24-13_fr0_prof_shortened.nc`: reference FR0 profile subset.
- `ct96-24-13_fr1_prof_shortened.nc`: reference FR1 profile subset for later TLC/FR validation.

Notes:
- The shortened NetCDF references preserve metadata from the original production files.
- In particular, the embedded `instr_id` attribute still reads `12661` in the provided reference NetCDFs,
  while the corresponding raw HR text filename and deployment catalog entry use `12664`.
  The current Python tests therefore compare core profile content rather than that specific attribute.
