# CFD/MIT Scoring Tables

The library ships with baked-in placeholder weights for CFD and MIT. You can override them at runtime from a JSON file via:

- C++: `load_cfd_tables("path/to/cfd.json");`
- Python: `import crispr_gpu as cg; cg.load_cfd_tables("path/to/cfd.json")`
- CLI: add `--cfd-table path/to/cfd.json` to `crispr-gpu score`.

## JSON schema

```
{
  "mm_scores": {"AG": 0.48, "AC": 0.97, ...},   // guide base + genome base
  "position": [0.61, 0.61, ... up to guide_len],  // positional multipliers
  "mit_position_penalty": [0.1, 0.1, ... up to guide_len] // optional
}
```

- `mm_scores` keys are two-letter strings (guide base, genome base) using A/C/G/T.
- `position` and `mit_position_penalty` arrays are truncated/extended to the guide length in use.
- Missing entries fall back to built-in defaults.

## Example

A minimal file overriding two mismatch types and the first two positions:

```
{
  "mm_scores": {"AG": 0.40, "CT": 0.65},
  "position": [0.70, 0.72]
}
```

## Notes
- JSON is loaded on the host and propagated to GPU runs automatically.
- Loading new tables is thread-unsafe; call it during initialization before launching searches.
- Guide-length specific arrays beyond the supplied length reuse the default values.
