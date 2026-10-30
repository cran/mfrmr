# Private Wald result fixtures

These three files contain the `result` objects from the source-bound
`validation-results/gmfrm-result-20261011/qualified/` epoch:

- `gmfrm-wald-result.rds`: `live_binary.rds`, the deterministic N20 binary witness.
- `gmfrm-wald-N20-observed.null.rds`: the original constrained N20 result.
- `gmfrm-wald-N20-b09.free.rds`: the original unselected boundary result.

Extraction verified the preceding replay's result-file hashes. No point estimate,
stage, matrix, admission flag or fingerprint was altered. The new tests check
full-information/contrast calculations, log-scale endpoints and refusals without
running an optimizer or regenerating responses. The independent Louis kernel is
also compared with the existing scalar reference and a directional derivative
at nonstationary binary/four-category points. The current validation manifest is
`validation-results/gmfrm-wald-20261011/manifest.rds`.

These numerical fixtures are not new independent sampling-performance evidence.
