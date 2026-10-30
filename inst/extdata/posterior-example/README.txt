Synthetic external GPCM posterior example

This example demonstrates reading an exported summary, retaining fixed anchors,
and connecting credible intervals and sampling diagnostics to figures and reports.
It uses synthetic rating data and two short chains (300 retained draws each);
the sampling screen correctly requires review. The numbers are not findings
about a real assessment and do not qualify posterior recovery or coverage.

model-data.rds contains the response model, priors, constraints and synthetic
rating design. posterior-summary.rds is bound to those exact model-data bytes
and contains marginal summaries, not joint draws. No Stan/brms installation
or sampling is needed to read it with read_mfrm_posterior().
