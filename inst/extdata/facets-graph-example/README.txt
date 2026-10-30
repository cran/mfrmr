Synthetic FACETS Graphfile example (FACETS 4.5.0)

graph.txt is an unmodified native CSV Graphfile from a small four-scale PCM
analysis. Each of 4,004 rows is one reported measure-grid coordinate on a
Table 8 scale; Prob:0 through Prob:3 give category probabilities. Expected
and ExpCat are distinct exported fields.

All six Person measures, four Task measures and four sets of thresholds are
fixed in input.txt. There are 21 synthetic observed ratings with unequal
numbers of observations per Person. This is a format/plotting example, not
evidence of model fit, estimation recovery, interval coverage or fairness.
It is a different analysis from extdata/facets-example (the binary example).

To reproduce with a licensed FACETS installation, copy input.txt into an empty
directory and run: Facets.exe BATCH=YES input.txt output.txt
The specification requests named CSV Graphfile output with headings. Keep the
source categories and coordinates; scale numbers do not themselves identify
items when importing an arbitrary analysis.
