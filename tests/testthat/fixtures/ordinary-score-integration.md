# Ordinary retained-prior scoring regression

`ordinary-score-integration.rds` reuses the synthetic fixed-N(0,1) PCM-MML
and ordinary PCM-JML fits from the October 2 A/B comparison. Calibration has
20 Persons, Rater 3 x Criterion 3, four categories and six responses per
Person. The independent, balanced held-out panel has 60 Persons and two
events per assigned cell (12 responses each). It is the retained F-exposure
panel, not a new calibration dataset or a new independent study replicate.

Before the repair, NULL-prior scoring returned 31-node EAP/SD without a batch
integration check. The explicit, identical N(0,1) prior required 121 scoring
nodes on this panel. Maximum EAP/SD differences were .005816587/.008506647
for PCM-MML and .001878822/.003752990 for PCM-JML. These numerical differences
motivate a regression for check parity, not a statistical accuracy claim.

Original calibration phase checksums, held-out panel fingerprint and execution
source hashes are stored in the fixture. The test holds calibration and rows
fixed, checks refusal/review behavior at 31 nodes, and compares retained and
explicit identical priors at 121 nodes. No fitting or data generation is needed.
