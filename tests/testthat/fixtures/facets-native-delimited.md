# Native FACETS delimited Scorefiles

`facets-native-comma-score.txt` and `facets-native-semicolon-score.txt` are
unchanged Person exports from FACETS 4.5.0 (64 bit), executed on 2026-10-24.
Executable SHA256:
`dfb0afb0faa18f026d1b3b4175f22e42cc3764430eb83cbd368c7a572b3593a1`.

| File | SHA256 |
| --- | --- |
| comma | `706ce9226623fb96606ea6e9bd6895b4a7448a1b94b1f550a2abc6660378164d` |
| semicolon | `3b4223f596e78a597246dcbe2429281b29f47efcf1baf1eb6ce5e4b0c43113bd` |

The binary model reuses the synthetic weighted-anchor responses and saved
all-fixed calibration underlying the bundled FACETS score example. All Person,
Rater and Task locations are fixed; no new responses were simulated. The two
controls differ only in `CSV=Commas` / `CSV=Semicolons`. Both request
`Heading lines=Yes`, `QM=Double`, `Udecimals=6` and `Scorefile=score.txt`.
All 28 source fields agree for each record across formats. The native leading
facet line and dangling final header quote are retained exactly.

Question marks are intentional ASCII replacements in the input controls,
not evidence of Japanese labels surviving export. Duplicate labels and the
numeric-looking label `003` exercise element identity separately from labels.
The parser tests with Japanese UTF-8/CP932 text remain constructed input tests.
These fixtures qualify descriptive text import, not calibration or inference.
