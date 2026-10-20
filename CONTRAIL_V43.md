# TASOPT used by the four-lambda v4.3 campaign

This preservation branch starts at upstream `main` commit
`884389cd` (fetched 2026-10-20). Its `src/`, `Project.toml`, and
`Manifest.toml` are exact copies of the frozen TASOPT package used by jobs
668712, 668713, 668714, and 668715 (lambda 90, 180, 30, and 45).

Source snapshot on Hex:
`/home/aditeya/Contrail-Matrices/artifacts/v43/n50000-four-lambdas-20261019/TASOPT.jl`

The snapshot is authoritative for reproducing those jobs. Consequently this
commit includes differences from newer upstream main beyond the deviation
feature: silently retaining newer upstream scientific changes would not
reproduce the running campaign. Documentation, examples, and tests inherited
from upstream main are not a frozen campaign regression suite.

`contrail-v43-source-sha256.json` records every preserved package file.
The original campaign uses Julia 1.12.2 on Linux. The manifest is retained
even though the usual repository policy ignores it; its TASOPT entry has
the relative path `"."`. Native dependencies must be instantiated for the
destination platform. Linux Julia depots and Python environments should
not be copied as macOS installations.

This branch preserves the TASOPT dependency, not the complete optimization
campaign. The Contrail-Matrices package must separately preserve the frozen
`code/` and `workflow/`, shared weather/sample matrices, configs, warm-start
aircraft, and saved results. The campaign's Julia driver also invokes Python
audits. Its absolute paths and worker count need adjustment for local use.

The user's existing `Jenv_with_deviation` checkout and all running jobs are
unmodified by preparation of this branch.
