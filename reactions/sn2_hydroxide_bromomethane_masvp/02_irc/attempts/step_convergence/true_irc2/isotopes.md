# Isotope check: deuterium on all four H against normal masses, same settings

Beyond the first step, each deuterium frame against the normal-mass path at the same key bond (Å).

| Path | Side | frames | RMSD median | RMSD largest | RMSD last | largest atom shift | other key bond | first step vs normal / D4 TS mode (°) |
|---|---|---|---|---|---|---|---|---|
| ORCA (Cartesian steepest descent) forward | reactant | 29 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.20 / 1.76 |
| ORCA (Cartesian steepest descent) backward | product | 29 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.20 / 1.76 |
| ORCA, D masses also in the Hessian file forward | reactant | 29 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.20 / 1.76 |
| ORCA, D masses also in the Hessian file backward | product | 29 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.20 / 1.76 |
| pysisyphus EulerPC (mass-weighted IRC) forward | reactant | 30 | 0.0027 | 0.0139 | 0.0139 | 0.0446 | 0.0021 | 1.83 / 0.18 |
| pysisyphus EulerPC (mass-weighted IRC) backward | product | 30 | 0.0028 | 0.0049 | 0.0049 | 0.0103 | 0.0008 | 1.83 / 0.18 |
