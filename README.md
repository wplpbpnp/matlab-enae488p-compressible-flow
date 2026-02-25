# MATLAB ENAE488P Hypersonic Flow

MATLAB coursework and project scripts for ENAE488P hypersonic flow / high-temperature gas dynamics work, including a hypersonic caret waverider project.

## Course Summary (Testudo)

ENAE488P is a Special Topics in Aerospace Engineering course. The Testudo listing for the offering that matches this material is titled `Hypersonics and High-Temperature Gas Dynamics`, covering hypersonic flow, aerodynamic heating, aerothermodynamics, analysis methods, and engineering applications.

Testudo (topic offering reference): `https://app.testudo.umd.edu/soc/search?courseId=ENAE488P&sectionId=&termId=202408`

## Contents

- Homework scripts (`hw1.m` to `hw5.m`)
- Project/problem scripts (`final_project.m`, `tbm.m`, `oshock.m`, `equation.m`, `m2.m`, `p2p1.m`)
- Supporting images (`*.jpg`, `*.tif`)
- Published outputs in `html/` (`*.pdf`)

## Getting Started (MATLAB)

1. Open MATLAB in the repo root.
2. Start with the function pair `tbm.m` (theta-beta-mach solver) and `oshock.m` (oblique shock relations), which are reused by `final_project.m`.
3. Run smaller standalone scripts such as `equation.m`, `m2.m`, or `p2p1.m` to validate the environment and inspect intermediate methods.
4. Run `final_project.m` for the most complete project example (it calls `tbm(...)` and `oshock(...)` in a hypersonic caret waverider analysis workflow).
5. Review matching published outputs in `html/` (for example `html/oshock.pdf`, `html/tbm.pdf`, `html/final_project.pdf`) to compare expected report-style results.

## Dependencies / Compatibility Notes

- `tbm.m` uses Symbolic Math (`sym`, `vpasolve`); Symbolic Math Toolbox is required for workflows that call it (including `final_project.m`).
- Scripts assume helper functions remain on the MATLAB path (running from the repo root is the safest default).

## Highlighted Files

- `final_project.m`
- `oshock.m`
- `tbm.m`
- `html/final_project.pdf`

## Project Highlight

- `final_project.m`: hypersonic caret waverider analysis/design study using oblique-shock and theta-beta-mach calculations.
