# RMMC — Robust Multi-Matrix Completion via Hybrid Ordinary–Welsch Loss

MATLAB implementation, the six baselines, the S1–S4 experiment programs and the
ORL dataset, for

> **Robust Multi-Matrix Completion via Hybrid Ordinary–Welsch Loss**
> Hankuan Gao, Hao Nan Sheng, Hing Cheung So, Zhiyong Wang
> submitted to ICASSP 2027

*Multi-matrix completion* (MMC) restores a collection of matrices that share a row and a
column subspace and whose missing entries are **structural** (whole rows or columns can be
absent), by combining the shared low-rank factors with a first-difference smoothness
penalty.  Its fidelity term is quadratic, so a few corrupted observations can dominate the
shared factors.  **RMMC** keeps MMC's shared subspaces and smoothness regularisation and
replaces the quadratic fidelity by the **hybrid ordinary–Welsch (HOW)** loss.  The
half-quadratic form of the HOW loss turns the outlier into a per-matrix variable that
absorbs a corrupted entry instead of fitting it, so a collection can be recovered when
rows/columns are missing *and* the observed entries are corrupted.  Under explicit penalty
and nondegeneracy conditions the iterates are bounded and every accumulation point is a
critical point of the limiting-threshold problem.

## Repository layout

```
<repo root>/
├── RMMC.m                 the proposed solver (Algorithm 1 of the paper)
├── rmmc_defaults.m        λ and (ρ₁,ρ₂) per data family — the values used in the paper
├── MMC.m                  reference MMC implementation (released code)
├── MMC_half.m             the MMC row actually run in the comparison
├── utils/
│   └── GLRAM.m            the warm-started subspace sweeps used by both solvers
├── baselines/
│   ├── RMC_HOW/           HOW_RMC.m, L0_IQR.m
│   ├── RTCUR/             upstream RTCUR release + tools/ + LICENSE   (MIT, third party)
│   ├── TriTD/             upstream triple-decomposition release      (see PROVENANCE.txt)
│   └── RGNMR/             upstream RGNMR release + LICENSE           (MIT, third party)
├── data/
│   └── orl/orl_face/      ORL face database — the S4 input (40 subjects × 10 images)
├── third_party/
│   └── tensor_toolbox/    Tensor Toolbox v3.1 for MATLAB (BSD-style, LICENSE.txt)
├── experiments/
│   ├── main_run.m         one (scenario, setting, seed) cell, every method
│   ├── main_plan.m        build the missing-work schedule
│   ├── main_worker.m      run the cells assigned to one worker
│   ├── main_report.m      aggregate cells/ into tables/
│   ├── main_timing.m      the dedicated single-core runtime pass
│   └── run_main_6p.ps1    launcher: six pinned workers, then main_report
└── results/               everything the programs write (created on first run, ignored)
```

## Addressing

Every input and output path is derived from the location of the program, so a clone runs
without editing anything:

```
<repo root>/experiments/main_*.m   the programs
<repo root>/                       RMMC.m, MMC*.m, utils/, baselines/, third_party/
<repo root>/data/orl/orl_face/     the S4 input (ORL, 40 subjects × 10 images)
<repo root>/results/               created by main_run: schedule.csv, cells/, metrics/,
                                   traces/, mats/, timing/, logs/, tables/
```

`main_run.m` resolves `codeRoot` from its own filename, reads the S4 ground truth from
`data/orl`, puts the baselines and the Tensor Toolbox on the path, and writes into
`results/`.  `main_plan.m`, `main_worker.m` and `main_report.m` do the same.
`run_main_6p.ps1` derives the package root from `$PSScriptRoot` and takes the MATLAB
executable as `-MatlabExe`.

## The solver

Each RMMC round is

1. **robust block** (first in the round, because it depends only on the previous iterate):
   thresholds from the signed IQR of the *observed* residuals,
   `c_k^i = min(c_k^{i-1}, η₁ d_k^i)`, `σ_k^i = min(σ_k^{i-1}, η₂ d_k^i)`, then
   `O_k^i = p_{c,σ}(P_{Ω_k}(X̃_k − M_k^{i-1}))`;
2. **factor block**: three warm-started GLRAM sweeps on `A_k^i = M_k^{i-1} + G_k^{i-1}/ρ₁`;
3. **quadratic blocks**: `M`, `X`, `N` in closed form through the precomputed
   `H_L = (1+ρ₁+ρ₂)I + 2λLᵀL` and `H_R = ρ₂I + 2λRRᵀ`, then plain dual ascent on `G`, `W`.

Deterministic (no random numbers), driven only by `λ`, `ρ₁`, `ρ₂`, `η₁`, `η₂`, the ranks
and the iteration cap.  `rmmc_defaults.m` returns the paper's settings:
`(ρ₁,ρ₂) = (9.05, 1.45)` for both families, `λ = 5e-3` (synthetic) or `0.1` (ORL),
`(η₁,η₂) = (2.5, 3.5)`.

`MMC_half.m` is the same alternating scheme with the outlier channel removed and MMC's
observed-entry calibration kept — the MMC row of the tables.  `MMC.m` is the released
reference implementation, kept for comparison.

## Baselines

| Row in the tables | File | Provenance |
|---|---|---|
| MMC | `MMC_half.m` (`MMC.m` = released reference) | TPAMI 2026 |
| RMC-HOW | `baselines/RMC_HOW/HOW_RMC.m` | our transcription of the released HOW method, TSP 2023 |
| RMC-ℓ₀ | `baselines/RMC_HOW/L0_IQR.m` | our transcription of the ℓ₀-norm robust MC method, TAES 2023/2025 |
| RTCUR | `baselines/RTCUR/RTCUR_fc.m` | upstream release, MIT |
| TriTD | `baselines/TriTD/triple_decomp_ADMM.m` | upstream release, see `PROVENANCE.txt` |
| RGNMR | `baselines/RGNMR/RGNMR.m` | upstream release, MIT; ships prebuilt Windows mex |

## Running S1–S4

MATLAB R2025b.  From the repository root:

```matlab
addpath('experiments'); cd('experiments')

main_run('s1',   9,  5)      % S1 at 9 dB, seed 5, every method
main_run('s2', 0.5,  5)      % S2 at 50 % missing
main_run('s3',   9,  5)      % S3 at 9 dB
main_run('s4sp',10,  5)      % S4 at q = 10 %
```

| Scenario | Setting swept | Fixed |
|---|---|---|
| `s1` | SNR ∈ {0,3,6,9,12,15} dB | random missing R = 0.2 |
| `s2` | R ∈ {0, 0.1, …, 0.9} | 10 dB |
| `s3` | SNR as S1, with 10 rows and 10 columns removed per frame | R = 0.2 |
| `s4sp` | salt-and-pepper q ∈ {5, 10, 15} % | R = 0.2 plus six structural bands |

Five seeds (5–9) per setting give the 125 cells of the paper.  The whole sweep is run by

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File experiments\run_main_6p.ps1
```

which plans the missing work, runs it on six pinned single-threaded workers and finally
calls `main_report`.  The dedicated runtime pass is `main_timing(1..4)` (one call per
scenario, one process per physical core); it writes `results/timing/`.

`main_run` writes, per cell: one row in `results/cells/<scenario>/`, per-frame metrics in
`metrics/`, the per-iteration trace in `traces/`, the full float reconstruction `Xh_<M>` in
`mats/<scenario>/<setting>_seed<k>.mat`, and the worker log.  A cell already present in the
CSV is skipped, so a sweep can be interrupted and resumed — a fresh clone has no
`results/`, so it runs all 125 cells from scratch.

Settings used for every number in the paper:

| | synthetic (S1–S3) | ORL faces (S4) |
|---|---|---|
| data | 100 × 200 × 20, Tucker rank (20,20,20) | 400 images, 112 × 92 |
| λ | 5 × 10⁻³ | 0.1 |
| (ρ₁, ρ₂) | (9.05, 1.45) | (9.05, 1.45) |
| (η₁, η₂) | (2.5, 3.5) | (2.5, 3.5) |
| rank | r = s = 20 | 90 |
| iteration cap | 200 | 400 |

## Plotting

Read the reconstructions from `results/mats/<scenario>/<setting>_seed<k>.mat` (`X0`, `In`,
`Om`, `Xh_<METHOD>`) or the per-iteration CSVs in `results/traces/`, then plot what you
need.

Note on `traces/`: the column meaning differs per method — `relerr` for RMMC and MMC;
`err` (an internal sub-block residual) for RTCUR; `residual` (an ADMM residual) for TriTD;
per-frame records for HOW, L0IQR and RGNMR.

## Licence

Our code is MIT (`LICENSE`).  The package also redistributes third-party research code
with its original licence; see *Third-party notices* in `LICENSE` and the
`PROVENANCE`/`LICENSE` files next to each component.  The ORL face database
(`data/orl/`) is redistributed for research use with attribution to AT&T Laboratories
Cambridge.

## Citation

```bibtex
@inproceedings{Gao2027RMMC,
  title     = {Robust Multi-Matrix Completion via Hybrid Ordinary--Welsch Loss},
  author    = {Gao, Hankuan and Sheng, Hao Nan and So, Hing Cheung and Wang, Zhiyong},
  booktitle = {IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
  year      = {2027}
}
```
