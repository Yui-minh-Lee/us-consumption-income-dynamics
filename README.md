# Frequency-Dependent U.S. Consumption-Income Dynamics

> A MATLAB state-space study of how U.S. consumption co-moves with market income at long and short business-cycle frequencies, with unemployment and fiscal transfers as additional signals.

## 项目简介

本项目研究美国居民消费如何响应不同频率的收入波动。传统总量模型往往用单一消费敏感度概括所有冲击；本项目使用带失业率锚定的状态空间模型，将市场收入分解为长、短周期，并分别估计消费对两类周期的响应。

研究使用美国月度宏观数据、Hamilton filter、HP filter、不可观测成分模型和约束极大似然估计，展示从经济问题、数据处理、模型识别到稳健性检验和结果可视化的完整实证流程。

## Main findings

| Component | Estimated period | Consumption sensitivity | Okun coefficient |
|---|---:|---:|---:|
| Long cycle | 86.4 months | 1.031*** | -0.368*** |
| Short cycle | 24.0 months | 0.296*** | 0.034 (not significant) |

- Long-cycle consumption-income co-movement is materially stronger than short-cycle co-movement.
- The baseline state-dependence coefficient is not statistically significant (`p = 0.706`), motivating a frequency decomposition rather than a crisis-dummy interpretation.
- The transfer-factor loading is positive (`gamma = 0.541`) in the dual-cycle model. It is a conditional association, not a causal fiscal multiplier.
- The short-cycle Okun coefficient is not significant in the main specification, so the short component is not assigned a definitive structural-shock label.

Full archived estimates and robustness checks are reported in [`results/RESULTS.md`](results/RESULTS.md).

![Long- and short-cycle extraction](results/figures/Fig1_Cycles.png)

![Consumption-cycle decomposition](results/figures/Fig4_Cons_Decomp.png)

## Research workflow

1. Construct real consumption, real disposable income, market income, and a transfer factor.
2. Inspect the consumption-income relation and test adjusted residual stationarity.
3. Estimate a trivariate baseline state-space model using income, consumption, and the unemployment gap.
4. Extend the model to separate long and short stochastic cycles.
5. Compare alternative period bounds and measurement-noise restrictions.
6. Export recruiter-readable figures and parameter tables.

See [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md) for the model equations and interpretation.

## Repository structure

```text
data/
  US_Data.csv              Monthly source dataset
src/
  prepare_data.m           Transformations and stationarity check
  baseline_model.m         Single-cycle trivariate model
  dual_cycle_model.m       Long/short-cycle extended model
  export_figures.m         Publication-quality figure export
results/
  figures/                 Selected 300-DPI outputs
  RESULTS.md               Estimates and robustness summary
docs/
  DATA_DICTIONARY.md       FRED series, units, and transformations
  METHODOLOGY.md           Model specification and limitations
```

Development notes, downloaded papers, autosave files, and superseded model variants are retained locally but excluded from the public repository.

## Reproduction

Tested for static compatibility with MATLAB R2024a. Required products are MATLAB, Econometrics Toolbox, Optimization Toolbox, and Statistics and Machine Learning Toolbox.

Open MATLAB at the repository root and run:

```matlab
run('src/prepare_data.m')
run('src/baseline_model.m')      % baseline estimates
run('src/dual_cycle_model.m')    % extended estimates
run('src/export_figures.m')      % run after the extended model
```

The estimation scripts use `fminsearch` followed by `fminunc` and can take several minutes. `prepared_data.mat` is generated under `data/processed/` and is not version-controlled.

## Data and sample

The source CSV contains monthly observations from January 1959 to September 2025. To reproduce the archived analysis while excluding the pandemic regime, the scripts select the first 733 observations (January 1959 through January 2020). The one-month difference from a strict December 2019 cutoff is documented rather than silently changed because the archived figures use this specification.

The model uses public FRED series including DSPI, PCE, PCEPI, W875RX1, and UNRATE. See the [data dictionary](docs/DATA_DICTIONARY.md) for definitions and source links.

## Limitations

- The cycle decomposition is model-dependent and is not structural causal identification.
- Fiscal transfers are endogenous to economic conditions; the estimated loading is not a causal policy multiplier.
- The exact historical FRED vintage/retrieval date was not preserved in the original workflow.
- Results are an in-sample macroeconomic decomposition, not an asset-return forecasting strategy.

## Disclaimer

This repository is a research and portfolio project. It is not investment advice. Any errors are my own.
