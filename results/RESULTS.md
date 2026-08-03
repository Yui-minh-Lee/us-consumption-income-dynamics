# Archived Results

These tables reproduce the estimates reported in the completed course paper so that a reviewer can inspect the main evidence without running MATLAB.

## Baseline single-cycle model

| Parameter | Estimate | Std. error | p-value |
|---|---:|---:|---:|
| Trend volatility | 0.0109 | 0.0010 | <0.001 |
| Cycle damping | 0.9904 | 0.0054 | <0.001 |
| Cycle period (months) | 95.9999 | 8.6520 | <0.001 |
| Cycle volatility | 0.0099 | 0.0012 | <0.001 |
| Trend-cycle correlation | -0.8671 | 0.0275 | <0.001 |
| Consumption sensitivity | 1.0411 | 0.0228 | <0.001 |
| State dependence | -0.0107 | 0.0284 | 0.7056 |
| Okun coefficient | -0.3141 | 0.0362 | <0.001 |
| Transfer loading | 0.7161 | 0.0277 | <0.001 |

## Extended dual-cycle model

| Component | Parameter | Estimate | Std. error | p-value |
|---|---|---:|---:|---:|
| Long | Damping | 0.9838 | 0.0069 | <0.001 |
| Long | Period (months) | 86.4251 | 7.9597 | <0.001 |
| Long | Volatility | 0.0054 | 0.0006 | <0.001 |
| Long | Consumption sensitivity | 1.0311 | 0.0477 | <0.001 |
| Long | Okun coefficient | -0.3680 | 0.0250 | <0.001 |
| Short | Damping | 0.8495 | 0.0287 | <0.001 |
| Short | Period (months) | 24.0000 | - | Bound |
| Short | Volatility | 0.0080 | 0.0004 | <0.001 |
| Short | Consumption sensitivity | 0.2964 | 0.0616 | <0.001 |
| Short | Okun coefficient | 0.0338 | 0.0344 | 0.3262 |
| Other | Consumption transfer loading | 0.5412 | 0.0354 | <0.001 |
| Other | Unemployment transfer loading | 0.0602 | 0.0189 | 0.0015 |

## Robustness summary

- Alternative period bounds leave long-cycle sensitivity (1.027) above short-cycle sensitivity (0.338).
- Relaxing measurement-noise restrictions increases short-cycle sensitivity, but it remains below long-cycle sensitivity.
- Adding a noise floor produces sensitivities of 0.972 (long) and 0.564 (short).
- The short-cycle Okun coefficient is sensitive to measurement-noise restrictions, so its structural interpretation is deliberately limited.
