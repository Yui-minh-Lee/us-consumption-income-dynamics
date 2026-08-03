# Data Dictionary

The public CSV is a locally compiled monthly dataset. The final models use the five series below; other columns are retained for provenance but are not used by the portfolio scripts.

| CSV column | FRED series | Definition | Unit in source | Model transformation |
|---|---|---|---|---|
| `Income` | [DSPI](https://fred.stlouisfed.org/series/DSPI) | Disposable personal income | Billions of dollars, SAAR | `log(DSPI / PCEPI x 100)` |
| `Price` | [PCEPI](https://fred.stlouisfed.org/series/PCEPI) | Personal consumption expenditures price index | Index, 2017 = 100 | Deflator |
| `Consumption` | [PCE](https://fred.stlouisfed.org/series/PCE) | Personal consumption expenditures | Billions of dollars, SAAR | `log(PCE / PCEPI x 100)` |
| `Real_income_excptTransfer` | [W875RX1](https://fred.stlouisfed.org/series/W875RX1) | Real personal income excluding current transfer receipts | Billions of chained 2017 dollars, SAAR | Natural log |
| `Unemploy_rate` | [UNRATE](https://fred.stlouisfed.org/series/UNRATE) | Civilian unemployment rate | Percent, seasonally adjusted | Divide by 100; extract gap |

## Constructed variable

```text
transfer_factor = log(real disposable personal income)
                - log(real personal income excluding transfers)
```

Because transfers respond automatically to economic conditions, the associated coefficient is descriptive rather than a causal fiscal multiplier.

## Sample and provenance

- Raw coverage: 1959-01 through 2025-09.
- Archived estimation window: first 733 observations, 1959-01 through 2020-01.
- Frequency: monthly.
- The exact download date and historical FRED vintage were not recorded in the original workflow.
- Blank observations in unused auxiliary columns do not affect the final model inputs.

