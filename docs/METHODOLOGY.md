# Methodology

## Economic question

Does aggregate consumption respond differently to persistent business-cycle movements and shorter-lived income fluctuations?

## Data preparation

Nominal consumption and disposable personal income are deflated by PCEPI. Market income is measured using real personal income excluding current transfer receipts. A transfer factor is constructed as the log difference between real disposable income and market income.

The scripts use Hamilton filtering to remove slow-moving components of the consumption-income ratio and transfer factor. The unemployment rate is decomposed into trend and gap components and enters the state-space system as an external labor-market anchor.

## Baseline model

```text
y_t     = tau_t + psi_t
c_t     = tau_t + lambda_t psi_t + gamma Tr_t + u_c,t
u_gap,t = theta psi_t + u_u,t
```

The latent cycle follows a damped stochastic cycle. The baseline permits `lambda_t = lambda_0 + beta S_(t-1)`. In the archived estimate, beta is not statistically significant.

## Dual-cycle extension

```text
y_t     = psi_L,t + psi_S,t
c_t     = lambda_L psi_L,t + lambda_S psi_S,t + gamma_c Tr_t + u_c,t
u_gap,t = theta_L psi_L,t + theta_S psi_S,t + gamma_u Tr_t + u_u,t
```

Each component follows a damped rotation with constrained period bounds. The reported specification targets 60-96 months for the long cycle and 8-24 months for the short cycle.

## Estimation

Parameters are mapped from an unconstrained optimization space into economically admissible ranges. `fminsearch` supplies a simplex pre-optimization; MATLAB's state-space `estimate` routine and `fminunc` then refine the solution and provide an approximate Hessian.

## Interpretation boundaries

- Extracted cycles are statistical latent components, not independently identified structural shocks.
- An insignificant Okun coefficient should not be assigned a structural-shock label.
- The transfer factor is endogenous to macroeconomic conditions; its coefficient is conditional association.
- Standard errors rely on local curvature in a constrained nonlinear model.

