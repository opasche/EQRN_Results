Keywords for EQRN model:
- `bestval`: best EQRN model
- `bestval_nopen`: best unpenalized EQRN model (L2 weight penalty = 0)

Keywords for competitor intermediate model (separate grid-search has been performed in each case for gbex):
- `comprnn`: competitors use QRN as intermediate quantile model (same as EQRN, fairest comparison)
- `compgrf`: competitors use GRF as intermediate quantile model (most consistent/likely choice in practice)
- `comporacle`: competitors use true quantiles as intermediate quantiles (gives them a strong initial advantage)

