# JSON Import Branch Status

Findings from comparing `json_import` with `main`.

1. **Grouped imports lose the grouping-column name**

   `R/import.R:657` sets `group.var` to `TRUE` instead of restoring the original
   grouping-column name. The imported object therefore does not preserve the
   grouped model metadata exactly.

2. **Bayesian imports do not reconstruct native Bayesian fields**

   Bayesian imports preserve `original_model` as parsed JSON, but do not rebuild
   native fit fields such as `marginals`, `marginal.quantiles`, `mlik`, or
   `used.INLA`.

3. **`include_network` is ignored during export**

   The argument is documented and passed through the export functions, but it
   does not affect the output. `variables` and `arcs` are always exported.

4. **Fitted-model JSON validation is shallow**

   Validation checks required top-level fields and link-function presence, but
   does not enforce all documented uniqueness, reference, state, distribution,
   or DAG validity rules.

5. **Raw-data JSON validation is narrower than documented**

   Validation does not currently enforce all rules described in the vignette,
   including category membership, count-value validity, role/distribution
   consistency, and required metadata fields.

6. **Raw-data imports can mishandle categorical columns without R adapters**

   When R-specific adapter metadata is absent, categorical label columns can be
   reconstructed incorrectly, especially when category metadata is missing.

7. **Special numeric values are not encoded consistently with metadata**

   `NaN` and infinite values are not explicitly serialized according to the
   declared `missing_values` metadata. In particular, `data_json_scalar_value()`
   treats `NaN` as missing.

8. **The branch contains trailing whitespace**

   `git diff --check` reports trailing whitespace in newly added implementation,
   vignette, and test lines.

## Scope and Test Coverage

The branch has broad tests for standard MLE, grouped MLE, multinomial states,
Bayesian export metadata, JSON determinism, name irrelevance, raw-data round
trips, and summary repair.

The main remaining test gaps correspond to findings 1 through 7, especially
malformed JSON semantics, `include_network`, grouped metadata preservation,
Bayesian native-object reconstruction, and portability without R adapter
metadata.
