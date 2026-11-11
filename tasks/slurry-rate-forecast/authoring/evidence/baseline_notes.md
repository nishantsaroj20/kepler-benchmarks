Every naive answer was run through the verifier image on the shipped truth. The reward was 0 in each case, and no charge passed both its rate test and its fraction test.

Constant 0, constant 1, the mean of the latent values, and the median of the latent values each failed 16 tests. A constant equal to one latent rate, copied into the rate and into all three fractions, failed 15 tests: the matching charge can meet the rate tolerance and still misses the fraction tolerances. Copying `bath_C`, `P_total_barg`, `N_rps`, `m_cat_g`, or `dp_um` into every numeric field failed 16 tests.

The do-nothing run, with `/app/output` empty, wrote reward 0. Pytest collected 19 tests in junitxml (19 failures, 0 errors). The CTRF file is the pytest-json-ctrf plugin output; that plugin merges the parametrized tests, so the count used here is the junit count.

Symlink, empty `cases`, a duplicate `case_id`, an extra case, NaN, and Infinity each wrote reward 0.

Per-baseline verifier output is summarized in `local_runs.json`. Per-case tolerance arithmetic is in `baseline_failures.json`.
