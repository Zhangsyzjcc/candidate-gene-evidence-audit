# Input schema

One row represents one claim in one cohort. Required fields are `layer`,
`claim_id`, `cohort`, `cohort_role`, `same_estimand`, `direction`, `p_value`,
`q_value`, `leading_edge_overlap`, `sensitivity_pass`, `measurement_adequate`,
and `notes`. See the anonymous example for valid encodings. Use `NA` for an
undefined numeric value; never replace missing statistics with zero.
