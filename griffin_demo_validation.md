Merge fix/feature-window-fencepost

Feature windows now match Griffin's grid-snapped half-open construction
(132/4/128 bins for save/center/fft). Validated against Griffin's published
demo output: profile r=0.9904, features within 4%.

KNOWN: the stored golden artifacts predate this fix and encode the old
closed-interval windows. golden_test.sh will FAIL on mean_coverage,
central_coverage and amplitude at --feature-tol 0. This is expected.
Do NOT resolve it by loosening the tolerance. Regenerate the artifacts on
LeoMed as a separate commit.
