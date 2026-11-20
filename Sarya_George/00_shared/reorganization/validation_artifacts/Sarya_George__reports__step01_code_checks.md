# Step 1 code checks

10 focused checks passed.

- PASS: Negative TSB-vs-SN means higher in SN
- PASS: Positive stk-vs-anc means higher in stk
- PASS: The same medium response in both genotypes gives zero interaction
- PASS: A larger rise in SN in the ancestor gives a positive interaction
- PASS: The stp interaction uses stp and ancestor in the correct order
- PASS: The stp-vs-stk interaction preserves that mutant ordering
- PASS: Pooled contrasts are not assigned assumed weights
- PASS: Shortened pulse names recover the canonical sample ID
- PASS: Medium labels distinguish SNpulse from raw pulseSN spelling
- PASS: Reinjection flags follow the acquisition filename

These synthetic checks use no input measurements. The original testthat suite is separate.
