# SCAN-B sequencing-platform sensitivity before association fitting

Original analysis plan SHA256: `f174e513f1c3091a19d0ab3233512048c64e0b02f564a46101fd324ab550a7ae`.

The source clinical metadata identifies126 eligible primary samples sequenced on HiSeq2000 and17 on NextSeq500. The cohort analyst confirmed that the expression download/model implementation were still in progress and no IFI44L value or expression-outcome association had been inspected or fitted when this amendment was made.

Retain all original seven Cox specifications. Add two platform-stratified sensitivity models on the common143 age/IFN/IFI44L/outcome complete cases: age/10+IFN-alpha+strata(instrument) and age/10+IFN-alpha+IFI44L+strata(instrument). This permits different baseline hazards without estimating a potentially unstable coefficient for the small instrument group. Report platform-specific case/event counts, all coefficients and diagnostics; a stratum with no events remains explicitly visible and does not generate event information.

Add the one-df LRT for this nested pair to the two previously fixed LRTs. Use a single fixed three-test BH family for this additional cohort; preserve all nominal P values. Thus the cohort extension has nine Cox models and three exploratory LRTs. No specification is selected by significance, and this amendment changes neither receptor eligibility, originalF/replicate handling, expression transformation nor any primary-three-cohort hypothesis family.
