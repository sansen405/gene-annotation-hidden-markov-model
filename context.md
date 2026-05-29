1. RF as a second-stage splice-site rescorer (highest payoff)
Right now DONOR_1/2/3 and ACCEPTOR_1/2/3 get a single CNN log-odds per GT/AG candidate, and that's it. The CNN alone is producing too many confident false positives. Instead of replacing it, stack an RF on top that decides "is this GT/AG a real splice site?" using features the CNN doesn't see directly:

the CNN donor/acceptor score itself (your strongest single feature)
local sequence composition (GC%, k-mer counts in ±20–60 bp windows)
distance to nearest in-frame stop, distance to last predicted donor/acceptor
PSSM score of the splice window (you already compute PSSMs)
polypyrimidine tract / branch-point signal strength upstream of acceptors
reading-frame / phase context
The RF outputs a calibrated probability p, and you feed log(p/(1-p)) into Emission_Model exactly where the CNN log-odds goes today. Because the HMM only needs a log-odds number at each canonical candidate, this is a near drop-in swap of the emission value — your sparse-scoring design (only GT/AG rows) already matches the "one tabular example per candidate" shape an RF wants.

This directly attacks the precision collapse: the RF learns to suppress the CNN's confident-but-wrong sites using context the CNN ignores.

2. RF-based recalibration instead of the manual scale/bias sweep
Your results file literally says "bias values carried over from old BCE model — calibration sweep in progress" and you have a hand-tuned --tune-cnn-calibration grid over donor_scale/bias, acceptor_scale/bias. A small RF (or even gradient-boosted trees) trained on {CNN score → is_true_site} is a strictly more flexible calibrator than a single scale+bias affine transform, and it removes the stale-bias problem you're fighting right now. This is a lighter-weight version of #1 if you only want to use the CNN score as the feature.

3. RF as a predicted-gene / exon post-filter
After Viterbi, you predict 3,993 gene intervals vs 3,694 gold. An RF can act as a final filter on each predicted gene/exon, classifying real vs spurious from:

gene/exon/intron lengths and counts
mean Forward-Backward posterior confidence (you already compute this per base)
ORF integrity, GC content, codon-usage scores
number of marginal-quality splice sites in the gene
You drop or flag low-probability predictions. This cleans up precision after decoding without touching the HMM internals.