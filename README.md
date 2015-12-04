# Ensemble SNV Calling

The included script, `ensemble_caller.py`, implements the recommendation described and supported in Ewing at al. (doi: 10.1038/nmeth.3407), namely to perform ensemble single nucleotide variant (SNV) calling with a voting scheme. Specifically, given sets of SNV calls from various algorithms, this script will match up corresponding variants and output SNV calls that are supported by a majority of algorithms.
