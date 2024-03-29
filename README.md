# codon-usage

[![DOI](https://zenodo.org/badge/128677597.svg)](https://zenodo.org/badge/latestdoi/128677597)

This collection includes a number of functions and models that I have been using repeatedly in many sequence analysis and engineering projects, and are often not available out-of-the-box in any of the standard tools. They are organized into the following categories:

## models

* _CUFS_: The Codon Usage Frequency Similarity (Diament, Pinter & Tuller, 2014) defines a distance metric between gene pairs based on their coding sequence, and specificlly their codon and amino acid usage. Files: `calc_CUFS`, `calc_synCUFS`
* _CAI_: The Codon Adaptation Index (Sharp & Li, 1987) defines an optimality score for codons based on their frequency of appearance in a given reference set of genes. This is one of the most widely-used model for estimating the translational efficiency of a coding sequence. The weights assigned to codons can also be used to generate vectors for estimating the local optimality of the sequence. Files: `calc_CAI_weights`, `calc_score_from_weights`, `calc_vec_from_weights`
* _tAI_: The tRNA Adaptation Index (dos Reis, Savvy & Wernisch, 2004) defines an optimality score for codons based on a simple biophysical model of codon recognition during the translation process. Files: `calc_tAI_weights`, `calc_score_from_weights`, `calc_vec_from_weights`

## next-generation sequencing (ngs)

* _Resmpling_: Read count data is heavily dependent on the number of samples that were available when the dataset was generated. Downsampling (e.g., when comparing two experiments with widely different sample sizes), as well as sampling from a theoretical / simulated distribution, are often used in analysis and result validation. Files: `resample_reads`, `resample_profiles`

## optimize

*  _Codon usage bias_: The most popular methods for coding sequence optimization involve the selection of the most optimal codon to encode each amino acid in the protein according to some model. The converse can also be easily applied. A 'balanced', randomly generated profile based on a desired codon distribution may also be useful at times. Files: `maximize_CUB`, `minimize_CUB`, `.random/randseq_CUB`
* _Combinatorical sequences_: More elaborate optimization methods search the space of all possible sequences based on some objective function. The following functions are dedicated to the exhaustive  generation of synonymous sequences (producing the same protein), non-synonymous sequences, and their combinations. Files: `all_synonym_seq`, `all_nt_seq`, `all_combined_seq`

## random

* _Sampling synonymous sequences_ from a given codon distribution can be very useful (see above) and has been surprisingly missing from matlab's bioinformatics toolbox, so here it is. Files: `randseq_CUB`
* _Permuting sequences synonymously_ allows one to preserve the exact same codon usage while generting many samples. Codon permutations can be performed within each gene independently - changing their order while preserving the codon composition, and any optimality score that is determined by this composition of the gene - or globally, where codons may be exchanged between genes. Files: `shuffle_codons`
* _Generalized feature permutation_: Permutations of features along a vector based on their similarity / equivalence can be used as a control test. The following function returns a permutation for any vector of features with some distance defined between the features. For example, in the above scenario multiple codons that encode the same amino acid would have a small defined distance between them. But the approach here can further be extended to amino acids with similar biophysical properties, to features other than coding sequences, and so on. Files: `randperm_conserv_feat`
* _Generalized feature sampling_: Similarly, we can generalize the sampling approach to arbitrary features. Files: `rand_conserv_feat`

