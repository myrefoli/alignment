# alignment
A Python implementation of the Smith-Waterman and Needleman-Wunsch algorithms.

We often want to compare the DNA or protein sequences from different species to examine evolutionary relationships or to look at shared functionality. In order to make these comparisons, we need to examine which portions of these sequences align. Smith-Waterman and Needleman-Wunsch algorithms are classical sequence alignment algorithms using dynamic programming.

These algorithms are implemented using the affine gap penalty and can be used for either global ends-free or local alignment (the input file has an indicator on which alignment it wants).
