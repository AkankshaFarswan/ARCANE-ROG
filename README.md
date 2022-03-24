# ARCANE-ROG

Step 1: Run denoise_RGL_2times- It will denoise and impute the noisy and incomplete matrix.


Input: Noisy and incomplete matrix


Output: Denoised and complete matrix, Adjacency graph


Step 2: Run_leiden_RGL.R- It will generate the clonal tree using the adjacency matrix and denoised matrix learnt at step 1


Citation policy:
Please cite the following paper if you use the code in your analysis-
Akanksha Farswan, Ritu Gupta, and Anubha Gupta, "ARCANE-ROG: Algorithm for Reconstruction of Cancer Evolution from single-cell data using Robust Graph Learning," Journal of Biomedical Informatics, Elsevier, March 03, 2022.

