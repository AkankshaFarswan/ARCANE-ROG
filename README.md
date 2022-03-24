# ARCANE-ROG

Step 1: Run denoise_RGL_2times- It will denoise and impute the noisy and incomplete matrix.


Input: Noisy and incomplete matrix


Output: Denoised and complete matrix, Adjacency graph


Step 2: Run_leiden_RGL.R- It will generate the clonal tree using the adjacency matrix and denoised matrix learnt at step 1
