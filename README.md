# SCW
SCW (Single-Cell Whole-genome) is a computational method to reconstruct high-resolution 3D whole-genome structures based on the zero-inflated single-cell Hi-C data.
## 1. Usage
Compilation command:

g++ -o scw -pthread -std=c++11 scw.cpp

Execution command:

./scw -i NXT.896 -res 0.5 -len mm9_mouse_chr_length -or chromosome_order -np 40 -o3 model_scl/coords_scl -o2 model_rearrangement/coords_reaarange -o model_refinement/coords_refine

-i NXT.896: Specifies the input file as NXT.896, which is a Hi-C file containing information about interactions between chromosomes.

-res 0.5: Sets the resolution parameter to 0.5 Mb (megabases). This parameter determines the granularity or scale at which the Hi-C data will be processed.

-len mm9_mouse_chr_length: Uses the mm9_mouse_chr_length file to provide the lengths of each chromosome. The file format consists of one chromosome per line in the format chrX N, where N represents the length of chromosome X.

-or chromosome_order: Orders the chromosomes based on the number of inter-chromosome Hi-C contacts. This option determines the order of the chromosomes based on the strength of interactions between them.

-np 40: Specifies the number of threads as 40. This is an optional parameter used to control the number of threads for parallel computation, with a default value of 40.

-o3 model_scl/coords_scl: Specifies the output file coords_scl in the model_scl directory to store the coordinates of the generated SCL model.

-o2 model_rearrangement/coords_reaarange: Specifies the output file coords_reaarange in the model_rearrangement directory to store the coordinates of the generated model during the rearrangement process.

-o model_refinement/coords_refine: Specifies the output file coords_refine in the model_refinement directory to store the coordinates of the generated model during the refinement process.
## 2. Hi-C contacts file

The Hi-C files for this project must follow the format specified in the provided sample file "NXT.896". The first and third columns in the file represent the chromosome names, the second and fourth columns represent the pair ends, and the fifth column represents the number of contacts.

Ensure that your Hi-C files adhere to this format for proper compatibility and accurate processing.

## 3. GAE_package

In the GAE_package folder, you can generate an imputed Hi-C matrix using the following instructions:

python gae.py --nodes 5289 --fin NXT896_500kb.edg --fout gae_imputed_matrix --lr 0.01

Explanation:

python gae.py: Runs the Python script gae.py to perform the imputation process using the GAE (Graph Autoencoder) package.

--nodes 5289: Specifies the number of nodes (or dimensions) for the GAE model.

--fin NXT896_500kb.edg: Specifies the input edge file NXT896_500kb.edg containing the Hi-C data. Please refer to the readme file in the folder for instructions on generating this edge file.

--fout gae_imputed_matrix: Specifies the output file name as gae_imputed_matrix to store the generated imputed Hi-C matrix.

--lr 0.01: Sets the learning rate parameter to 0.01, controlling the speed of the learning process for the GAE model.

Please refer to the readme file in the GAE_package folder for detailed instructions on how to generate the required edge file.
