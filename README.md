SCW has 5 versions aimed at different situations. 

We released two folders. One for one-step, the other for two-step.

In one-step, there are two subfolders. One for GAE, one for 2D Guassian.

In two-step folder, there are three subfolders. One for GAE, one for 2D Gaussian low resolution, and the other is 2D Gaussian high-resolution (recommended for resolutions larger than 100 Kbp, like 50 Kbp, or 20 Kbp). 

1. One step 2D Gaussian version is the simplest one. For running this, we need one Hi-C contact pair file, one chromosome length file.

One contact pair file is in the example_data file, which is GM12878 cell12. The format for Hi-C is "1_0 766645 1_0 16241940", where 1_0 means the paternal of chromosome 1, 1_1 means the maternal of chromosome 1. The format for length is "1_1 249250621", 1_1 is the maternal chromosome 1, and 249250621 is the reference length that matched your data.

"g++ -o so2d -pthread -std=c++11 one_step_2D_Gaussian.cpp";#so2d is the compiled file, you can name it anything you like.

"./so2d -i ../../example_data/gm12.hic -res 1 -len ../../example_data/human_list -np 40 -o out_g12_1mb";# -i followed by the Hi-C file, -res followed by the resolution, unit is Mbp, double value used here (0.5, 0.05, or 1, is the right style), -np followed by the thread numbers you will assigned, value type is integer, -o followed by the output file name.

2. One step GAE version is requested for one more parameter, which is the imputation file. We need to go to "GAE_package" in the main folder to generate that impuatation file.

To generate an imputation file, we need to use GAE to generate the imputation matrix, then convert the matrix file for SCW. 

For example, go to "GAE_package" folder, run "perl generate_edge.pl human_list 1 ../example_data/gm12.hic >gm12_edge" to generate the node edge file first for GAE input, where 1 here represents the 1 Mbp resolution. 

Then run GAE for gernerating the inputation matrix by using "python gae.py --nodes 6106 --fin gm12_edge --fout gm12_mat --lr 0.001", where --nodes is the number of nodes, --fin is the edge file, and --fout is the output file for storing the reconstructed matrix. --lr is to set the learning rate. The more details you can see in the readme file in GAE_package file. 

Here, we can convert the GAE imputed matrix to the SCW input we need by using "perl mat_to_hic.pl human_list 1 gm12_mat >gm12_imp.hic"; where 1 is the resolution, human_list is the chromosome length file gm12_mat is the GAE output matrix.

Similar to one-step GAE, we run the following to generate the modeling:

"g++ -o sog -pthread -std=c++11 one_step_GAE.cpp";

"./sog -i ../../example_data/gm12.hic -res 1 -len ../../example_data/human_list -np 40 -imp ../../GAE_package/gm12_imp.hic -o out_g12_1mb"; where -imp is followed by the GAE-imputed Hi-C file, the rest are the same with one-step-2D-Gaussian.

3. Two step 2D Gaussian is the same as one step version.

"g++ -o st2d -pthread -std=c++11 two_step_2D_Gaussian.cpp";

"./st2d -i ../../example_data/gm12.hic -res 1 -len ../../example_data/human_list -np 40 -o out_g12_1mb";

4. Two step 2D Gaussian high resolution is also the same.

"g++ -o st2dh -pthread -std=c++11 scw_two_step_high_resolution.cpp";

"./st2dh -i ../../example_data/gm12.hic -res 0.05 -len ../../example_data/human_list -np 40 -o out_g12_50kb";  

5. Two step GAE is similar to one step GAE version, but needs two GAE-imputed Hi-C file, one at 10 Mb resolution, the other is the target resolution.

"./stg -i ../../example_data/gm12.hic -res 1 -len ../../example_data/human_list -np 40 -imp0 ../../GAE_package/gm12_imp_10mb.hic -imp ../../GAE_package/gm12_imp_1mb.hic -o out_g12_1mb";# -imp0 is followed by the 10 mb GAE imputated Hi-C, -imp is followed by the target that same with -res followed resolution.

6. Visualization. The output from our tool is 3D coordinates like "1497.45 1502.96 1499.08 1_1 0-10000000", the first three columns are X Y Z coordinates, 1_1 is the chromosome, last column is the beads position. In the output file, the last 46 lines (this number depends on how many chromomoes in your hic file, could be 20, 23, 40, or 46)are the chromosome summary with the bead numbers at the target resolution, which shown as "1_1 230, 1_0 228, 2_1 242, 2_0 242 ...". We provide the tool for coverting the orignal 3D coords to pdb file. Before using it, you need to ignore the last chromosome order lines. We can simply use the following "head -n -46 out_g12_1mb > out_g12_temp"; then go to the example_data folder, using "perl convert_single.pl perl convert_single.pl ../two-step-SCW/two-step-2D-Gaussian/out_g12_temp g12.pdb"; to generate the pdb file for Pymol (not included, you should download it yourself) to visualize.
