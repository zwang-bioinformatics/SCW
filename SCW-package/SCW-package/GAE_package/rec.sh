perl read.pl human_list4 ../dipc/gm02.hic >gm2_edg
perl read.pl human_list4 ../dipc/gm03.hic >gm3_edg
perl read.pl human_list4 ../dipc/gm05.hic >gm4_edg
perl read.pl human_list4 ../dipc/gm06.hic >gm6_edg
perl read.pl human_list4 ../dipc/gm07.hic >gm7_edg
perl read.pl human_list4 ../dipc/gm09.hic >gm9_edg
perl read.pl human_list4 ../dipc/gm11.hic >gm11_edg
perl read.pl human_list4 ../dipc/gm12.hic >gm12_edg
perl read.pl human_list4 ../dipc/gm13.hic >gm13_edg
perl read.pl human_list4 ../dipc/gm14.hic >gm14_edg
perl read.pl human_list4 ../dipc/gm15.hic >gm15_edg
perl read.pl human_list4 ../dipc/gm17.hic >gm17_edg
perl read.pl ../human_list5 ../dipc/pbmc01.hic >p1_edg
perl read.pl ../human_list5 ../dipc/pbmc02.hic >p2_edg
perl read.pl ../human_list5 ../dipc/pbmc03.hic >p3_edg
perl read.pl ../human_list5 ../dipc/pbmc04.hic >p4_edg
perl read.pl ../human_list5 ../dipc/pbmc05.hic >p5_edg
perl read.pl ../human_list5 ../dipc/pbmc06.hic >p6_edg
perl read.pl ../human_list5 ../dipc/pbmc07.hic >p7_edg
perl read.pl ../human_list5 ../dipc/pbmc08.hic >p8_edg
perl read.pl ../human_list5 ../dipc/pbmc09.hic >p9_edg
perl read.pl ../human_list5 ../dipc/pbmc10.hic >p10_edg
perl read.pl ../human_list5 ../dipc/pbmc11.hic >p11_edg
perl read.pl ../human_list5 ../dipc/pbmc12.hic >p12_edg
perl read.pl ../human_list5 ../dipc/pbmc13.hic >p13_edg
perl read.pl ../human_list5 ../dipc/pbmc14.hic >p14_edg
perl read.pl ../human_list5 ../dipc/pbmc15.hic >p15_edg
perl read.pl ../human_list5 ../dipc/pbmc16.hic >p16_edg
perl read.pl ../human_list5 ../dipc/pbmc17.hic >p17_edg
perl read.pl ../human_list5 ../dipc/pbmc18.hic >p18_edg




python gae.py --nodes 6106 --fin gm3_edg --fout g3_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm5_edg --fout g5_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm6_edg --fout g6_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm7_edg --fout g7_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm9_edg --fout g9_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm11_edg --fout g11_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm12_edg --fout g12_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm13_edg --fout g13_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm14_edg --fout g14_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm15_edg --fout g15_out --lr 0.001 &
python gae.py --nodes 6106 --fin gm17_edg --fout g17_out --lr 0.001 &

python gae.py --nodes 5950 --fin p1_edg --fout p1_out --lr 0.001 &
python gae.py --nodes 5950 --fin p2_edg --fout p2_out --lr 0.001 &
python gae.py --nodes 5950 --fin p3_edg --fout p3_out --lr 0.001 &
python gae.py --nodes 5950 --fin p4_edg --fout p4_out --lr 0.001 &
python gae.py --nodes 5950 --fin p5_edg --fout p5_out --lr 0.001 &
python gae.py --nodes 5950 --fin p6_edg --fout p6_out --lr 0.001 &
python gae.py --nodes 5950 --fin p7_edg --fout p7_out --lr 0.001 &
python gae.py --nodes 5950 --fin p8_edg --fout p8_out --lr 0.001 &
python gae.py --nodes 5950 --fin p9_edg --fout p9_out --lr 0.001 &
python gae.py --nodes 5950 --fin p10_edg --fout p10_out --lr 0.001 &
python gae.py --nodes 5950 --fin p11_edg --fout p11_out --lr 0.001 &
python gae.py --nodes 5950 --fin p12_edg --fout p12_out --lr 0.001 &
python gae.py --nodes 5950 --fin p13_edg --fout p13_out --lr 0.001 &
python gae.py --nodes 5950 --fin p14_edg --fout p14_out --lr 0.001 &
python gae.py --nodes 5950 --fin p15_edg --fout p15_out --lr 0.001 &
python gae.py --nodes 5950 --fin p16_edg --fout p16_out --lr 0.001 &
python gae.py --nodes 5950 --fin p17_edg --fout p17_out --lr 0.001 &
python gae.py --nodes 5950 --fin p18_edg --fout p18_out --lr 0.001 &
