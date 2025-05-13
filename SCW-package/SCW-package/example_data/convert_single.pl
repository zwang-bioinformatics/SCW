my $file = $ARGV[0];
my $file1 = $ARGV[1];
my $biggest = 100;
my $convert = $ARGV[3];
open(I,$file);
my $big=0.0;
my @a;
my $l=0;
#print $convet."\n";

my $temp="temp.txt";
open(OUT1,">$temp");
while(<I>){
	my $line = $_;
	my @item = split(/\s+/,$line);
	push @a,[@item];
	for(my $i = 0; $i<4;$i++){
		my $temp = $item[$i];
		if($big < $temp){
			$big = $temp;
		}
	}
	$l++;
}
close(I);
my $ratio=0.0;
$ratio = int($big/$biggest);
#print $ratio."\n";
if($ratio<1 ){
	$ratio=($big/$biggest);
}
if($convert eq "no"){
	$ratio=1;#ratio is 1, dont convert, use the original coordinates 
}
#print $ratio."\n";

my @a1;
for(my $i=0;$i<$l;$i++){
	my @a2;
	$a2[0]=($a[$i][0]/$ratio);
	$a2[1]=($a[$i][1]/$ratio);
	$a2[2]=($a[$i][2]/$ratio);
	push @a1,[@a2];
	print OUT1 $a2[0]." ".$a2[1]." ".$a2[2]."\n";
	#print $a2[0]." ".$a2[1]." ".$a2[2]."\n";	
}
`perl coords_2_pdb_linear.pl temp.txt $file1`;
`rm temp.txt`;
close(I);
close(OUT1);
