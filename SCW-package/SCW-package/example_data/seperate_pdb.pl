my $file = $ARGV[0];
open(IN,$file);
my %hash;
my $l = 0;
while(<IN>){
	my $line = $_;
	chomp $line;
	#my @item = split(/\s+/,$line);
	$hash{$l} = $line;
	$l++;	
}
close(IN);
open(IN2,$ARGV[1]);
my $dir = $ARGV[2];
`mkdir $dir`;
my $p = 0;
while(<IN2>){
	my $line = $_;
	chomp $line;
	my @item = split(/\s+/,$line);
	open(OUT,">$dir/$item[0].pdb");
	#print OUT $hash
	#close(OUT);
	my $m = 0;
	my $n = 0;
	$m = $p;
	#print "$m ";	
	$p = $p + $item[1];
	$n = $p;
	#print $m." $n\n";
	for(my $i = $m; $i < $n; $i++){
		print OUT $hash{$i}."\n";			
	}
	close(OUT);		
}
close(IN2);
