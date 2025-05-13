my $file = $ARGV[0];
my $res = $ARGV[1];
$res = $res * 1000000;
open(IN,$file);
my $l = 0;
my %hash;
while(<IN>){
	my $line = $_;
	chomp $line;
	my @item = split(/\s+/,$line);
	my $len = int($item[1]/$res) + 1;
	#print $item[0]." $len\n";
	for(my $i = 0; $i < $len; $i++){
		$hash{$item[0]." ".$i} = $l;
		$l++;
	}
}
#print $l."\n";
close(IN);
open(IN2,$ARGV[2]);
my %hash2;
while(<IN2>){
	my $line = $_;
	chomp $line;
	my @item = split(/\s+/,$line);
	my $s = int($item[1]/$res);
	my $e = int($item[3]/$res);
	my $p1 = $hash{$item[0]." ".$s};
	my $p2 = $hash{$item[2]." ".$e};
	if($p1 ne $p2){
		$hash2{$p1." ".$p2} += 1;
		$hash2{$p2." ".$p1} += 1;
	}else{
		$hash2{$p1." ".$p2} += 1;
	}
}
close(IN2);
foreach my $key (keys %hash2){
	print "$key $hash2{$key}\n";
}
