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
		$hash{$l} = $item[0]." $i";
		$l++;
	}
}
#print $l."\n";
close(IN);
open(IN2,$ARGV[2]);
my $l1 = 0;
while(<IN2>){
	my $line = $_;
	chomp $line;
	my @item = split(/\s+/,$line);
	my $l2 = 0;
	foreach my $it (@item){
		print $hash{$l1}." $hash{$l2} $it\n";
		$l2++;	
	}
	
	$l1++;
}
close(IN2);
