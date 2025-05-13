#!/usr/bin/env perl
use warnings;
use strict;

my $f_coordinates = $ARGV[0];
my $f_output = $ARGV[1];

my @coords;
open IN, "$f_coordinates" or die $!;
while (my $line = <IN>) {
	chomp $line; 
  	$line =~ s/^\s+//; 
  	$line =~ s/\s+$//;
  	my @items = split(/\s+/, $line);
  	push(@coords, \@items);
}
close IN;

my $num_points = @coords;

my $pdb_id = 0;
open OUT, ">$f_output";
for (my $i = 0; $i < $num_points; $i++) {
  	if ($i == 0) {
    		$pdb_id++;
    		print_pdb_atom($coords[0][0], $coords[0][1], $coords[0][2], $pdb_id);
  	} else {
    		$pdb_id++;
    		print_pdb_atom($coords[$i][0], $coords[$i][1], $coords[$i][2], $pdb_id);
  	}
}
close OUT;

sub print_pdb_atom {
  	my $x = $_[0];
	my $y = $_[1];
  	my $z = $_[2];
  	my $curr_index = $_[3];

  	$x = sprintf "%0.3f", $x;
  	$y = sprintf "%0.3f", $y;
  	$z = sprintf "%0.3f", $z;
	#$x = $x."00";
	
  	printf OUT "ATOM";
  	printf OUT "%7s", $curr_index;
  	printf OUT "  ";
  	printf OUT "%-3s", "CA";
  	printf OUT "%4s", "LIG";
  	printf OUT "%6s", $curr_index;

  	printf OUT "%12s", $x;
  	printf OUT "%8s", $y;
  	printf OUT "%8s", $z;
  	printf OUT "%6s", "1.00";
  	printf OUT "%6s", "75.00";
  	printf OUT "\n";
}

