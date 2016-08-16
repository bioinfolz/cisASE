#!/usr/bin/perl
use strict;
use warnings;
my $sample = 50000;  ## $sample reads will read to determine the format

if( !(scalar(@ARGV)) || ($ARGV[0] eq "-h") || ($ARGV[0] =~ m/-help/) ){
	&usage;
	}
	
open (IN,$ARGV[0]) or die "\n can not open file $ARGV[0] \n";
my $counter = 0;
my (%Hash,$baseQ,@bq,$min,$max);
print "\nsampling..\n";
	while(<IN>){
	$_ = <IN>;
	$_ = <IN>;
	$_ = <IN>;

	s/\n//;
	s/\r//;
		if($counter >= $sample){
		last;
		}

	my @a1 = split('',$_);
		for(my $i = 0;$i<scalar(@a1);++$i){
		$baseQ = ord($a1[$i]); 
		$Hash{$baseQ} = 1;
		}
	$counter  += 1;
	}  
close IN;


	for my $key (sort {$a <=> $b}keys %Hash){
	push(@bq,$key);
	}

$min = $bq[1];  
$max = $bq[scalar(@bq)-1]; 
print "Quality range: ",$bq[0],"-",$bq[scalar(@bq)-1],"\n";

print "format: ",format_decision($min,$max),"\n\n";

	sub format_decision{
	my($min,$max)  = @_;
	my $qualFormat;
		if( ($min >= 33) && ($max <= 76) ){
		$qualFormat = "Sanger or Illumina 1.8+ (offset by 33)";
		}
		elsif( ($min >= 67) && ($max <= 106) ){
		$qualFormat = "Illumina 1.5+ (offset by 64)";
		}
		
		elsif( ($min >= 64) && ($max <= 106) ){
		$qualFormat = "Illumina 1.3+ (offset by 64)";
		}
		
		elsif( ($min >= 59) && ($max <= 106) ){
		$qualFormat = "Solexa (offset by 64)";
		}
		else{
		$qualFormat = "Unknown";
		}
	return($qualFormat);
	}  

sub usage {
  die(qq/
Usage:   perl GetFastqQualityEncoding.pl <input fastq> \n
possible formats: Sanger( or Illumina 1.8+) , Solexa, Illumina 1.3+ and Illumina 1.5+
\n/);
}

