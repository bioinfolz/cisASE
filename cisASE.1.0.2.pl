#!/usr/bin/perl
####report time for each step:
 # ***************************************************************************
 #
 #                                   cisASE
 #      A likelihood-based method for detecting putative cis-regulated allele
 #		            specific expression in RNA sequencing data 
 #
 # ***************************************************************************
 #
 # Copyright Shanghai Institutes for biological Sciences,CAS. 2015 
 # Author :  Zhi Liu
 # Correspondings: Xiao Dong; Yixue Li 
 # liuzhi2011@sibs.ac.cn                 
 # biosinodx@gmail.com                 
 # 
 # This software is a computer program whose purpose is to detect allele specific 
 # expression events from RNA-seq data.                   

use Getopt::Std;
use strict;
&main;
exit;

sub main {
  &usage if (@ARGV < 1);
  my $command = shift(@ARGV);
  my %func = (SNV=>\&SNV,EXON=>\&EXON,GENE=>\&GENE);
  die("Unknown command \"$command\".\n") if (!defined($func{$command}));
  &{$func{$command}};
}

sub SNV(){
my ($SNPlist,$RNApileup,$DNApileup,$output,$SIMU,$Dphred,$Rphred,$DNAbias,$Method,$dD,$dR,$dA,$rD,$rR,$rA);
my $usage = "Usage: perl $0 SNV \n\n\t"
			."-L <List of SNVs>\n\t"
			."-X <RNA pileup file >\n\t"
			."-Y <DNA pileup file>\n\t"
			."-o <output file >\n\t"
			."-S <Times for simulation; Default:2000>\n\t"
			."-B <DNA bias; Default: simulated from RNA sequencing data>\n\t"
			."-F <DNA fastq encoding:64= Solexa,Illumina 1.3+,Illumina 1.5+; 33= Sanger,Illumina 1.8+ >\n\t"
			."-f <RNA fastq encoding: 64= Solexa,Illumina 1.3+,Illumina 1.5+;  33= Sanger,Illumina 1.8+ >\n\t"		
			."-M <1:chisquare test; 2: Likelihood model; 3:Both; Default:2>\n\t"
			."-D <minimum DNA read depth Default:20>\n\t"
			."-R <minimum DNA Reference Allele counts Default:5>\n\t"
			."-A <minimum DNA Alternative Allele counts Default:5>\n\t"
			."-d <minimum RNA read depth Default:20 >\n\t"
			."-r <minimum RNA Reference Allele counts Default:0>\n\t"
			."-a <minimum RNA Alternative Allele counts Default:0>\n\t"
			."\n";
my %opts=();

getopts('L:X:Y:o:S:B:F:f:M:D:R:A:d:r:a',\%opts);
if(defined $opts{L}){  $SNPlist = $opts{L};} else {die "$! \n$usage";}
if(defined $opts{X}){  $RNApileup = $opts{X};} else {die "$! \n$usage";}
if(defined $opts{Y}){  $DNApileup = $opts{Y};}
if(defined $opts{o}){  $output = $opts{o};} else {die "$! \n$usage";}
if(defined $opts{S}){  $SIMU = $opts{S};}else{$SIMU=2000;}
if(defined $opts{B}){  $DNAbias = $opts{B};}
if(defined $opts{F}){  $Dphred = $opts{F};}
if(defined $opts{f}){  $Rphred = $opts{f};} else {die "$! \n$usage";}
if(defined $opts{M}){  $Method = $opts{M};} else {$Method="2";}
if(defined $opts{D}){  $dD = $opts{D};} else {$dD=20;}
if(defined $opts{R}){  $dR = $opts{R};} else {$dR=5;}
if(defined $opts{A}){  $dA = $opts{A};} else {$dA=5;}
if(defined $opts{d}){  $rD = $opts{d};} else {$rD=20;}
if(defined $opts{r}){  $rR = $opts{r};} else {$rR=0;}
if(defined $opts{a}){  $rA = $opts{a};} else {$rA=0;}

my %hashs=();
my %hashsite=();
my ($Dall,$Df,@out);
my @oa=("SNP_id","Ref_allele","Alt_allele","DNA_cnt(Ref:Alt:Other)");
my @ob=("DNAref_freq","RNAref_freq","Adjusted_RNA_freq");
if(!defined $opts{Y}){pop @oa; shift @ob;}
if($Method ==1){
	push(@out,join("\t",@oa,"RNA_cnt(Ref:Alt:Other)",@ob,"Chisq-p")."\n")
}elsif($Method==2){
	push(@out,join("\t",@oa,"RNA_cnt(Ref:Alt:Other)",@ob,"LLR")."\n")
}elsif($Method==3){
	push(@out,join("\t",@oa,"RNA_cnt(Ref:Alt:Other)",@ob,"Chisq-p","LLR")."\n")
}else{die "Wrong Method parameter!";}


my (%SNP,%SNP2,%SNP3,$RNAbias,$DNAfreq);
my %SNP=&ListRead($SNPlist);
($RNAbias,%SNP2)=&RNARead($RNApileup,$rR,$rA,$rD,\%SNP);
if(defined $opts{Y}){
	($DNAfreq,%SNP3)=&DNARead($DNApileup,$dR,$dA,$dD,\%SNP2,\%SNP);
	$DNAbias="Cal";
}else{%SNP3=%SNP2;
	if(!defined $opts{B}){$DNAbias=$RNAbias;}
	$DNAfreq=$DNAbias;
}


my (@srall,@sdall,@srbq);
foreach my $key( sort keys %SNP3){
	my ($tdref,$tdalt,$rmut,$dmut,$DNAratio);
	my @sep=split/\t/,$SNP3{$key};
	my $REF=$sep[0];my $ALT=$sep[1];
	my $SNPid=$key;
	my @R_allele=split//,$sep[6];
	my @Rbq_v=split//,$sep[7];
	my $rnaalt=$sep[3];my $rnaall=$sep[5];
	my $rmut=$rnaalt/$rnaall;
	my ($trref,$tralt)=&theory_counts($Rphred,$REF,$ALT,\@R_allele,\@Rbq_v);
	if($#sep> 7){ ##have information at DNA
		$dmut=$sep[9]/$sep[11];
		my @D_allele=split//,$sep[12];
		my @Dbq_v=split//,$sep[13];
		($tdref,$tdalt)=&theory_counts($Dphred,$REF,$ALT,\@D_allele,\@Dbq_v);
	}else{
		$dmut=1-$DNAfreq;
		$tdref=($trref+$tralt)*$DNAfreq; $tdalt=($trref+$tralt)*$dmut;
	}
	$DNAratio=$tdref/($tdalt+$tdref);
	my ($Adj_ratio,$RNAratio)=&Adj($DNAratio,$trref,$tralt);
	my @m=($tdref,$tdalt,$trref,$tralt,$dmut,$rmut,$Rphred);
	if($sep[11] >0){@sdall=(@sdall,$sep[11]);}
	if($sep[5] > 0){@srall=(@srall,$sep[5]);}
	if($sep[7] ne ""){@srbq=(@srbq,$sep[7]);}
	$hashsite{1}=join("\t",$sep[0],$sep[1],$sep[6],$sep[7]);

	
	my $Result=&Method(\%hashsite,\@m,$Method);
	if($DNAbias ne "Cal"){
		push(@out,join("\t",$SNPid,$sep[0],$sep[1])."\t".join(":",$sep[2],$sep[3],$sep[4])."\t",join("\t",$RNAratio,$Adj_ratio,$Result)."\n");
	}
	else{
		push(@out,join("\t",$SNPid,$sep[0],$sep[1])."\t".join(":",$sep[8],$sep[9],$sep[10])."\t",join(":",$sep[2],$sep[3],$sep[4])."\t",join("\t",$DNAratio,$RNAratio,$Adj_ratio,$Result)."\n");
	}		
}

if(!defined $opts{Y}){@sdall=@srall;}
if($Method ne "1"){
my ($fdr);
	if($SIMU > 0){
		my @couts=(1,1);
		my($fdr1,$fdr2)=&Simulation(\@sdall,\@srall,\@srbq,\@couts,$DNAfreq,$Rphred,$SIMU);
		$fdr="##LLR cutoff for significant level of 0.05:".$fdr1." under ".$SIMU." times simulation; LLR cutoff for significant level of 0.01:".$fdr2." under ".$SIMU." times simulation\n";
	}else{
		$fdr="##Recommanded LLR cutoff for significant level of 0.05 : 0.84; for significant level of 0.01 : 1.44 \n";
	}
unshift(@out,$fdr);
}
open(OUT,">$output");
print OUT @out;
}


sub EXON(){
my ($SNPlist,$RNApileup,$annofile,$Phase,$DNApileup,$annotation,$output,$SIMU,$Dphred,$Rphred,$DNAbias,$Method,$dD,$dR,$dA,$rD,$rR,$rA);
my $usage = "Usage: perl $0 EXON \n\n\t"
			."-L <List of SNVs>\n\t"
			."-X <RNA pileup file >\n\t"
			."-Y <DNA pileup file>\n\t"
			."-o <output file >\n\t"
			."-N <Annotation file>\n\t"
			."-S <Times for simulation; Default:2000>\n\t"
			."-B <DNA bias; Default: simulated from RNA sequencing data>\n\t"
			."-F <DNA fastq encoding:64= Solexa,Illumina 1.3+,Illumina 1.5+; 33= Sanger,Illumina 1.8+ >\n\t"
			."-f <RNA fastq encoding: 64= Solexa,Illumina 1.3+,Illumina 1.5+;  33= Sanger,Illumina 1.8+ >\n\t"		
			."-M <1:chisquare test; 2: Likelihood model; 3:Both; Default:2>\n\t"
			."-D <minimum DNA read depth Default:20>\n\t"
			."-R <minimum DNA Reference Allele counts Default:5>\n\t"
			."-A <minimum DNA Alternative Allele counts Default:5>\n\t"
			."-d <minimum RNA read depth Default:20 >\n\t"
			."-r <minimum RNA Reference Allele counts Default:0>\n\t"
			."-a <minimum RNA Alternative Allele counts Default:0>\n\t"
			."-P <phasing file >\n\t"			
			."\n";
my %opts=();

getopts('L:X:Y:o:N:S:B:F:f:M:P:D:R:A:d:r:a',\%opts);
if(defined $opts{L}){  $SNPlist = $opts{L};} else {die "$! \n$usage";}
if(defined $opts{X}){  $RNApileup = $opts{X};} else {die "$! \n$usage";}
if(defined $opts{Y}){  $DNApileup = $opts{Y};} else {}
if(defined $opts{N}){  $annofile = $opts{N};} else {die "$! \n$usage";}
if(defined $opts{o}){  $output = $opts{o};} else {die "$! \n$usage";}
if(defined $opts{S}){  $SIMU = $opts{S};}else{$SIMU=2000;}
if(defined $opts{B}){  $DNAbias = $opts{B};}
if(defined $opts{F}){  $Dphred = $opts{F};}
if(defined $opts{f}){  $Rphred = $opts{f};} else {die "$! \n$usage";}
if(defined $opts{M}){  $Method = $opts{M};} else {$Method="2";}
if(defined $opts{P}){  $Phase = $opts{P};} else {$Phase="NO";}
if(defined $opts{D}){  $dD = $opts{D};} else {$dD=20;}
if(defined $opts{R}){  $dR = $opts{R};} else {$dR=5;}
if(defined $opts{A}){  $dA = $opts{A};} else {$dA=5;}
if(defined $opts{d}){  $rD = $opts{d};} else {$rD=20;}
if(defined $opts{r}){  $rR = $opts{r};} else {$rR=0;}
if(defined $opts{a}){  $rA = $opts{a};} else {$rA=0;}

my ($Dall,$Df,%hashs,%hashsite,%hash,%hashlab,@out);
my (%SNP,%SNP2,%Anno,%SNP3,@SNP4,%phase);
if(defined $opts{P}){ 
	%phase=&PhaseRead($Phase)
}

open(IN4,"$annofile")or die $!;
while (<IN4>){
	chomp;
	my $sort++;
	my @sep=split/\t/;
	my $chr="chr".$sep[0];
	if($sep[4]<$sep[5]){
		$Anno{$chr}->{$sep[2]}->{start}=$sep[4];
		$Anno{$chr}->{$sep[2]}->{end}=$sep[5];
	}else{
		$Anno{$chr}->{$sep[2]}->{start}=$sep[5];
		$Anno{$chr}->{$sep[2]}->{end}=$sep[4];
	}		
}

my %sort;
foreach my $chr(keys %Anno){
	my @sort = sort { $Anno{$chr}{$a}{start} <=> $Anno{$chr}{$b}{start} } keys %{$Anno{$chr}};
	$sort{$chr}=join(",",@sort);
	}

my (%SNP,%SNP2,%SNP3,$RNAbias,$DNAfreq);
my %SNP=&ListRead($SNPlist);
($RNAbias,%SNP2)=&RNARead($RNApileup,$rR,$rA,$rD,\%SNP);
if(defined $opts{Y}){
	($DNAfreq,%SNP3)=&DNARead($DNApileup,$dR,$dA,$dD,\%SNP2,\%SNP);
	$DNAbias="Cal";
}else{%SNP3=%SNP2;
	if(!defined $opts{B}){$DNAbias=$RNAbias;}
	$DNAfreq=$DNAbias;
}

SNPMark:
foreach my $key(sort keys %SNP3){
	my @sep=split/\:/,$key;
	my @sort=split/\,/,$sort{$sep[0]};
	foreach my $exons(@sort){
		if(exists $Anno{$sep[0]}{$exons}){
			if(($sep[1] >=$Anno{$sep[0]}{$exons}{start})&&($sep[1] <=$Anno{$sep[0]}{$exons}{end})){
				push(@SNP4,join("\t",$key,$exons,$SNP3{$key}));
			}if($sep[1] < $Anno{$sep[0]}{$exons}{start}){
				next SNPMark;
			}
		}
	}
}

if(scalar @SNP4 == 0) { die "No snv is annotated!\n";}

my @oa=("DNAhap1_freq","RNAhap1_freq","Adj_RNAhap1_freq");
my @ob=("heterogeneity p-value");
if($DNAbias ne "Cal"){shift @oa;}
if($Method ==1){
	push(@out,join("\t","EXONid",join("\t",@oa),"Chisq-p",join("\t",@ob))."\n");
}elsif($Method==2){
	push(@out,join("\t","EXONid",join("\t",@oa),"LLR",join("\t",@ob))."\n");
}elsif($Method==3){
	push(@out,join("\t","EXONid",join("\t",@oa),"Chisq-p","LLR",join("\t",@ob))."\n");
}else{die "Wrong Method parameter!";}

my (@srall,@sdall,@srbq);
foreach my $key( @SNP4){
	my ($tdref,$tdalt,$rmut,$dmut,$DNAratio);	
	my @counts=split/\t/,$key;
	my $SNPid=shift (@counts);
	my $EXONid=shift(@counts);
	if(defined $opts{P}){
		if($phase{$SNPid}{ref} eq $counts[1]){
			@counts=&Change(\@counts);				
		}
	}else{
		if($counts[2] < $counts[3]){
		@counts=&Change(\@counts);	
		}
	}
	my @R_allele=split//,$counts[6];
	my @Rbq_v=split//,$counts[7];
	my ($trref,$tralt)=&theory_counts($Rphred,$counts[0],$counts[1],\@R_allele,\@Rbq_v);
	if($#counts> 7){ ##have information at DNA
		my @D_allele=split//,$counts[12];
		my @Dbq_v=split//,$counts[13];
		($tdref,$tdalt)=&theory_counts($Dphred,$counts[0],$counts[1],\@D_allele,\@Dbq_v);
	}else{
		$tdref=($trref+$tralt)*$DNAfreq; $tdalt=($trref+$tralt)-$tdref;
		}
	$DNAratio=$tdref/($tdref+$tdalt);
	my ($tAdj_ratio,$RNAratio)=&Adj($DNAratio,$trref,$tralt);
	my $adj_rfirst=($trref+$tralt)*$tAdj_ratio;
	my $adj_rsecond=($trref+$tralt)-$adj_rfirst;
	$hashlab{$EXONid}->{fisherref}=$hashlab{$EXONid}->{fisherref}."\t".$adj_rfirst;
	$hashlab{$EXONid}->{fisheralt}=$hashlab{$EXONid}->{fisheralt}."\t".$adj_rsecond;
	$hashlab{$EXONid}->{counts} +=1;
	my $i=$hashlab{$EXONid}->{counts};
	if($counts[11] >0){@sdall=(@sdall,$counts[11]);}
	if($counts[5] > 0){@srall=(@srall,$counts[5]);}
	if($counts[7] ne ""){@srbq=(@srbq,$counts[7]);}
	$hash{$EXONid}{$i}=join("\t",$counts[0],$counts[1],$counts[6],$counts[7]);
	$hashlab{$EXONid}->{rf} += $counts[2];
	$hashlab{$EXONid}->{rs} += $counts[3];
	$hashlab{$EXONid}->{Rall}+= $counts[5];
	$hashlab{$EXONid}->{trf} += $trref;
	$hashlab{$EXONid}->{trs} += $tralt;
	$hashlab{$EXONid}->{df} += $counts[8];
	$hashlab{$EXONid}->{ds} += $counts[9];
	$hashlab{$EXONid}->{Dall} += $counts[11];
	$hashlab{$EXONid}->{tdf} += $tdref;
	$hashlab{$EXONid}->{tds} += $tdalt;
}


foreach my $exon (sort keys %hash){
	my ($dref,$dalt,$dmut,$rmut,$DNAratio)=(0,0,0,0,0);
	my $rref=int(0.5+$hashlab{$exon}{trf}); 
	my $ralt=int(0.5+$hashlab{$exon}{trs});
	my $dref=int(0.5+$hashlab{$exon}{tdf});
	my $dalt=int(0.5+$hashlab{$exon}{tds});
	my %hashexon=%{$hash{$exon}};
	if($hashlab{$exon}{Dall} == 0){
		$dmut=1-$DNAfreq;
	}else{
		$dmut=$hashlab{$exon}{ds}/$hashlab{$exon}{Dall};
	}
	$DNAratio=$hashlab{$exon}{tdf}/($hashlab{$exon}{tds}+$hashlab{$exon}{tdf});
	my $trref=$hashlab{$exon}{trf};
	my $tralt=$hashlab{$exon}{trs};
	my $tdref=$hashlab{$exon}{tdf};
	my $tdalt=$hashlab{$exon}{tds};
	my ($Adj_ratio,$RNAratio)=&Adj($DNAratio,$trref,$tralt);	
	$rmut=$hashlab{$exon}{rs}/$hashlab{$exon}{Rall};
	my $fref=$hashlab{$exon}{fisherref};
	my @fisherref=split/\t/,$fref;
	my $falt=$hashlab{$exon}{fisheralt};
	my @fisheralt=split/\t/,$falt;
	my ($stat,$pfisher);
	if($hashlab{$exon}{counts}==1){ $pfisher = "NA";}
	else{
		shift  @fisherref;
		shift  @fisheralt;
		$pfisher=&heteroTest(\@fisherref,\@fisheralt);
	}
	my @m=($dref,$dalt,$rref,$ralt,$dmut,$rmut,$Rphred);
	my $Res=&Method(\%hashexon,\@m,$Method);
	if($DNAbias ne "Cal"){
		push(@out,join("\t",$exon,$RNAratio,$Adj_ratio,$Res,$pfisher)."\n");	
	}else{
		push(@out,join("\t",$exon,$DNAratio,$RNAratio,$Adj_ratio,$Res,$pfisher)."\n");
	}
}

if(!defined $opts{Y}){@sdall=@srall;}
if($Method ne 1){
	my ($fdr);
	if($SIMU > 0){	
		my($Drefractio,%simu,%hashloci,@cnts);
		my $i=0;
			foreach my $exon(keys %hashlab){
				push(@cnts,$hashlab{$exon}{counts});
			}
		$Drefractio=$DNAfreq;
		my ($fdr1,$fdr2)=&Simulation(\@sdall,\@srall,\@srbq,\@cnts,$Drefractio,$Rphred,$SIMU);
		$fdr="##LLR cutoff for significant level of 0.05:".$fdr1." under ".$SIMU." times simulation; LLR cutoff for significant level of 0.01:".$fdr2." under ".$SIMU." times simulation\n";
	}else{
		$fdr="##Recommanded LLR cutoff for significant level of 0.05 : 0.82; for significant level of 0.01 : 1.48 \n";
	}
unshift(@out,$fdr);
}

open(OUT, ">$output");
print OUT @out;	
}

sub GENE(){
my ($SNPlist,$RNApileup,$Cons,$annofile,$Phase,$DNApileup,$annotation,$output,$SIMU,$Dphred,$Rphred,$DNAbias,$Method,$dD,$dR,$dA,$rD,$rR,$rA);
my $usage = "Usage: perl $0 GENE \n\n\t"
			."-L <List of SNVs>\n\t"
			."-X <RNA pileup file >\n\t"
			."-Y <DNA pileup file>\n\t"
			."-o <output file >\n\t"
			."-N <Annotation file>\n\t"
			."-C <Consider exon statue (1) or not(0); Default:0>\n\t"
			."-S <Times for simulation; Default:2000>\n\t"
			."-B <DNA bias; Default: simulated from RNA sequencing data>\n\t"
			."-F <DNA fastq encoding: 64= Solexa,Illumina 1.3+,Illumina 1.5+; 33= Sanger,Illumina 1.8+ >\n\t"
			."-f <RNA fastq encoding: 64= Solexa,Illumina 1.3+,Illumina 1.5+; 33= Sanger,Illumina 1.8+ >\n\t"		
			."-M <1:chisquare test; 2: Likelihood model; 3:Both Default:2>\n\t"
			."-D <minimum DNA read depth Default:20>\n\t"
			."-R <minimum DNA Reference Allele counts Default:5>\n\t"
			."-A <minimum DNA Alternative Allele counts Default:5>\n\t"
			."-d <minimum RNA read depth Default:20 >\n\t"
			."-r <minimum RNA Reference Allele counts Default:0>\n\t"
			."-a <minimum RNA Alternative Allele counts Default:0>\n\t"
			."-P <phasing file >\n\t"			
			."\n";
my %opts=();

getopts('L:X:Y:o:N:C:S:B:F:f:M:P:D:R:A:d:r:a',\%opts);
if(defined $opts{L}){  $SNPlist = $opts{L};} else {die "$! \n$usage";}
if(defined $opts{X}){  $RNApileup = $opts{X};} else {die "$! \n$usage";}
if(defined $opts{Y}){  $DNApileup = $opts{Y};} else {}
if(defined $opts{N}){  $annofile = $opts{N};} else {die "$! \n$usage";}
if(defined $opts{C}){  $Cons = $opts{C};}else{$Cons=0;}
if(defined $opts{o}){  $output = $opts{o};} else {die "$! \n$usage";}
if(defined $opts{S}){  $SIMU = $opts{S};}else{$SIMU=2000;}
if(defined $opts{B}){  $DNAbias = $opts{B};}
if(defined $opts{F}){  $Dphred = $opts{F};}
if(defined $opts{f}){  $Rphred = $opts{f};} else {die "$! \n$usage";}
if(defined $opts{M}){  $Method = $opts{M};} else {$Method="2";}
if(defined $opts{P}){  $Phase = $opts{P};} else {$Phase="NO";}
if(defined $opts{D}){  $dD = $opts{D};} else {$dD=20;}
if(defined $opts{R}){  $dR = $opts{R};} else {$dR=5;}
if(defined $opts{A}){  $dA = $opts{A};} else {$dA=5;}
if(defined $opts{d}){  $rD = $opts{d};} else {$rD=20;}
if(defined $opts{r}){  $rR = $opts{r};} else {$rR=0;}
if(defined $opts{a}){  $rA = $opts{a};} else {$rA=0;}

my (@out,%hashC,%hashClab,%hashA,%hashAlab,%hash_alternative,%hash_constitutive);
my (%SNP,%SNP2,%SNP3,%SNP4,%Anno,%phase);
if(defined $opts{P}){ 
	%phase=&PhaseRead($Phase);
}else{}

open(IN4,"$annofile")or die $!;
while (<IN4>){
	chomp;
	my @sep=split/\t/;
	my $chr="chr".$sep[0];
	$Anno{$chr}->{$sep[2]}->{C}=$sep[3];
	$Anno{$chr}->{$sep[2]}->{gene}=$sep[1];	
	if($sep[4]<$sep[5]){
		$Anno{$chr}->{$sep[2]}->{start}=$sep[4];
		$Anno{$chr}->{$sep[2]}->{end}=$sep[5];		
	}else{
		$Anno{$chr}->{$sep[2]}->{start}=$sep[5];
		$Anno{$chr}->{$sep[2]}->{end}=$sep[4];
	}		
}

my %sort;
foreach my $chr(keys %Anno){
	my @keys = sort { $Anno{$chr}{$a}{start} <=> $Anno{$chr}{$b}{start} } keys %{$Anno{$chr}};
	$sort{$chr}=join(",",@keys);
}

my (%SNP,%SNP2,%SNP3,$DNAfreq,$RNAbias);
my %SNP=&ListRead($SNPlist);
($RNAbias,%SNP2)=&RNARead($RNApileup,$rR,$rA,$rD,\%SNP);
if(defined $opts{Y}){
	($DNAfreq,%SNP3)=&DNARead($DNApileup,$dR,$dA,$dD,\%SNP2,\%SNP);
	$DNAbias="Cal";
}else{%SNP3=%SNP2;
	if(!defined $opts{B}){$DNAbias=$RNAbias;}
	$DNAfreq=$DNAbias;
}

my @snp4;
SNPMark:
foreach my $key(keys %SNP3){
	my @sep=split/\:/,$key;
	my @sort=split/\,/,$sort{$sep[0]};
	foreach my $exons(@sort){
		if(exists $Anno{$sep[0]}{$exons}){
			if(($sep[1] >=$Anno{$sep[0]}{$exons}{start})&&($sep[1] <=$Anno{$sep[0]}{$exons}{end})){
					if($Cons==1){
						my $con=$Anno{$sep[0]}{$exons}{C};
						push(@snp4,join("\t",$key,$exons,$Anno{$sep[0]}{$exons}{gene},$con,$SNP3{$key}));
						next SNPMark;
					}else{
						push(@snp4,join("\t",$key,$exons,$Anno{$sep[0]}{$exons}{gene},1,$SNP3{$key}));
						next SNPMark;
					}

			}if($sep[1] < $Anno{$sep[0]}{$exons}{start}){
				next SNPMark;
			}
		}
	}	
}

if(scalar @snp4 == 0) { die "No snv is annotated!\n";}
my @oa=("DNAhap1_freq","RNAhap1_freq","Adj_RNAhap1_freq");
my @ob=("heterogeneity p-value","Gene:exon heterogeneity");
if($Cons==0){pop @ob;}
if($DNAbias ne "Cal"){ shift @oa;}
if($Method ==1){
	push(@out,join("\t","Geneid",join("\t",@oa),"Chisq-p",join("\t",@ob))."\n");
}elsif($Method==2){
	push(@out,join("\t","Geneid",join("\t",@oa),"LLR",join("\t",@ob))."\n");
}elsif($Method==3){
	push(@out,join("\t","Geneid",join("\t",@oa),"Chisq-p","LLR",join("\t",@ob))."\n");
}else{die "Wrong Method parameter!";}


my (@srall,@sdall,@srbq);
foreach my $key(@snp4){
	my ($tdref,$tdalt,$rmut,$dmut,$DNAratio);	
	my @counts=split/\t/,$key;
	my $SNPid=shift(@counts);
	my $EXONid=shift(@counts);
	my $GENEid=shift(@counts);
	my $statue=shift(@counts);	
	if(defined $opts{P}){
		if($phase{$SNPid}{ref} eq $counts[1]){
			@counts=&Change(\@counts);				
		}
	}else{
		if($counts[2] < $counts[3]){
			@counts=&Change(\@counts);				
		}
	}
	my @R_allele=split//,$counts[6];
	my @Rbq_v=split//,$counts[7];
	my ($trref,$tralt)=&theory_counts($Rphred,$counts[0],$counts[1],\@R_allele,\@Rbq_v);
	if($#counts> 7){ ##have information at DNA
		my @D_allele=split//,$counts[12];
		my @Dbq_v=split//,$counts[13];
		($tdref,$tdalt)=&theory_counts($Dphred,$counts[0],$counts[1],\@D_allele,\@Dbq_v);
	}else{
		$tdref=($trref+$tralt)*$DNAfreq; $tdalt=($trref+$tralt)-$tdref;
	}	
	if($statue == "1"){		
		$hash_constitutive{$GENEid}=1;
		$hashClab{$GENEid}->{counts}+=1;
		my $i=$hashClab{$GENEid}->{counts};	
		if($counts[11] >0){@sdall=(@sdall,$counts[11]);}
		if($counts[5] > 0){@srall=(@srall,$counts[5]);}
		if($counts[7] ne ""){@srbq=(@srbq,$counts[7]);}		
		$hashC{$GENEid}{$i}=join("\t",$counts[0],$counts[1],$counts[6],$counts[7]);
		$hashClab{$GENEid}->{rf} += $counts[2] ;
		$hashClab{$GENEid}->{rs} += $counts[3];
		$hashClab{$GENEid}->{Rall} += $counts[5];
		$hashClab{$GENEid}->{df} +=$counts[8];
		$hashClab{$GENEid}->{ds} +=$counts[9];
		$hashClab{$GENEid}->{Dall} +=$counts[11];
		$hashClab{$GENEid}->{trf} += $trref;
		$hashClab{$GENEid}->{trs} += $tralt;
		$hashClab{$GENEid}->{tdf} += $tdref;;
		$hashClab{$GENEid}->{tds} += $tdalt;
		$DNAratio=$tdref/($tdref+$tdalt);
	my ($tAdj_ratio,$RNAratio)=&Adj($DNAratio,$trref,$tralt);
	my $adj_rfirst=($trref+$tralt)*$tAdj_ratio;
	my $adj_rsecond=($trref+$tralt)-$adj_rfirst;
	$hashClab{$GENEid}->{fisherref}=$hashClab{$GENEid}->{fisherref}."\t".$adj_rfirst;
	$hashClab{$GENEid}->{fisheralt}=$hashClab{$GENEid}->{fisheralt}."\t".$adj_rsecond;	
	}
	elsif($statue == "0"){
		$hash_alternative{$GENEid}=1;		
		$hashAlab{$GENEid}{$EXONid}->{trf} += $trref;
		$hashAlab{$GENEid}{$EXONid}->{trs} += $tralt;
		$hashAlab{$GENEid}{$EXONid}->{tdf} += $tdref;;
		$hashAlab{$GENEid}{$EXONid}->{tds} += $tdalt;
	}
}

foreach my $gene (sort keys %hash_constitutive){
	my ($dref,$dalt,$dmut,$rmut,$DNAratio,$Adj_ratio,$RNAratio)=(0,0,0,0,0,0,0);
	my $rref=int(0.5+$hashClab{$gene}->{trf}); 
	my $ralt=int(0.5+$hashClab{$gene}->{trs});
	my $dref=int(0.5+$hashClab{$gene}->{tdf}); 
	my $dalt=int(0.5+$hashClab{$gene}->{tds});
	my $counts=$hashClab{$gene}{counts};
	my %hashgene=%{$hashC{$gene}};
	if($hashClab{$gene}{Dall} == 0){
		$dmut=1-$DNAfreq;
	}else{
		$dmut=$hashClab{$gene}{ds}/$hashClab{$gene}{Dall};
	}
	$DNAratio=$hashClab{$gene}{tdf}/($hashClab{$gene}{tdf}+$hashClab{$gene}{tds});
	my $trref=$hashClab{$gene}{trf};
	my $tralt=$hashClab{$gene}{trs};
	($Adj_ratio,$RNAratio)=&Adj($DNAratio,$trref,$tralt);	
	$rmut=$hashClab{$gene}{rs}/$hashClab{$gene}{Rall};
	my $fref=$hashClab{$gene}{fisherref};
	my @fisherref=split/\t/,$fref;
	my $falt=$hashClab{$gene}{fisheralt};
	my @fisheralt=split/\t/,$falt;
	my ($stat,$pfisher);
	use Statistics::Descriptive;
	$stat = Statistics::Descriptive::Full->new();
	if($hashClab{$gene}{counts}==1){ $pfisher = "NA";}
	else{
	    shift  @fisherref;
		shift  @fisheralt;
		$pfisher=&heteroTest(\@fisherref,\@fisheralt);
}
	my @m=($dref,$dalt,$rref,$ralt,$dmut,$rmut,$Rphred);
	my $Res=&Method(\%hashgene,\@m,$Method);
	if ($DNAbias ne "Cal"){
		push(@out,join("\t",$gene,$RNAratio,$Adj_ratio,$Res,$pfisher));
	}else{
		push(@out,join("\t",$gene,$DNAratio,$RNAratio,$Adj_ratio,$Res,$pfisher));
	}
	if($Cons==0){
		push(@out,"\n");
		}
	else{
		if(exists($hash_alternative{$gene})){
			my ($Exon_gene,@pca);
			foreach my $exon ( keys %{$hashAlab{$gene}}){
				my $rall=$hashClab{$gene}{trf}+$hashClab{$gene}{trs};
				my $generef=$Adj_ratio*$rall;
				my $genealt=(1-$Adj_ratio)*$rall;
				my $exonrref=$hashAlab{$gene}{$exon}{trf};
				my $exonralt=$hashAlab{$gene}{$exon}{trs};
				my $exondref=$hashAlab{$gene}{$exon}{tdf};
				my $exondalt=$hashAlab{$gene}{$exon}{tds};
				my $exonDNAratio=$exondref/($exondref+$exondalt);
				my ($exonAdj_ratio,$exonRNAratio)=&Adj($exonDNAratio,$exonrref,$exonralt);	
				my $exonref=$exonAdj_ratio*($exonrref+$exonralt);
				my $exonalt=($exonrref+$exonralt)-$exonref;
				my @REF=($generef,$exonref);
				my @ALT=($genealt,$exonalt);
				my ($pCA,$od,$conf1,$conf2)=&fisherexact(\@REF,\@ALT);
				$Exon_gene=join(":",$exon,$pCA,$od,join("-",$conf1,$conf2));
				push(@pca,$Exon_gene);
			}
			push(@out,"\t".join("|",@pca)."\n");
		}else{push(@out,"\tNA\n");}
	}
}
if(!defined $opts{Y}){@sdall=@srall;}
if($Method ne 1){
	my ($fdr);
	if($SIMU > 0){
		my($Dall,$Dref,%simu,%hashloci,@cnts);
		my $i=0;
		foreach my $gene(keys %hashClab){
			push(@cnts,$hashClab{$gene}{counts});
		}
		my $Drefratio=$DNAfreq;
		my ($fdr1,$fdr2)=&Simulation(\@sdall,\@srall,\@srbq,\@cnts,$Drefratio,$Rphred,$SIMU);
		$fdr="##LLR cutoff for significant level of 0.05:".$fdr1." under ".$SIMU." times simulation; LLR cutoff for significant level of 0.01:".$fdr2." under ".$SIMU." times simulation\n";
		print OUT "P=0.05: ".$fdr1."\nP=0.01: ".$fdr2."\n";
	}else{
		$fdr="##Recommanded LLR cutoff for significant level of 0.05 : 0.82; for significant level of 0.01 : 1.48 \n";
	}
unshift(@out,$fdr);
	}
open(OUT, ">$output");
print OUT @out;	
}
	
#################################	
sub PhaseRead(){
	my ($file)=@_;
	my %phase;
	open(IN2,"$file");
	while(<IN2>){
		chomp;
		my @sep=split/\t/;
		my $name=join(":",$sep[0],$sep[1]);
		$phase{$name}{ref}=$sep[2];
		$phase{$name}{alt}=$sep[3];
	}
	return %phase;
}

sub ListRead(){
	my ($file)=@_;
	my %SNPa;
	open(IN1,"$file")or die $!;
	while (<IN1>){
		chomp;
		my @sep=split/\t/;
		my $name=join(":",$sep[0],$sep[1]);
		$SNPa{$name}=join("\t",$sep[2],$sep[3]);
		
	}
	return %SNPa;
}

sub RNARead(){
	my ($file,$r,$a,$d,$snp1)=@_;
	my %SNPl=%$snp1;
	my (%SNPr,$Total,$REF);
	open(IN2,"$file")or die $!;
	while (<IN2>){
		chomp;
		my @sep=split/\t/;
		my $name=join(":",$sep[0],$sep[1]);
		if(exists $SNPl{$name}){
			my @b=split/\t/,$SNPl{$name};
			my @Rallele_v=split//,$sep[4];
			my @Rbq_v=split//,$sep[5];
		    my($CurRefCnt,$CurNRefCnt,$OtherAlleleCnt,$all,$charout,$bqual)= &pileup_count($b[0],$b[1],\@Rallele_v,\@Rbq_v);
			if(($CurRefCnt >= $r)&&($CurNRefCnt >= $a)&&(($CurRefCnt+$CurNRefCnt) > 0)&($all >=$d)){
				$Total=$Total+$all;
				$REF=$REF+$CurRefCnt;				
				$SNPr{$name}=join("\t",$SNPl{$name},$CurRefCnt,$CurNRefCnt,$OtherAlleleCnt,$all,$charout,$bqual);
			}
		}
	}
	if ($Total==0){die "No SNVs on RNA-seq passed filter!\n";}
	my $RNAbias=$REF/$Total;
	return ($RNAbias,%SNPr);
}

sub DNARead(){
	my ($file,$R,$A,$D,$snp2,$snp1)=@_;
	my (%SNPf,%SNPd,$Total,$REF);
	my %SNPr=%$snp2;my %SNPlist=%$snp1;
	open(IN3,"$file")or die $!;
	while (<IN3>){
		chomp;
		my @sep=split/\t/;
		my $name=join(":",$sep[0],$sep[1]);
		if(exists $SNPlist{$name}){
			my @b=split/\t/,$SNPlist{$name};
			my @Dallele_v=split//,$sep[4];
			my @Dbq_v=split//,$sep[5];
			my ($CurRefCnt,$CurNRefCnt,$OtherAlleleCnt,$all,$charout,$bqual)= &pileup_count($b[0],$b[1],\@Dallele_v,\@Dbq_v);
			if(($CurRefCnt >= $R)&&($CurNRefCnt >= $A)&&(($CurRefCnt*$CurNRefCnt) > 1)&&($all >=$D)){
				$SNPd{$name}=join("\t",$CurRefCnt,$CurNRefCnt,$OtherAlleleCnt,$all,$charout,$bqual);
				$Total=$Total+$all;
				$REF=$REF+$CurRefCnt;	
			}
		}
	}
	if ($Total ==0){die "No SNVs on DNA-seq passed filter!\n";}
	my $DNAfreq=$REF/$Total;
	foreach my $keys(keys %SNPr){
		$SNPf{$keys}=join("\t",$SNPr{$keys},$SNPd{$keys});

	}		
	return ($DNAfreq,%SNPf);
}


sub pileup_count(){
	my($REF,$ALT,$pile,$pbq)=@_;
	my @ReadChar = @$pile;
	my @Readbq=@$pbq;
	my ($charout,$bqout,$OtherAllele);
	my ($Pos,$CurRefCnt,$CurNRefCnt,$all,$OtherAlleleCnt, $InsertCnt,$tempc)=(0,0,0,0,0,0,0,0);
	while($Pos <= $#ReadChar){
		my $Char = uc($ReadChar[$Pos]);
		if(($Char eq "\." )|| ($Char eq "\," )) {$CurRefCnt++;$Pos++;$charout=$charout.$REF;$bqout=$bqout.$Readbq[$tempc];$tempc++;}
		elsif($Char=~ /[atcgnACGTN]/) {
			$charout=$charout.uc($Char);
			$bqout=$bqout.$Readbq[$tempc];$tempc++;
			if($Char eq $ALT ){
				$CurNRefCnt++;
			}else {
				$OtherAllele .= $Char;
				$OtherAlleleCnt++;
			}
			$Pos++;
		} 
		elsif($Char eq "\$")	{$Pos++;}	# indicates that the end of a read maps to the prev pos
		elsif($Char eq "^")	{$Pos+=2;}	# beginning after a cigar operation
		elsif(($Char eq "+")||($Char eq "-"))  	{
			$Pos++; 
			$InsertCnt = "";
			while($ReadChar[$Pos] =~ /[0-9]/) {
				$InsertCnt .= $ReadChar[$Pos]; 
				$Pos++;
			}
			$Pos+= $InsertCnt; 
		}
		elsif(($Char eq ">")||($Char eq "<")||($Char eq "*")){ $Pos++;$tempc++;}
		else {print "unable to process character $ReadChar[$Pos]\n";}
	}
	$all=$CurRefCnt+$CurNRefCnt+$OtherAlleleCnt;
	return ($CurRefCnt,$CurNRefCnt,$OtherAlleleCnt,$all,$charout,$bqout);
}

sub theory_counts(){
	my($phred,$REF,$ALT,$alleles,$quality)=@_;		
	my ($i,$e,$dREF,$dALT,$fc)=(0,0,0,0,0);
	my @baseq=@$quality;
	my @alleles=@$alleles;
	for ($i=0;$i<=$#baseq;$i++){
		my $bq=ord($baseq[$i])-$phred;
		$e=10**(-$bq/10);
		if((uc($alleles[$i]) eq $REF)){
			$dREF=$dREF + (1-$e);
		}elsif(uc($alleles[$i]) eq $ALT){
			$dALT =$dALT + (1-$e);
		}
	}
	return($dREF,$dALT);
}
	
sub model(){
	use POSIX qw(log10);
	my($REF,$ALT,$allele,$bq,$mut,$mode)=@_;
	my @allele=split//,$allele;
	my @bq=split//,$bq;
	my ($p,$Likho,$e)=(0,0,0);
	if($#allele == $#bq){
		for(my $i=0;$i<=$#allele;$i++){
			$bq=ord($bq[$i])-$mode;
			$e=10**(-$bq/10);
			if($allele[$i] eq $REF){
				$p=log10($mut*$e/3+(1-$mut)*(1-$e));
			}elsif ($allele[$i] eq $ALT){
				$p=log10($mut*(1-$e)+(1-$mut)*$e/3);
			}else {
				$p=log10($e/3)
			}
			$Likho=$Likho+$p;
		}
	} else{print  "Invalid input!";}
	return $Likho;
}
	

sub M1(){
	my ($sta)=@_;
	my @stat=@$sta;
	my $dref=shift @stat; my $dalt=shift @stat;
	my $rref=shift @stat; my $ralt=shift @stat;
	my @ref=($dref,$rref);
	my @alt=($dalt,$ralt);	
	my $chisq=&siteChitest(\@ref,\@alt);
	return $chisq;
}

sub M2(){
	use POSIX qw(log10);
	my ($hash,$stat)=@_;
	my %hash=%$hash; my @stat=@$stat;
	my $Rphred= pop @stat;
	my $rmut=pop @stat; my $dmut=pop @stat;
	my ($LOD,$cnt)=(0,0);
	foreach my $snp(sort keys %hash){
		$cnt++;
		my @sep=split/\t/,$hash{$snp};
		my $RefA=$sep[0]; my $AltA=$sep[1];
		my $RNA_allele=$sep[2];
		my $Rbq=$sep[3];
		my $sLOD=&BBASED($RefA,$AltA,$dmut,$rmut,$RNA_allele,$Rbq,$Rphred);
		$LOD=$LOD+$sLOD;
		}
	my $aLOD=$LOD/$cnt;
	return $aLOD;
}	

sub Method(){
	my ($hash,$stat,$Method)=@_;
	my @stat=@$stat;my %hash=%$hash;
	my ($Re);
	if($Method eq "1"){
		my $chisq=&M1(\@stat);
		$Re=$chisq;
	}elsif($Method eq "2"){
		my $lod=&M2(\%hash,\@stat);
		$Re=$lod;
	}elsif($Method eq "3"){
		my $chisq=&M1(\@stat);
		my $lod=&M2(\%hash,\@stat);
		$Re=join("\t",$chisq,$lod);
	}
	return $Re;
}
	
sub Adj(){
	my ($DNAratio,$rref,$ralt)=@_;
	my $dratio=$DNAratio/(1-$DNAratio); 
	my $tAdj_ratio;
	my $RNAratio =$rref/($rref+$ralt);
	if($ralt ==0){$tAdj_ratio=1;}
	elsif($rref==0){$tAdj_ratio=0;}
	else{
		my $rratio=$rref/$ralt;
		my $tfc=$rratio/$dratio;
		$tAdj_ratio=$tfc/(1+$tfc);
		}
	return ($tAdj_ratio,$RNAratio);
}

sub Change(){
	my ($cnt)=@_;
	my @counts=@$cnt;
	@counts[0,1]=@counts[1,0];
	@counts[2,3]=@counts[3,2];
	if($#counts > 7){
	@counts[8,9]=@counts[9,8];
	}
	return @counts;
}

sub BBASED(){
	my($RefA,$AltA,$d_mut,$r_mut,$Rallele,$Rbq,$Rphred)=@_;	
	my $DNA_LLD=&model($RefA,$AltA,$Rallele,$Rbq,$d_mut,$Rphred);
	my $RNA_LLD=&model($RefA,$AltA,$Rallele,$Rbq,$r_mut,$Rphred);
	my $LLR=$RNA_LLD-$DNA_LLD;
	return $LLR;
}

sub siteChitest(){
	use Statistics::R;
	my ($ref,$alt)=@_ ;
	my ($p)=(0);
	my ($R);
	my @ref=@$ref;my @alt=@$alt;
	my $R = Statistics::R->new();
	$R->set('ref',\@ref); 
	$R->set('alt',\@alt);
	$R->run(q`df <- rbind(ref,alt);chi<-chisq.test(df)`);
	$p=$R->get('chi$p.value');
	$R->stop();
	return($p);
	}

sub heteroTest(){
	use Statistics::R;
	my ($ref,$alt)=@_ ;
	my ($rc,$ac,$p)=(0,0,0);
	my ($R);
	my @ref=@$ref;my @alt=@$alt;
	foreach my $refc(@ref){
		$rc=$rc+$refc;
	}
	foreach my $altc(@alt){
		$ac=$ac+$altc;
	}
	if(($rc==0)||($ac==0)){ $p =1;}
	else{
		my $R = Statistics::R->new();
		$R->set('ref',\@ref); 
		$R->set('alt',\@alt);
		$R->run(q`df <- round(rbind(ref,alt));chisq<-chisq.test(df)`);
		$p=$R->get('chisq$p.value');
		$R->stop();
		}
	return($p);
}
	
sub fisherexact(){
	use Statistics::R;
	my ($ref,$alt)=@_ ;
	my ($p,$odds,$conf1,$conf2)=(0,0,0,0);
	my ($R);
	my @ref=@$ref;my @alt=@$alt;
	my $R = Statistics::R->new();
	$R->set('ref',\@ref); 
	$R->set('alt',\@alt);
	$R->run(q`df <- round(rbind(ref,alt));fish<-fisher.test(df,conf.level=0.99)`);
	$p=$R->get('fish$p.value');
	$odds=$R->get('fish$estimate');
	$conf1=$R->get('fish$conf.int[1]');
	$conf2=$R->get('fish$conf.int[2]');
	$R->stop();
	return($p,$odds,$conf1,$conf2);
	}

sub Simulation(){
	use Statistics::R;
	my ($ddep,$rdep,$rbq,$cnt,$dfreq,$Rphred,$times)=@_;
	my @ddep=@$ddep;my @rdep=@$rdep,my @rbq=@$rbq;
	my (@LODD,%hash_sim,%hash_simlab);
	my $cnts=&rand($times,$cnt);
	my @cnts=@$cnts;
	my $length=$#ddep;
	for (my $i=0;$i<$times;$i++){
		my ($DD,$RD,$Df,$Rf)=(0,0,0,0);
		for(my $j=1;$j<=$cnts[$i];$j++){
			my $samp=int(rand($length));
			my $Ddepth=$ddep[$samp];
			my $Rdepth=$rdep[$samp];
			my $Rbq=$rbq[$samp];
			$Rbq=~s/[\\]//g,$Rbq;
			my $Dref=&binorm($Ddepth,$dfreq,1);
			my $rfreq=$Dref/$Ddepth;
			my $Rref=&binorm($Rdepth,$rfreq,1);
			my $Ralt=$Rdepth-$Rref;
			$DD=$DD+$Ddepth;$RD=$RD+$Rdepth;
			$Df=$Df+$Dref;$Rf=$Rf+$Rref;
			my $R = Statistics::R->new();
			$R->set( 'Refcnt', $Rref );	
			$R->set( 'Altcnt', $Ralt );
			$R->set( 'bq', $Rbq );
			$R->set( 'Rdepth', $Rdepth );			
			$R->run(q`REF<-c(rep("A",2997),"G","T","C");
					ALT<-c(rep("G",2997),"T","C","A");
					Allele<- paste(c(sample(REF,Refcnt,replace = T),sample(ALT, Altcnt,replace = T)),collapse="");
					bq_v<-unlist(strsplit(bq, "",fixed = TRUE));
					Rbq<-sample(bq_v,Rdepth,replace=T);`);
			my $allele=$R->get('Allele');
			my $bq=$R->get('Rbq');
			my @bq=@$bq;
			my $rbq=join("",@bq);
			$hash_sim{$i}->{$j}=join("\t","A","G",$allele,$rbq);
		}
		my $dmut=($DD-$Df)/$DD;
		my $rmut=($RD-$Rf)/$RD;
		$hash_simlab{$i}->{dmut}=$dmut;
		$hash_simlab{$i}->{rmut}=$rmut;
	}		
	foreach my $gene (keys %hash_simlab){
		my $DNAmut=$hash_simlab{$gene}{dmut};
		my $RNAmut=$hash_simlab{$gene}{rmut};
		my @si=($DNAmut,$RNAmut,$Rphred);
		my %hashm=%{$hash_sim{$gene}};
		my $LLR=&Method(\%hashm,\@si,2);
		if($LLR >0){
		push(@LODD,1/$LLR);
		}
	}
	my $R = Statistics::R->new();
	$R->set("LOD",\@LODD);
	$R->run(q`qt <-quantile(LOD,c(0.05,0.01));fdr1=qt[1];fdr2=qt[2];`);
	my $fdr1=$R->get('fdr1');
	my $fdr2=$R->get('fdr2');
    return(1/$fdr1,1/$fdr2);
}

sub rand(){
	use Statistics::R;
	my ($c,$array)=@_;
	my @arry=@$array;
	my $R = Statistics::R->new();
	$R->set('key',\@arry); 
	$R->set('time',$c);	
	$R->run(q`samp <- sample(key,time,replace=T)`);
	my $p=$R->get('samp');
	$R->stop();
    return $p;
}
	
sub binorm(){
	use Statistics::R;
	my ($dep,$freq,$n)=@_;
	my $R = Statistics::R->new();
	$R->set("dep",$dep);
	$R->set("f",$freq);
	$R->set("n",$n);
	$R->run(q`ref <-rbinom(n,dep,f)`);
	my $ref=$R->get('ref');
	return ($ref);
}
	
sub usage {
  die(qq/
Version : 1.0.2beta (2016-08-15)
Contact: Zhi Liu  <liuzhi2011\@sibs.ac.cn>\n
Usage:   cisASE.pl <command> [<arguments>]\n
Command: SNV       SNV level ASE detection
         EXON      EXON level ASE detection
         GENE      GENE level ASE detection

\n/);
}
	
