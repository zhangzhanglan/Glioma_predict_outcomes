use strict;
use warnings;

#!/usr/bin/perl -w
use strict;
use warnings;

# https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-GBM.htseq_fpkm.tsv.gz
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LGG.htseq_fpkm.tsv.gz
# my $file = "TCGA-GBM.htseq_fpkm.tsv.gz";
my $file = "TCGA-LGG.htseq_fpkm.tsv.gz";
my %hash=();
my @normalSamples=();
my @tumorSamples=();
my @sampleArr=();

open(RF,"gzip -dc $file|") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		for(my $i=1;$i<=$#arr;$i++){
			my @idArr=split(/\-/,$arr[$i]);
			if($idArr[3]=~/^0/){
				push(@tumorSamples,$arr[$i]);
			}
			else{
				push(@normalSamples,$arr[$i]);
			} 
		}
		@sampleArr=@arr;
		next;
	}
	else{
		for(my $i=1;$i<=$#arr;$i++){
			${$hash{$arr[0]}}{$sampleArr[$i]}=$arr[$i];
		}
	}
}
close(RF);

open(WF,">tcgaEnsembl.txt") or die $!;
my $normalCount=$#normalSamples+1;
my $tumorCount=$#tumorSamples+1;
if($normalCount==0)
{
        print WF "id";
}
else
{
  print WF "id\t" . join("\t",@normalSamples);
}
print WF "\t" . join("\t",@tumorSamples) . "\n";
foreach my $key(keys %hash)
{
        print WF $key;
        foreach my $normal(@normalSamples)
        {
                print WF "\t" . ${$hash{$key}}{$normal};
        }
        foreach my $tumor(@tumorSamples)
        {
                print WF "\t" . ${$hash{$key}}{$tumor};
        }
        print WF "\n";
}
close(WF);

print "normal count: $normalCount\n";
print "tumor count: $tumorCount\n";

# https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap
my $gtfFile="gencode.v22.annotation.gene.probeMap";
my $expFile="tcgaEnsembl.txt";
# my $outFile="tcgaSymbol_GBM.txt";
my $outFile="tcgaSymbol_LGG.txt";

my %hash_da=();
open(RF,"$gtfFile") or die $!;
while(my $line=<RF>)
{
        chomp($line);
        next if $line =~ /id/;
        my @arry = (split /\t/, $line);
        $hash_da{$arry[0]}=$arry[1];
}

close(RF);

open(RF,"$expFile") or die $!;
open(WF,">$outFile") or die $!;
while(my $line=<RF>)
{
	if($.==1)
	{
		print WF $line;
		next;
	}
	chomp($line);
	my @arr=split(/\t/,$line);
	# $arr[0]=~s/(.+)\..+/$1/g;
	if(exists $hash_da{$arr[0]})
	{
		$arr[0]=$hash_da{$arr[0]};
		print WF join("\t",@arr) . "\n";
	}
}
close(WF); 
close(RF);
