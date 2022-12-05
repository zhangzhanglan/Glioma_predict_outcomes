use strict;
use warnings;

## https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap
## https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_fpkm.gz
my $gtfFile="probeMap_gencode.v23.annotation.gene.probemap";
my $expFile="gtex_RSEM_gene_fpkm.gz";
my $outFile="GTExSymbol.txt";

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

open(RF, "gzip -dc $expFile|") or die $!;
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
                for(my $i=1;$i<=$#arr;$i++){
                        my $fpkm=2**$arr[$i]-0.001;
                        if($fpkm<0){
                                $fpkm=0;
                        }
                        $arr[$i]=log($fpkm+1)/log(2);
                }
                print WF join("\t",@arr) . "\n";
        }
}
close(WF); 
close(RF);


my %hash=();
my $normalFlag=0;
## https://toil.xenahubs.net/download/GTEX_phenotype.gz
my $sampleFile="GTEX_phenotype.gz";

open(RF,"gzip -dc $sampleFile|") or die $!;
while(my $line=<RF>){
    chomp($line);
    next if $line =~ /^Sample/;
    my @line = split /\t/, $line;
    if ($line[2] eq "Brain") {
        # print $line[0], "\n";                                                                                         
        $hash{$line[0]}=1;
    }
}
close(RF);

my @indexs=();
open(RF,"GTExSymbol.txt") or die $!;
open(WF,">GTExNormalBrainExp.txt") or die $!;
my @samples=();
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		for(my $i=1;$i<=$#arr;$i++){
					my $sampleName=$arr[$i];
					if(exists $hash{$sampleName}){
						  push(@indexs,$i);
						  push(@samples,$arr[$i]);
						  delete($hash{$sampleName});
					}
			}
			print WF "ID\t" . join("\t",@samples) . "\n";
	}
	else{
			my @sampleData=();
			foreach my $col(@indexs){
				  push(@sampleData,$arr[$col]);
			}
			print WF "$arr[0]\t" . join("\t",@sampleData) . "\n";
	}
}
close(WF);
close(RF);

print "sample: " . ($#samples+1) . "\n";
