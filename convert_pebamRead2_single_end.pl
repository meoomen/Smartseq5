#!/usr/bin/perl -w
use English;
use warnings;
use strict;

sub outputWrapper($;$) {
    # required
    my $outputFile=shift;
    # optional
    my $outputCompressed=0;
    $outputCompressed=shift if @_;
               
    $outputCompressed = 1 if($outputFile =~ /\.gz$/);
    $outputFile .= ".gz" if(($outputFile !~ /\.gz$/) and ($outputCompressed == 1));
    $outputFile = "| gzip -c > '".$outputFile."'" if(($outputFile =~ /\.gz$/) and ($outputCompressed == 1));
    $outputFile = ">".$outputFile if($outputCompressed == 0);
               
    return($outputFile);
}
 

sub inputWrapper($) {
     my $inputFile=shift;
               
     $inputFile = "samtools view -Sh '".$inputFile."' | " if(($inputFile =~ /\.bam$/) and (!(-T($inputFile))));
               
     return($inputFile);
}


print "hi - start of the program\n";
print "\n";


# get the input file
my $samFile=$ARGV[0];

print "samFile = $samFile\n";

open(IN,inputWrapper($samFile)) or die "cannot open < $samFile - $!";

open(OUT,outputWrapper($samFile.'.single_end.sam.gz'));

my $lineNum=0;

my %chrCounterHash=();
my %lengthCounterHash=();


while(my $line = <IN>) { 
	chomp($line);
	
 	if($line =~ /^@/) {
		print OUT "$line\n";
		next;
	}
	
	print "$lineNum ... \n" if(($lineNum % 100000) == 0);
	$lineNum = $lineNum + 1;

	
	my @array=split(/\t/,$line);
	my $numColumns=@array;
	#print "\tfound $numColumns columns\n";
	
	my $readID=$array[0];
	my $flag=$array[1];
	my $flag_new= $flag - 131;


	print OUT "$readID\t$flag_new\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t$array[9]\t$array[10]\n";
		
		

		
}	

close(IN);
close(OUT);

print "done\n";
