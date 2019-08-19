#!/usr/bin/env perl
=head1 USAGE

	perl esomCodonMod.pl -lrn file.lrn -o outputFile.lrn

=cut

use strict;
use Getopt::Long;

my ($lrn, $tri);
my $out=$$."_".$lrn;
my $version=0.2.0;

GetOptions(
	'lrn=s'=>\$lrn,
	'o|out:s'=>\$out,
	'tri'=>\$tri,
	'v|version'=>sub{print STDERR $0."\tversion:".$version."\n";},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

&help  if (! $lrn);
sub help{system('perldoc', $0); exit;}

my @removeCodons=qw (ATG TAG TAA TGA);
my @nucl=qw(A T G C);

my %removeTetra;
foreach my $c(@removeCodons){
    foreach my $n(@nucl){
	$removeTetra{$n.$c}++;
	$removeTetra{$c.$n}++;		
    }
}
if($tri){	
    foreach (@removeCodons){	$removeTetra{$_}++; }
}

#print "Possible Tetramers that can be Removed:\t".keys(%removeTetra)."\n";

open(LRN, $lrn) || die $!;
my (@codonOrder);
my ($cols, $secondPart, $firstLine, $removed);
while(my $line=<LRN>){
    chomp $line;
    next unless $line;
    if ($line=~ /^\% Key/){
	@codonOrder=split(/\t/, $line);
	my @thisLine;
	foreach (@codonOrder){
	    if ($removeTetra{$_}){
		$removed++;
		next;
	    }
	    push @thisLine, $_;
	    $cols++;
	}
	$secondPart.= join("\t",@thisLine);
    } elsif($line=~ /^\d/){
	my @thisLine;
	my @frequencies=split(/\t/, $line);
	my $pos=-1;
	foreach my $freq(@frequencies){
	    $pos++;
	    next if ($removeTetra{$codonOrder[$pos]});
	    push @thisLine, $freq;
	}
	$secondPart.= join("\t",@thisLine);
    } elsif($.==1){
	$firstLine=$line;
    }
}
close LRN;

print "Tetramers Removed:\t".$removed."\n";

open(OUT, ">".$out) || die $!;
print OUT "$firstLine\n","% $cols\n";

my @thisLine ="% 9"; # can fix this to use an array instead of string append
for (my $i=$cols; $i > 1; $i--) {
    push @thisLine, "1";
}

print OUT join("\t",$thisLine),"\n";
print OUT $secondPart;
close OUT;

exit 0;
