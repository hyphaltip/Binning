#!/usr/bin/perl

=head1 DESCRIPTION

	Give me:
		- Path to folders that has the fasta files.

	I will:
		- Create an annotation file and a concatenated fasta for ESOM binning;
		- Run the tetramer frequency script on the files;

=head2 Dependencies

	As such there are no external perl module dependancies but this script is a wrapper which formats inputs runs other scripts, they are:
	tetramer_freqs_esom.pl - To calculate the tetramer frequencies of your contigs [ REQUIRED ]
	esomCodonMod.pl - To remove the tetramers containing stop codons from the analysis [ OPTIONAL ]

=head1 USAGE

	perl esomWrapper.pl -path Folder_Path -ext extension_of_files

=head2 Options

	-p or path	[characters]	path to folder containing fasta files; use "." (dot, without the quotes) for current folder.
	-e or ext	[characters]	file extension to look for in folder; default= fasta
	-prefix		[characters]	prefix filename for annotation and concatenated file; default=esom
	-scripts	[characters]	location of all the wrapped scripts; default="current working directory" OR "/geomicro/data1/COMMON/scripts/wrappers/ESOM/"
	-DIR or dir	[characters]	name of the output directory; default= ESOM
			
	-min	[integer]	Optional	default=2500; Minimal length (in nt) of input contig to be included in output
	-max	[integer]	Optional	default=5000
	Note:	The script will split sequence after each 'max' nt; join last part, if remaining seq shorter than 'max', with second-last part
			eg: in default settings, a sequence of 14 kb will be split into a 5 kb and a 9 kb fragment if window_size = 5 kb.
			
	-mod	[BOOLEAN]	[Experimental Feature] Use the codon mod script which removes the kmers containing stop-codons.(output=*.mod.lrn)
	-h	this page.

=head3 Example 1: Required Options

	perl esomWrapper.pl -path .

=head3 Example 2: Other Options

	perl esomWrapper.pl -path . -ext fa -dir MyESOM -prefix esomOutput -min 2000 -max 5000

=head1 Suggestions/Corrections/Feedback/Beer

	Sunit Jain, sunitj@umich.edu
	January 2013

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Spec;
use File::Basename;
use Pod::Usage;
use Bio::SeqIO;
#use POSIX ":sys_wait_h";

my $scripts = $Bin;
my $version=join("\t",qw(esomWrapper.pl v0.2.91));
my $path; # Fasta Folder path
my $ext="fasta";
my $prefix="esom";
my $outDir="ESOM";
my $kmer = 4;
my $mod;
my $force = 0;
my $train;
my $info;
my $min_length = 2500; #Minimal length (in nt) of input contig to be included in output
my $window_size = 5000; #split sequence after each window_size nt, 
                     #join last part, if shorter than window_size, 
                     #with second-last part (a sequence of 
                     #14 kb will be split into a 5 kb and a
                     #9 kb fragment if window_size = 5 kb)

my $man = 0;
my $help = 0;

GetOptions(
    'p|path:s'   =>  \$path,
    'e|ext:s'    =>  \$ext,
    'k|kmer:i'   =>  \$kmer,
    'prefix:s'   =>  \$prefix,
    'DIR|dir:s'  =>  \$outDir,
    'min:i'      =>  \$min_length,
    'max:i'      =>  \$window_size,
    'mod'        =>  \$mod,
    'train:s'    =>  \$train,
    'info:s'     =>  \$info,
    'scripts:s'  => \$scripts,
    'h|help|?'   => \$help, 
    'man'        => \$man,
    'force!'     => \$force,
    'v|version'  => sub{print "#".$version."\n"; exit;}
    );

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

print "## $version ##\n";
die "[ERROR: $0] Folder Path Required! See $0 -h for help on the usage" if !$path;

opendir(PATH,$path) || die "cannot open directory $path: $!";
if( $ext =~ s/^\.//) {
    warn("stripping leading . from extension since will be added later");
}
my @files;
foreach my $file ( readdir(PATH) ) {
    next unless $file =~ /\.\Q$ext\E$/;
    push @files, File::Spec->catfile($path,$file);
}

unless( scalar @files) {
    die ("[ERROR] Can't find any files with .$ext in \"$path\"\n".
	 "Please check that the path exist or that you have sufficient privilages.\n") 
}

my $annotationFile   = $prefix.".ann";
my $concatenatedFasta= $prefix.".".$ext;
my $logFile          = $prefix.".log";

if (-d $outDir && ! $force ){
    die "[ERROR: $0]$outDir already exists!\n";
} else{
    mkdir($outDir, 0755);
}

my $ann     = File::Spec->catfile( $outDir, $annotationFile);
my $catFasta= File::Spec->catfile( $outDir, $concatenatedFasta);
my $log     = File::Spec->catfile( $outDir, $logFile);

open(LOG, ">$log") || die "cannot open $log for writing: $!";

if ( ! $scripts || ! -d $scripts ){
    if(-e "tetramer_freqs_esom.pl"){
	$scripts=`pwd`;
	chomp $scripts;
    } else{
	die "[ERROR: $0] Could not locate helper scripts: 'tetramer_freqs_esom.pl', 'esomCodonMod.pl' and 'esomTrain.pl', please provide the location using '-scripts' flag\n";
    }
}

print "# Setting scripts folder location as:\n$scripts\n";
print LOG "# Setting scripts folder location as:\n$scripts\n";

### Wrapped Scripts ###
my $tetramerScript= File::Spec->catfile($scripts,"tetramer_freqs_esom.pl");
my $codonModScript= File::Spec->catfile($scripts,"esomCodonMod.pl");
my $esomTrain     = File::Spec->catfile($scripts,"esomTrain.pl");
###

unless (-e $tetramerScript) {
    die "Can't find the required scripts, please provide the location using the '-scripts' flag\n";

}
#$|++;

open(my $fastafh => ">$catFasta") || die "cannot open $catFasta for writing: $!";
open(my $annfh   => ">$ann") || die "cannot open $ann for writing: $!";

my $class = 0;
my $filesProcessed=0;

warn "# FileName\tNumber of Sequences found\n";
print LOG "# FileName\tClass Assigned\tNumber of Sequences\n";
print $annfh "# Contig\tAnnotation\tClass\n";
my %cls;
foreach my $file (@files) {
	my $countSeqs =	parseFasta($file,$fastafh,$annfh);

	my $fileName = basename($file, ".".$ext);
	$fileName =~ s/\s+/\_/g;

	warn $fileName."\t".$countSeqs."\n";
	print LOG $fileName."\t".$class."\t".$countSeqs."\n";
	
	$cls{$class++} = $fileName;
	$filesProcessed++;
}
close($fastafh);
close($annfh);

warn "\n# Files processed:\t $filesProcessed\n";
print LOG "\n# Files processed:\t $filesProcessed\n";
close LOG;

print "# Calculating Tetramer Frequencies...\n";
chdir $outDir || die $!;
$log = $logFile;
open (LOG2, ">>".$log ) || die "$log not found\n$!\n";

system("perl $tetramerScript -f $concatenatedFasta -a $annotationFile -min $min_length -max $window_size -ext $ext -kmer $kmer >> $log");

my $lrnfile ="Tetra_".$prefix."_".$min_length.".lrn";
my $modLrnFile;
if($mod) {
    if(! -e $codonModScript){
	warn "[WARNING:] $codonModScript not found. Please use the '-scripts' flag to specify the script's location\n";
	warn "[WARNING:] The rest of the analysis will be done on the unmodified \*.lrn file.\n";
	$mod--;
    } else{
	warn "# Applying Codon Modification...\n";
	$modLrnFile = "Tetra_".$prefix."_".$min_length.".mod.lrn";
	system("perl $codonModScript -lrn $lrnfile -o $modLrnFile >> $log");
    }
}

print "# Adding class names and colors to the cls file\n";
my $clsFile= "Tetra_".$prefix."_".$min_length.".cls";
my $tmpCls="tmp.cls";
open(CLS, $clsFile)|| die "cannot open $clsFile: $!";
open(TMP, ">$tmpCls")|| die "cannot open $tmpCls for writing: $!";

while(my $line = <CLS>) {
    if ($. == 2 ){
	for(my $i=0; $i < $filesProcessed; $i++){
	    my $clsColor = randomColors();
	    print TMP "\%".$i."\t".$cls{$i}."\t".$clsColor."\n";
	}
    }
    print TMP $line;
}
close(CLS);
close(TMP);
# this should be using File::Copy qw(move) I think to be system indepdendent
# also need to close the files to be written before moving!
system("mv $tmpCls $clsFile");

=begin Training Commented
print "# Let the Training begin...\n";
my @dimensions=`grep "^>" $log`;
my ($rows, $cols);
foreach(@dimensions){
	chomp;
	my @blah=split(/\t/, $_);
	if(lc($_)=~ /^>rows/){
		$rows=$blah[-1];
#		print $rows."\n";
	}
	elsif(lc($_)=~ /^>cols/){
		$cols=$blah[-1];
#		print $cols."\n";	
	}
	else{
		print LOG2 "[ERROR] Script Borked! Try running the esomTrain.pl script independently\n";
	}
}
my %PIDs;
my $modLog="modTrain.log";

if (-e $esomTrain){
	if ($rows && $cols && $train){
		if($noMod){
			my $command="perl $esomTrain -lrn $lrnfile -cls $clsFile -rows $rows -columns $cols -norm $train".($info ? " -info $info" : "");
#			print LOG2 $command."\n";
			system ("$command >> $log");
		}
		else{
			my $command="perl $esomTrain -lrn $modLrnFile -cls $clsFile -rows $rows -cols $cols -norm $train".($info ? " -info $info" : "");
### Run this bit on a seperate thread.
#			print LOG2 "[Thread2:] ".$command."\n";
#			my $pid=&run("$command > $modLog");
#			$PIDs{$pid}++;
			system ("$command >> $log");			
		}
	}
}
else{
	print LOG2 "esomTrain.pl not found. Please specify the location of the script using the '-scripts' flag or try running the esomTrain.pl script independently\n";
}
=cut

#if (keys %PIDs){
#	&REAP;
#	print LOG2 "\n\n################## MOD FILE TRAINING LOG ##################\n\n";
#	system("cat $modLog >> $log");
#}
#unlink $modLog;

warn "\nAll done! please check the $log file for errors/warnings before you proceed\n";
warn "Also make sure that your class (.cls) files have values in both the columns\n";
print LOG2 "\nAll done! please check this file for errors/warnings before you proceed\n";
print LOG2 "Please make sure that your class (.cls) files have values in both the columns\n";
close LOG2;
exit 0;

#### Sub-routines ####
# Could be replaced with BioPerl ... Bio::SeqIO
sub parseFasta {
    my ($fileName,$fastafh,$annfh) = @_;
    my $in = Bio::SeqIO->new(-format => 'fasta',
			     -file   => $fileName);
    my $countSeqs = 0;
    while( my $seq = $in->next_seq ) {
	my $header = $seq->display_id;
	# if really want whole FASTA header then
	# my $header =  join(" ",$seq->display_id,$seq->description);
	my $seq    = $seq->seq;
	
	$header =~ s/[\s\W]+/\_/g;
	$header =~ s/\_+/\_/g;
	$header =~ s/\_+$//;
	$header =~ s/^\_+//;    
	$countSeqs++;
	print $fastafh ">".$header."\n".$seq."\n";
	print $annfh $header."\t".$header."\t".$class."\n";
    }
    return $countSeqs;
}

sub randomColors {
    my ($r, $g, $b) = map { int rand 256 } 1 .. 3;
    my $color= join("\t", $r, $g, $b);
    return ($color);
}

__END__

sub run{
	my $command=shift;
	my $pid = fork();

	if (!defined($pid)) {
    	die "unable to fork: $!";
	}
	elsif ($pid==0) { # child
		print "Executing:\t$command\n";
		exec($command) || die "unable to exec: [$?]\n$!\n";
		exit(0);
	}
	# parent continues here, pid of child is in $pid
	return($pid);
}

sub REAP{ ## Use this when you want to wait till the process ends before further processing.
	my $numPIDs= scalar(keys %PIDs);

#	print "in REAPER: ".$numPIDs."\n";
	while (scalar(keys %PIDs) > 0){
		my $pid= waitpid(-1, &WNOHANG);
		if ($pid > 0){
			print "in REAPER:$pid\n";
			if (WIFEXITED($?) && $PIDs{$pid}){
				`echo "Process ID: $pid\tFinished with status $?"`;
#				$numPIDs-- ;
				print "Process: ".$pid."\tStatus: ".$?."\nWaiting for ".$numPIDs." more processes...\n";
				delete $PIDs{$pid};
			}
		}
		else{
			sleep 30;
		}
	}
	return;
}

