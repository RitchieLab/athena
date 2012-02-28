#!/usr/bin/perl
use strict;
use Getopt::Std;

my $PRGM_NAME = 'makeathenaconfigs.pl';

# Subroutine: usage
# Displays current version date and usage
sub usage{
  print "\nusage:\t\t$PRGM_NAME <arguments>\n\n";
  print "example:\t$PRGM_NAME -p Sample_Dataset- -s .txt -d /my/full/path/to/sample_sets -a GENN -g /my/full/path/to/genn_add.gram\n\n";
  print "\tOptions:\n";
  print "\t\t-p Prefix for datasets\n";
  print "\t\t-s Suffix for datasets\n";
  print "\t\t-d Directory containing datasets\n";
  print "\t\t-a Algorithm to use (GENN or GESR, must match grammar!)\n";
  print "\t\t-g Full path to grammar file (must match algorithm!)\n";
  print "\t\t-o Output directory (optional -- set to dataset directory when omitted)\n\n";
}

# subroutine: check_args
# Checks number of arguments passed in and returns name of pedigree file.
sub check_args{
  my %opts;

  getopts('p:s:d:a:g:o:', \%opts);
  # Checks if atleast p,s,g,d arguments are present
  if (!$opts{'p'} or !$opts{'s'} or !$opts{'d'} or !$opts{'a'} or !$opts{'g'}){
    usage();
    exit();
  }
  # Checks for output argument
  if(!$opts{'o'})
  {
    $opts{'o'} = $opts{'d'};
  }

  return $opts{'p'}, $opts{'s'}, $opts{'d'}, $opts{'a'}, $opts{'g'}, $opts{'o'};
}


my ($prefix, $suffix, $datadirectory, $algo, $grammarfile, $outdirectory) = check_args();
print "outdirectory=$outdirectory\n";

#This script will create 100 configuration files that correspond to 100 genomeSIMLA datasets numbered 1-100 

# Write out the full path to the datasets before the number in the dataset name for prefix and after the number in the name for suffix

for (my $i = 1; $i<=100; $i++){
	open (CFG,">$outdirectory/$i.config");
	my $dataset = "$datadirectory/$prefix$i$suffix";
        my $outfile = "$outdirectory/$prefix$i$suffix";
	print CFG "
# ATHENA CONFIG FILE
#Algorithm specifies whether GESR or GENN will be run. Need to change grammarfile accordingly.
";
	if ($algo eq "GENN") {
		print CFG "
ALGORITHM GENN
MINSIZE 20
MAXSIZE 2000
MAXDEPTH 10
SENSIBLEINIT FALSE
POPSIZE 100
PROBCROSS 0.9
PROBMUT 0.01
";
	} elsif ($algo eq "GESR") {
		print CFG "
ALGORITHM GESYMBREG
MINSIZE 20
MAXSIZE 2000
MAXDEPTH 10
SENSIBLEINIT FALSE
POPSIZE 100
PROBCROSS 0.9
PROBMUT 0.01
GASELECTION DOUBLE 
DOUBTOURNF 7
DOUBTOURND 1.4
DOUBTOURNFITFIRST TRUE
";
	} else {
		print CFG "
ALGORITHM $algo
";
	}
	print CFG "
#write out the full path to the grammar file according to whether you are running GESR or GENN (grammar file shown for GENN)
GRAMMARFILE $grammarfile

CALCTYPE BALANCEDACC
#CALCTYPE RSQUARED
EFFECTIVEXO TRUE
GENSPERSTEP 25
INCLUDEALLSNPS TRUE
BLOCKCROSSGENS 0
END GENN

DATASET $dataset
IDINCLUDED FALSE
MISSINGVALUE -9
CONTINMISS -9
DUMMYENCODE STEPHEN
RANDSEED $i

#write out path to directory where you want the results to be written
OUT $outfile

CV 5
NUMSTEPS 4
WRITECV FALSE
SUMMARYONLY FALSE
STATUSADJUST NORMMAX
";
	close (CFG);
}
