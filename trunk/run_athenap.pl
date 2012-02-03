#!/usr/bin/perl

#This script will run athenap on all athena configuration files specified. 
#Write command in the pbs script as /path/to/run_athenap.pl /path/to/configs/*.configname
#IMPORTANT: Make sure -np # matches the number of processors you request in the .pbs script. This dictates number of demes for athenap.

@configs = @ARGV;
foreach $currentconfig (@configs) {
system("mpiexec -np 4 /path/to/athena/parallel/executable/ $currentconfig");
}
