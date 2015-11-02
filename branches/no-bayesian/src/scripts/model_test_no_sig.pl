#!/usr/bin/perl

use strict;

if(@ARGV < 2){
  print "model_test.pl train.txt test.txt\n";
  exit;
}

sub read_data{
  my $filename = shift;
  my @data;
  open(IN, $filename) or die "$filename:  $!\n\n";
  while(<IN>){
    chomp;
    my @info = split" ";
    push(@data, \@info);
  }
  close(IN);
  return @data;
}

#PM ( W ( ( ( 54.93 + 86.64 ) + ( ( 58.47 / 35.01 ) / 2.53 ) ) , V30 ) , W ( 1.79 , V30 ) , 2 ) 
#sub evaluate_ind{
#  my @ind_array = @_;
#  my $firstnode = 142.23 * $ind_array[29];

#  my $secondnode = 1.79 * $ind_array[29];

#  return ActivateSigmoid($firstnode * $secondnode);
#}

#PA ( W ( 0.0073 , V4 ) , W ( ( 0.0048 / 0.005 ) , PS ( W ( 0.0091 , V17 ) , W ( 0.0095 , V9 ) , 2 ) ) , 2 )
# sub evaluate_ind{
#   my @ind_array = @_;
#   my $firstnode = .0073 * $ind_array[3];
# 
#   my $secondnode = .0091 * $ind_array[16];
#   my $thirdnode = .0095 * $ind_array[8];
#   $secondnode = ActivateSigmoid(.0048/.005 * ($secondnode - $thirdnode));
#   
# #  if($secondnode == 0){
# #    return ActivateSigmoid(1.1);
# #  }
# #  else{
#     return ActivateSigmoid($firstnode+$secondnode);
# #  }
# 
# }


#PNoSigA( W(2.19,G5), W((0-0.61),G260),2)
sub evaluate_ind{
   my @ind_array = @_;
   my $firstnode = 2.19* $ind_array[5];
   my $secondnode = -0.61 * $ind_array[260];
   return $firstnode + $secondnode;
}


#PA ( W ( 2.05 , V22 ) , W ( 60.75 , PD ( W ( 26.58 , V30 ) , W ( ( 41.14 / 74.8 ) , V31 ) , W ( ( 62.19 + ( 62.01 / 35.61 ) ) , V12 ) , W ( ( 56.27 + 74.25 ) , V2 ) , 4 ) ) , 2 )
#sub evaluate_ind{
#  my @ind_array = @_;
#  my $firstnode = 2.05 * $ind_array[21];
  
#  my $secondnode;
  # compile the four nodes together
#  my $node1 = 26.58 * $ind_array[29];
#  my $node2 = 0.55 * $ind_array[30];
#  my $node3 = 63.931 * $ind_array[11];
#  my $node4 = 130.52 * $ind_array[1];

#  if($node2==0 or $node3==0 or $node4==0){
#   $secondnode = 1.1; 
#  }
#  else{
#    $secondnode = $node1/($node2*$node3*$node4);
#  }
  

#  $secondnode = 60.75 * ActivateSigmoid($secondnode);
  
#  return ActivateSigmoid($firstnode + $secondnode);
#}


sub ActivateSigmoid{
  my $x = shift;
  
  if($x < -709){ return -1.0;}
  if($x > 709){return 1.0;}
  return(1.0 / (1.0 + exp(-$x)));
}

my @data = read_data($ARGV[0]); # read training data
my $traintotal = 0;
my $trainerror = 0;

my $square_error_total=0;

foreach my $array(@data){
  my $sc = evaluate_ind(@$array);
  $traintotal++;
  my @ind = @$array;
  my $err = $ind[$#ind] - $sc;
  $square_error_total += $err * $err
}

print "square_error_total=$square_error_total\n";
print "training mean square error: ",$square_error_total/$traintotal, "\n";

my @data = read_data($ARGV[1]); # read training data
my $testtotal = 0;
my $testerror = 0;
$square_error_total = 0;
foreach my $array(@data){
  my $sc = evaluate_ind(@$array);
  my @ind = @$array;
print "$$array[19] $$array[8] score is $sc ";
print "status for $ind[$#ind] \n";
  $testtotal++;
  my $err = $ind[$#ind]- $sc;
#printf "adding %g\n", $err*$err;
  $square_error_total += $err * $err
}
print "square_error_total=$square_error_total\n";
print "testing mean square error:  ",$square_error_total/$testtotal, "\n";

#for(my $i=-100; $i<50; $i+=10){
#  printf "$i =>%e\n", ActivateSigmoid($i);
#}


