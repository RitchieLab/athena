#!/usr/bin/perl

use strict;
use Math::Trig;

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
#PA( W((0.75-8.81),G1), W((8-6),G17), W((0.75/9),G19),3)
# sub evaluate_ind{
#   my @ind_array = @_;
#   my $firstnode = (0.75-8.81) * $ind_array[1];
#   my $secondnode = (8-6) * $ind_array[17];
#   my $thirdnode = (0.75/9)* $ind_array[19];
#   my $result = $firstnode + $secondnode + $thirdnode;
#   
# print " $ind_array[1] $ind_array[17] $ind_array[19] => ", ActivateSigmoid($result);
# # print "\n";
#   
# #   my $result;
# #   if($secondnode == 0 || $thirdnode == 0){
# #     $result = 10000000;
# #   }
# #   else{
# #     $result = ($firstnode / $secondnode) / $thirdnode;
# #   }
# 
#   return ActivateSigmoid($result);
#   
# }

#PA( W(0,G5), W((0.35-(0.16-(3.8-4.4))),C2),2)
#PD( W((0.43-9.91),G175), W((7.4/(((6*5.63)+2)+(0.55-(8*4.5)))),G20),2)

# sub evaluate_ind{
#   my @ind_array = @_;
#   my $firstnode = (0.43-9.91)* $ind_array[175];
#   my $secondnode = (7.4/(((6*5.63)+2)+(0.55-(8*4.5)))) * $ind_array[20];
#   
#   my $result;
#   if($secondnode == 0){
#     $result = 10000000;
#   }
#   else{
#     $result = $firstnode / $secondnode;
#   }
#   
#  #  my $result = $firstnode / $secondnode;
#   
# # print " $ind_array[5] $ind_array[17] $ind_array[19] => ", ActivateSigmoid($result);
# 
#   return ActivateSigmoid($result);
#   
# }

sub isnan { ! defined( $_[0] <=> 9**9**9 ) }


sub AdjustResult
{
my $x = shift;
#         if(isinf($x) == 1) {
#                 return 1.0;
#         }
#         if(isinf($x) == -1) {
#                 return -1.0;
#         }
        if(isnan($x)) {
                return 0.0;
        }

        return $x;
}


sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub divide {
  my ($x, $y) = @_;
  if($y==0){
    return 1.0;
  }
  else{
    return $x/$y;
  }
}
# log ( Concat ( 0 . 8 3 ) ) / Concat ( . 4 4 3 )
#log ( tan ( tan ( Concat ( . 7 0 3 ) ) ) ) * G426
#cosin ( tan ( ( ( G446 - G307 + G426 / G257 ) - G187 + G172 - cosin ( G15 ) ) ) )
#( ( G447 ^ G330 ) + tan ( G73 ) ) - cosin ( G96 ) * G463 ^ G424 ^ ( G155 ^ G346 / sin ( G127 ) - ( log ( G168 ) - sin ( G396 ) ) ) + ( ( ( G78 * G466 ) - ( G260 ^ G349 ) ) + ( G359 * G245 ) / G212 ^ G74 * cosin ( log ( G483 / G412 ) ) )
#( G258 ^ Concat ( . 0 5 3 ) )


# "PA","(","W","(","(","Concat","(","7",".","2","1","4",")",
#       "/","(","Concat","(","0","1",")","-","Concat","(","0","1",")",")",")",",",
#       "G10",")",",","W","(","Concat","(","1","9",".","6","7","5",")",",","G19",
#       ")",",","2",")"
# PA((7.21 / (0-0) * G10), 19.67 * G19)
# PA (7.21/0*G10, 19.67*G19)
# PA (1*G10, 19.67*G19;
sub evaluate_ind{
  my @ind_array = @_;
#my $firstnode = 1 * $ind_array[10];
my $firstnode = 0.165306 * $ind_array[10];
my $secondnode = 19.67 * $ind_array[19];
#   my $thirdnode = -57.43 * $ind_array[1];


# my $result = ( $ind_array[258] ^ .05  );
# print "value is $ind_array[258]\n";
# my $result = $ind_array[258] ** .05;
# print "result=$result\n";
my $result = $firstnode + $secondnode;
#   if($secondnode == 0){
#     $result = 10000000;
#   }
#   else{
#     $result = $firstnode / $secondnode;
#   }
  
 #  my $result = $firstnode / $secondnode;
  
# print " $ind_array[5] $ind_array[17] $ind_array[19] => ", ActivateSigmoid($result);

  return ActivateSigmoid($result);
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
# print "ind status=$$array[0] ";
  my $sc = evaluate_ind(@$array);
  $traintotal++;
  my @ind = @$array;
  my $err = $ind[0] - $sc;
print "total=$square_error_total => ind status=$ind[0] and predicted status = $sc\n";  
  
  $square_error_total += $err * $err;
# print " $err $square_error_total\n";
}

print "square_error_total=$square_error_total\n";
print "training mean square error: ",$square_error_total/$traintotal, "\n";
print "training total inds=$traintotal\n";

my @data = read_data($ARGV[1]); # read training data
my $testtotal = 0;
my $testerror = 0;
$square_error_total = 0;
foreach my $array(@data){
# print "ind status=$$array[0] ";
  my $sc = evaluate_ind(@$array);
  my @ind = @$array;
# print "score is $sc ";
# print "status for $ind[$#ind] is ", ActivateSigmoid($ind[$#ind]), "\n";
  $testtotal++;
  my $err = $ind[0] - $sc;
# printf "adding %g\n", $err*$err;
  $square_error_total += $err * $err;
# print " $err $square_error_total\n";
}
print "square_error_total=$square_error_total\n";
print "testing mean square error:  ",$square_error_total/$testtotal, "\n";
print "testing total inds=$testtotal\n";

#for(my $i=-100; $i<50; $i+=10){
#  printf "$i =>%e\n", ActivateSigmoid($i);
#}


