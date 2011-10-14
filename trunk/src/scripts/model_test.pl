#!/usr/bin/perl

use strict;

if(@ARGV < 1){
  print "model_test.pl data.txt\n";
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

## put model here
#PM ( W ( ( Concat ( 1 . 7 9 4 ) - Concat ( 6 . 0 1 4 ) ) , G2 ) , 
#  W ( Concat ( 7 . 0 3 ) , G5 ) , 2 )

#PM(W((1.794-6.014), G2), W(7.03,G5),2)

#fitness = 0.5

# 1 in GEObjective: size= 124 fitness = 0.4925
# sub evaluate_ind{
# 
#   my @ind_array = @_;
# 
#  my $firstw = 1.794-6.014;
#  my $firstnode = $firstw * $ind_array[2];
#  my $secondw = 7.03;
#  my $secondnode = $secondw * $ind_array[5];
# 
# #  print "score  = "; 
# #  print $firstnode - $secondnode;
# #  print "status = $ind_array[0]\n";
#  
#  my $score = ActivateSigmoid($firstnode - $secondnode);
#  
#   print "score  = $score"; 
# #  print $firstnode - $secondnode;
#  print " status = $ind_array[0] predicted=";
#  my $val = (($score > 0.5)?1:0);
#  print "$val";
#  print "\n";
#  
#  return (($score > 0.5)?1:0);
# }¿


# PS ( W ( ( Concat ( 9 . 0 5 4 ) / Concat ( 5 . 8 9 4 ) ) , G10 ) , 
#   W ( Concat ( 8 . 7 5 4 ) , G7 ) , 
#   W ( Concat ( 1 1 ) , G11 ) , 3 )

# sub evaluate_ind{
# 
#   my @ind_array = @_;
# 
#  my $firstw = 9.05 / 5.89;
#  my $firstnode = $firstw * $ind_array[10];
#  my $secondw = 8.75;
#  my $secondnode = $secondw * $ind_array[7];
#  my $thirdw = 1;
#  my $thirdnode = $thirdw * $ind_array[11];
# 
# #  print "score  = "; 
# #  print $firstnode - $secondnode;
# #  print "status = $ind_array[0]\n";
#  
#  my $score = ActivateSigmoid($firstnode - $secondnode - $thirdnode);
#  
#   print "$firstnode-$secondnode-$thirdnode score  = $score"; 
# #  print $firstnode - $secondnode;
#  print " status = $ind_array[0] predicted=";
#  my $val = (($score > 0.5)?1:0);
#  print "$val";
#  
#  print " genotypes are $ind_array[10] $ind_array[7] $ind_array[11]";
#  print "\n";
#  
#  return (($score > 0.5)?1:0);
# }

# PM( W((1.79-6.01),G2), W(2.45,G11),2)
#PA( W(8.22,G1), W(((0.02/(((((0.02/0.24)/83.49)+((98.57-(0/2))/0.43))-((0*9.2)*(6.9-(3.5-0.3))))+4.2))*3),G23),2)
# sub evaluate_ind{
# 
#   my @ind_array = @_;
# 
#   my $firstw = 1.79-6.01;
#   my $firstnode = $firstw * $ind_array[2];
#   my $secondw = 2.45;
#   my $secondnode = $secondw * $ind_array[11];
# 
# #  print "score  = "; 
# #  print $firstnode - $secondnode;
# #  print "status = $ind_array[0]\n";
#  
#   my $score = ActivateSigmoid($firstnode * $secondnode);
#  
#   print "$firstnode*$secondnode score  = $score"; 
# #  print $firstnode - $secondnode;
#   print " status = $ind_array[0] predicted=";
#   my $val = (($score > 0.5)?1:0);
#   print "$val";
#  
#   print " genotypes are $ind_array[23] $ind_array[1] ";
#   print "\n";
#  
#   return (($score > 0.5)?1:0);
# }

# PD( W((9.72+7.1),G23), W(6.04,G18),2)
# sub evaluate_ind{
# 
#   my @ind_array = @_;
# 
#   my $firstw = 9.72+7.1;
#   my $firstnode = $firstw * $ind_array[23];
#   my $secondw = 6.04;
#   my $secondnode = $secondw * $ind_array[18];
# 
# #  print "score  = ";
# #  print $firstnode - $secondnode;
# #  print "status = $ind_array[0]\n";
#  
#   my $score;
#  
#   if($secondnode == 0){
#     $score = ActivateSigmoid($firstnode/1);
#   }
#   else{
#     $score = ActivateSigmoid($firstnode/$secondnode);
#   }
#  
#   print "$firstnode/$secondnode score  = $score"; 
# #  print $firstnode - $secondnode;
#   print " status = $ind_array[0] predicted=";
#   my $val = (($score > 0.5)?1:0);
#   print "$val";
#  
#   print " genotypes are $ind_array[23] $ind_array[1] ";
#   print "\n";
#  
#   return (($score > 0.5)?1:0);
# }

# PS( W(6.5,G21), W(0.36,G23),2)

# PM( W(0.88,G12), W(((0.1/63.97)*(8.4+((7+((8.3-62.55)/(4.8-(3*24.45))))-(1.7*0.74)))),G2),2)
# PA( W(8.22,G1), W(((0.02/(((((0.02/0.24)/83.49)+((98.57-(0/2))/0.43))-((0*9.2)*(6.9-(3.5-0.3))))+4.2))*3),G23),2)
sub evaluate_ind{

  my @ind_array = @_;

#   my $firstnode = 8.22 * $ind_array[1];
#   my $secondw =((0.02/(((((0.02/0.24)/83.49)+((98.57-(0/2))/0.43))-((0*9.2)*(6.9-(3.5-0.3))))+4.2))*3);
#   my $secondnode = $secondw * $ind_array[23];arejus


#   my $node1 = -0.7484 * $ind_array[2];
#   my $node2 = 4.22 * $ind_array[1];
#   my $node3 = -26.0346 * $ind_array[2];
#   my $node4 = 17.3286 * $ind_array[1];

  my $node1 = 0.40 * $ind_array[2];
  my $node2 = 3.75 * $ind_array[1];
  my $node3 = -37.1872 * $ind_array[2];
  my $node4 = 27.94 * $ind_array[1];  

#  print "score  = ";
#  print $firstnode - $secondnode;
#  print "status = $ind_array[0]\n";
 
  my $score = ActivateSigmoid($node1 + $node2 + $node3 + $node4);

print "G1=$ind_array[1] G2=$ind_array[2]=>$score\n";
 
#   if($secondnode == 0){
#     $score = ActivateSigmoid($firstnode/1);
#   }
#   else{
#     $score = ActivateSigmoid($firstnode/$secondnode);
#   }
 
#   print "$firstnode-$secondnode score  = $score"; 
#  print $firstnode - $secondnode;
#   print " status = $ind_array[0] predicted=";
  my $val = (($score > 0.5)?1:0);
#   print "$val";
 
#   print " genotypes are $ind_array[23] $ind_array[1] ";
#   print "\n";
 
  return (($score > 0.5)?1:0);
}


#PA ( W ( 76.95 , V20 ) , W ( ( ( 04.93 - 58.12 ) * 70.67 ) , V9 ) , 2 )

#sub evaluate_ind{
#  my @ind_array = @_;
#  
#  my $param1 = 85.74;
#  my $param2 = 99.01 + 57.05;
#  my $rawscore;
#  if($ind_array[8] ==0){
#    $rawscore = 1.1;
#  }
#  else{
#    $rawscore = $param1 * $ind_array[18] / $param2 * $ind_array[8];
#  }
#  my $score = ActivateSigmoid($rawscore);
# 
#  return (($score > 0.5)?1:0);
#}

#PD ( W ( 85.74 , V19 ) , W ( ( 99.01 + 57.05 ) , V9 ) , 2 )

#sub evaluate_ind{
#  my @ind_array = @_;

#  my $param1 = .68  *  $ind_array[19];
#  my $param2 = 32.39 * $ind_array[9]; 
#  my $rawscore = $param1 * $param2;
#  my $score = ActivateSigmoid($rawscore);
#  $param1 = ($param1 > 0)?1:0;
#  $param2 = ($param2 > 0)?1:0;


#  my $score = ($param1 != $param2);
  

#  return (($score > 0.5)?1:0);
#}
#PM ( W ( ( 3.3 * 83.11 ) , V20 ) , W ( 10.59 , V10 ) , 2 )
#PXOR ( W ( 99.48 , V20 ) , W ( 42.32 , V9 ) , 2 )


sub ActivateSigmoid{
  my $x = shift;
  
  if($x < -709){ return -1.0;}
  if($x > 709){return 1.0;}
  return(1.0 / (1.0 + exp(-$x)));
}

my @data = read_data($ARGV[0]); # read training data
my $traintotal = 0;
my $trainerror = 0;
my $correct = 0;

my $caseright =0;
my $casewrong = 0;
my $controlright = 0;
my $controlwrong = 0;

foreach my $array(@data){
  # status is first element
  my @ind = @$array; 

 if($ind[21] == 3 || $ind[23] == 3){
    print "skip missing\n";
    next;
  }
   my $sc = evaluate_ind(@$array);
  
 $traintotal++;
  
  if($ind[0]==1){
    if($sc==1){
      $caseright++;
    }
    else{
      $casewrong++;
    }
  }
  if($ind[0]==0){
    if($sc==0){
      $controlright++;
    }
    else{
      $controlwrong++;
    }
  }
  
  
  if($sc == $ind[0]){
    $correct++;
  }
  
  if($sc != $ind[0]){
    $trainerror++;
  }
}

print "caseright=$caseright casewrong=$casewrong controlright=$controlright controlwrong=$controlwrong\n";
print "balanced acc: ",($caseright /($caseright+$casewrong) 
  + $controlright/($controlright+$controlwrong))/2, "\n";
print "training acc: ",$correct/$traintotal, "\n";

# my @data = read_data($ARGV[1]); # read training data
# my $testtotal = 0;
# my $testerror = 0;
# foreach my $array(@data){
#   my $sc = evaluate_ind(@$array);
#   $testtotal++;
#   my @ind = @$array;
#   if($sc != $ind[$#ind]){
#     $testerror++;
#   }
# }
# print "testing error:  ",$testerror/$testtotal, "\n";
