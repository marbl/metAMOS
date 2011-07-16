#!/usr/bin/perl

open (IN1, $ARGV[0]);
open (IN2, $ARGV[1]);

my @list1 = ();
my @list2 = ();

while (<IN1>){

    chomp;
    my $name = $_;
    $name =~ s/^@//;
    $name =~ s/\/\d$//; 
    push @list1, $name;

    <IN1>; <IN1>; <IN1>;
}
while (<IN2>){

    chomp;
    my $name = $_;
    $name =~ s/^@//;
    $name =~ s/\/\d$//; 
    push @list2, $name;

    <IN2>; <IN2>; <IN2>;
}
my $i =0;
my $j = 0;

while ($i <= $#list1 && $j <= $#list2){
    if ($list1[$i] eq $list2[$j]){
	print $list1[$i], "\n";
	$i++; $j++;
    } else {
	# try to resynchronize
	my $k = 1; 
	my $l = 1;

	# advance in list 1
	while ($k + $i <= $#list1 || $l + $j <= $#list2){
	    if ($i + $k <= $#list1 && $list1[$i + $k] eq $list2[$j]){
              #  print "i + k = $i + $k, j = $j\n";
		for (my $m = 0; $m <= $k;  $m++){
		    print $list1[$m+$i], "\n";
		}
		$i = $i + $k + 1; $j++;
               # print "i = $i, j = $j\n";
		goto END;
	    } else {
		$k++;
	    }
	    if ($j + $l <= $#list2 && $list1[$i] eq $list2[$j + $l]){
                #print "i = $i, j + l = $j + $l\n";
		for (my $m = 0; $m <= $l; $m++){
		    print $list2[$m + $j], "\n";
		}
		$i++; $j = $j + $l + 1;
                #print "i = $i, j = $j\n";
		goto END;
	    } else {
		$l++;
	    }
	   # print $list1[$i++], "\n";
	   # print $list2[$j++], "\n";
	} # while i+k or j+l not at the end
	    print $list1[$i++], "\n";
	    print $list2[$j++], "\n";
      END:
    }
} # while both not at the end
while ($i <= $#list1){print $list1[$i++], "\n";}
while ($j <= $#list2){print $list2[$j++], "\n";}
