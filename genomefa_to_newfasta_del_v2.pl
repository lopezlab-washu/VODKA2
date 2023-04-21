#!/usr/bin/env perl

use warnings;
use strict;
my $USAGE = "perl genomefa_to_newfasta_del_v2.pl <input fasta> <bases_from_the_right> <size> <outputfile> <gap_size>

<input fasta>
<bases_from_the_right>: e.g. 3000
<size>: e.g. 100
<outputfile>
<gap_size>: e.g. 10

";
if (@ARGV<4){
    die $USAGE;
}
my $fasta = $ARGV[0];
my $bases_to_use = $ARGV[1];
$bases_to_use *= -1;
my $bp = $ARGV[2];
my $out = $ARGV[3];
my $gap = $ARGV[4];
my $ft = `file $fasta`;
chomp($ft);
my $gz = "";
if ($ft =~ /compressed/){
    $gz = "true";
}
else{
    $gz = "false";
}
if ($gz eq "true"){
    my $pipecmd = "gunzip -c $fasta";
    open(IN, "-|", $pipecmd) or die "Opening pipe [$pipecmd]: $!\n+";
}
else{
    open(IN, $fasta) or die "cannot open file '$fasta'\n";
}

my $onelineseq = "";
while(my $line = <IN>) {
    if($line =~ />/) {
	next;
    }
    else {
	chomp($line);
	$line = uc $line;
	$onelineseq .= $line;
    }
}
close(IN);
my $seql = length($onelineseq);

my $for_gc = $seql+$bases_to_use;
my $threekb = substr($onelineseq, $bases_to_use);

my $length = length($threekb);

my $s2_cnt = 0;
my $s2_start = $length - 15;

open(OUT, ">$out");
open(OT, ">$out.temp");

while ($s2_start > 0){
    my $s1_cnt = 0;
    my $s2 = substr($threekb, -15-$s2_cnt, $bp);
    for (my $i=0;$i<$bp-15;$i++){
        my $s1 = substr($threekb, 0, $s1_cnt+15);
        my $s1_end = $s1_cnt+15;
        if ($s2_start - $s1_end >= $gap){ #deletion size >= gap
            my $B = $s1_end; #A
            my $R = $s2_start + 1; #C
            my $GB = $B + $for_gc; #A
            my $GR = $R + $for_gc; #C
            print OUT ">$GB"."_"."$GR"."("."$B"."_"."$R)\n";
            print OUT "$s1$s2\n";
            print OT ">$GB"."_"."$GR"."("."$B"."_"."$R)\n";
            print OT "$s1|||$s2\n";
        }
        else{
            last;
        }
        $s1_cnt++;
    }
    $s1_cnt = 0;
    for (my $i=0;$i<$length;$i++){
	my $s1 = substr($threekb, $s1_cnt, $bp);
	my $s1_end = $s1_cnt+$bp;
	if ($s2_start - $s1_end >= $gap){ #deletion size >= gap
	    my $B = $s1_end; #A
	    my $R = $s2_start + 1; #C
	    my $GB = $B + $for_gc; #A
	    my $GR = $R + $for_gc; #C
	    print OUT ">$GB"."_"."$GR"."("."$B"."_"."$R)\n";
	    print OUT "$s1$s2\n";
	    print OT ">$GB"."_"."$GR"."("."$B"."_"."$R)\n";
	    print OT "$s1|||$s2\n";
	}
	else{
	    last;
	}
	$s1_cnt++;
    }
    $s2_cnt++;
    $s2_start = $length - 15 - $s2_cnt;
}
close(OUT);
close(OT);

print "got here\n";

