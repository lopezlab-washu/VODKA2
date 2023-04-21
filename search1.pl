use warnings;
use strict;
no warnings ('uninitialized', 'substr');
$|=1;

my $USAGE = "perl search1.pl <samfile> <search1 output> <bp> <readLength>

";

if (@ARGV<4){
    die $USAGE;
}

open(INFILE, $ARGV[0]);
open(OUTFILE, ">$ARGV[1]");

my $cutoff = 15;
my $bp = $ARGV[2];
my $readLength = $ARGV[3];
my (@a, $seq, $A, $B, $cigar, $d1, $d2, $t1, $t2, $x, $y, $str);
while(my $line = <INFILE>) {
    undef @a;
    undef $seq;
    undef $A;
    undef $B;
    undef $cigar;
    undef $d1;
    undef $d2;
    undef $t1;
    undef $t2;
    undef $x;
    undef $y;
    undef $str;
    if($line =~ /^@/) {
	next;
    }
    chomp($line);
    @a = split(/\t/,$line);
    if($a[5] =~ /\*/) {
	next;
    }
    my $refid = $a[2];
    $refid =~ /(\d+)\)$/;
    my $endid = $1;
    my $sizeA = $readLength;
    my $sizeB = $readLength;
    my $check = $bp - $readLength + 1;
    if ($endid > $check){
	$sizeA = $bp - $endid + 1;
    }
    $sizeA++;
    $line =~ /NM:i:(\d+)/;
    my $NM = $1;
    my $N = $a[3];
#    print "N:$N\n";
    $cigar = $a[5];
    my $cig = $cigar;
#    print "$cig\n";
    my (%DP,%IP);
    my $di_ct = 1;
    while ($cig =~ /(\d+)M(\d+)(D|I)(\d+)M/){
	my $str = $1 . "M" . $2 . $3 . $4 . "M";
	my $D_or_I = $3;
	if ($D_or_I eq "D"){
	    my $N = $1+$2+$4;
	    my $dpos = $1;
            my $dcnt = $2;
            $DP{"$di_ct.$dpos"} = "$dcnt";
            my $new_str = $N . "M";
            $cig =~ s/$str/$new_str/;
	}
	if ($D_or_I eq "I"){
            my $N = $1+$4;
            my $Icnt = $2;
            my $Ipos = $1;
            $IP{"$di_ct.$Ipos"} = "$Icnt";
            my $new_str = $N . "M";
            $cig =~ s/$str/$new_str/;
	}
	$di_ct++;
    }
    $cig =~ /[^M\d]*(\d+)M/;
    my $M = $1;
    my $A = $sizeA - $N;
#    print "M:$M\n";
    my $B = $N + $M - $sizeA;
#    print "A\tB\n$A\t$B\n";
    foreach my $dposl (keys %DP){
#	print "dpos:$dposl\n";
	my ($ctr,$dpos) = split(/\./,$dposl);
	my $dcnt = $DP{$dposl};
	if ($ctr == 1){
	    if ($dpos < $A){
		$A -= $dcnt;
	    }
	    else{
		$B -= $dcnt;
	    }
	}
	else{
	    if ($dpos <= $A){
                $A -= $dcnt;
            }
            else{
                $B -= $dcnt;
            }
	}
    }
    foreach my $iposl (keys %IP){
        my ($ctr,$ipos) = split(/\./,$iposl);
	my $icnt = $IP{$iposl};
#	print "iposl:$iposl\n";
#	print "ipos\t$ipos\n";
	if ($ctr == 1){
	    if ($ipos < $A){
		$A += $icnt;
	    }
	    else{
		$B += $icnt;
	    }
	}
	else{
            if ($ipos <= $A){
                $A += $icnt;
            }
            else{
                $B += $icnt;
            }
	}
    }
    $seq = $a[9];
    $cig =~ /^(\d+)(.)(\d+)(.)/;
    $d1 = $1;
    $t1 = $2;
    $d2 = $3;
    $t2 = $4;
    $str = "";
=comment #db
    print "$line\n";
    print "$seq\t$cig\n";
    print "d1:$d1\tt1:$t1\td2:$d2\tt2:$t2\n";
    print "A:$A\nB:$B\n";
=cut
    if($t1 eq "S" && $t2 eq "M") {
	$y = $d1;
	$str = $str . substr($seq,0,$y);
	$str = $str . ">>>";
	$str = $str . substr($seq,$y,$A);
	$str = $str . ":::";
	$x = $y + $A;
	$str = $str . substr($seq,$x,$B);
	$str = $str . "<<<";
	$x = $x + $B;
	$y = length($seq) - $x;
	$str = $str . substr($seq,$x,$y);
    } 
    else {
	$cig =~ /^(\d+)(.)/;
	$d1 = $1;
	$t1 = $2;
	$str = "";
	if($t1 eq "M") {
	    $str = ">>>";
	    $str = $str . substr($seq,0,$A);
	    $str = $str . ":::";
	    $x = $A;
	    $str = $str . substr($seq,$x,$B);
	    $str = $str . "<<<";
	    $x = $x + $B;
	    $y = length($seq) - $x;
	    $str = $str . substr($seq,$x,$y);
	}
    }
    if($A >= $cutoff && $B >= $cutoff) {
	print OUTFILE "$a[0]\t$a[2]\t$a[3]\t$a[5]\t$A\t$B\t$NM\t$str\n";
    }
}
close(INFILE);
close(OUTFILE);


print "got here\n";
