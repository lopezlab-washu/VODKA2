use warnings;
use strict;

my $USAGE = "perl organize.pl <outdir> <indexname> <samplenames>

<outdir>
<indexname>
<samplenames> comma separated list of samplenames

";

if (@ARGV<3){
    die $USAGE;
}

my $DIR = $ARGV[0];
my $indexname = $ARGV[1];
my @names = split(",",$ARGV[2]);

my %DATA;
my $str = "";
unless (-d "$DIR/results/"){
    `mkdir $DIR/results/`;
}
my $reads = "$DIR/results/$indexname" . "_RESULTS.fa";
open(OUTR, ">$reads");
for my $samplename (@names){
    $samplename =~ s/^\s+|\s+$//g;
    my $sfile = "$DIR/alignment/$samplename.$indexname.search_output.txt";
    unless (-e $sfile){
	die "HERE:file $sfile does not exist\n";
    }
    my $x = `cat $sfile`;
    $x =~ s/\n/\t$samplename\n/gs;
    $str = $str . $x;
}
chomp($str);
my @a = split(/\n/,$str);
for(my $i=0; $i<@a; $i++) {
    my @b = split(/\t/,$a[$i]);
    my $size = @b;
    if($b[7] !~ /\S/) {
	print "$str";
	exit();
    }
    my $r1 = $b[7];
    $r1 =~ />>>(.+):::(.+)<<</;
    my $t1 = $1;
    my $t2 = $2;
    print OUTR ">"."$b[0]:::::$b[8]\n$t1$t2\n";
    my $min;
    if ($b[4] >= $b[5]){
	$min = $b[5];
    }
    else{
	$min = $b[4];
    }
    if (exists $DATA{$b[1]}){
	$DATA{$b[1]} = $DATA{$b[1]} . "$b[1]\t$b[0]\t$b[4]\t$b[5]\t$min\t$b[2]\t$b[3]\t$b[6]\t$b[8]\t$b[7]\t$t1$t2\n";
    }
    else{
	$DATA{$b[1]} = "$b[1]\t$b[0]\t$b[4]\t$b[5]\t$min\t$b[2]\t$b[3]\t$b[6]\t$b[8]\t$b[7]\t$t1$t2\n";
    }
}
close(OUTR);
my $filename = "$DIR/results/$indexname" . "_RESULTS.txt";

open(OUTFILE, ">$filename");
print OUTFILE "A_C(A_C)\tREAD_ID\tLEFT\tRIGHT\tMIN\tSTART\tCIGAR\tNM\tSEQ_RUN\tREAD\tMAPPED_ONLY\n";
foreach my $key (keys %DATA) {
    print OUTFILE "$DATA{$key}";
}
close(OUTFILE);

print "got here\n";
