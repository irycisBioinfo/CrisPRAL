#!/usr/bin/perl
use Getopt::Long;
use Cwd 'abs_path';
$PATH = abs_path($0);
$PATH =~ s/\/pipeline.pl//;

@Options = (
		
		{OPT=>"r1=s",	VAR=>\$r1,	DESC=>"Reads 1"},
		{OPT=>"r2=s",	VAR=>\$r2,	DESC=>"Reads 2"},
		{OPT=>"min_len=s", VAR=>\$minLen, DEFAULT => "200" ,DESC=>"Minimun length for filtering"},
		{OPT=>"Np=s", VAR=>\$Np, DEFAULT => "8" ,DESC=>"N-percent maximum difference"},
		{OPT=>"Nm=s", VAR=>\$Nm, DEFAULT => "6" ,DESC=>"N-minimum overlap (nucleotides)"},
		{OPT=>"cov=s", VAR=>\$cov, DEFAULT =>"0.8" ,DESC=>"Coverage (clustering) (0-1)"},
		{OPT=>"id=s", VAR=>\$id, DEFAULT =>"1" ,DESC=>"Identity (clustering) (0-1)"}
		
	);


#Check options and set variables
(@ARGV < 1) && (usage()); 
GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.file
foreach (@Options) {
	if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	${$_->{VAR}} = $_->{DEFAULT};
	}
}
	




$r1Prinseq = $r1;
$r1Prinseq =~ s/.fastq/good.fastq/;
$r2Prinseq = $r2;
$r2Prinseq =~ s/.fastq/good.fastq/;

$rJoined = $r1;
$rJoined = s/.fastq/.join/;
$rCluster = $r1;
$rCluster = s/.fastq/.cluster/;

system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $r1 -fastq2 $r2 -min_len $minLen -out_good good -out_bad bad");
system("$PATH/bin/ea-utils/clipper/fastq-join -p $Np -m $Nm -o joined good_1.fastq good_2.fastq ");
system("$PATH/bin/seqtk/seqtk seq -A joinedjoin > joined.fasta");
system("$PATH/bin/cd-hit/cd-hit-est -bak 1 -c $id -s $cov -i joined.fasta -d 100 -o cluster -T 20 -M 1000 -g 1");
system("sed 's/at +\\///' cluster.bak.clstr | sed 's/\\.\\.\\./\\t/' | sed 's/nt, >/\\t/' > cluster.bak.tsv");


open(CLS,"cluster");
@fasta = <CLS>;
close CLS;

open(CLSout, ">cluster.tsv");

foreach $l (@fasta)
{
	if($l =~ />/)
	{
		$l =~ s/>//;
		@c = split(/ /,$l);
		print CLSout"$c[0]\t";
	} else {
		print CLSout $l;
	}
}

sub usage {
	foreach (@Options) {
		
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\n\n\n";
	exit(1);
}
