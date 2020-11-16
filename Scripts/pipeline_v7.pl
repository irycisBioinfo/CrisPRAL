#!/usr/bin/perl

use Getopt::Long;
use Cwd 'abs_path';
$PATH = abs_path($0);
$PATH =~ s/\/pipeline_v7.pl//;

@Options = (
		
		{OPT=>"r1=s",	VAR=>\$r1,	DESC=>"Reads 1"},
		{OPT=>"r2=s",	VAR=>\$r2,	DESC=>"Reads 2"},
		{OPT=>"single_end=s",	VAR=>\$single_end,	DESC=>"Boolean - Selects processing for single end or paired end data"},
		{OPT=>"trimA1=s", VAR=>\$trimA1, DEFAULT => "0",DESC=>"R1 Adapter trimming by length"},
		{OPT=>"trimA2=s", VAR=>\$trimA2, DEFAULT => "0",DESC=>"R2 Adapter trimming by length"},
		{OPT=>"Adapter_R1=s", VAR=>\$Adapter_R1, DEFAULT => "",DESC=>"R1 Adapter trimming by sequence"},
		{OPT=>"Adapter_R2=s", VAR=>\$Adapter_R2, DEFAULT => "",DESC=>"R2 Adapter trimming by sequence"},
		{OPT=>"trimP1=s", VAR=>\$trimP1, DEFAULT => "0",DESC=>"1st Primer trimming by length"},
		{OPT=>"trimP2=s", VAR=>\$trimP2, DEFAULT => "0",DESC=>"2nd Primer trimming by length"},
		{OPT=>"F_Primer1=s", VAR=>\$F_Primer1, DEFAULT => "",DESC=>"1st Forward Primer trimming by sequence"},
		{OPT=>"F_Primer2=s", VAR=>\$F_Primer2, DEFAULT => "",DESC=>"2nd Forward Primer trimming by sequence"},
		{OPT=>"R_Primer1=s", VAR=>\$R_Primer1, DEFAULT => "",DESC=>"1st Reverse Primer trimming by sequence"},
		{OPT=>"R_Primer2=s", VAR=>\$R_Primer2, DEFAULT => "",DESC=>"2nd Reverse Primer trimming by sequence"},
		{OPT=>"min_len=s", VAR=>\$minLen, DEFAULT => "200" ,DESC=>"Minimun length for filtering"},
		{OPT=>"Np=s", VAR=>\$Np, DEFAULT => "8" ,DESC=>"N-percent maximum difference"},
		{OPT=>"Nm=s", VAR=>\$Nm, DEFAULT => "6" ,DESC=>"N-minimum overlap (nucleotides)"},
		{OPT=>"cov=s", VAR=>\$cov, DEFAULT =>"0.8" ,DESC=>"Coverage (clustering) (0-1)"},
		{OPT=>"id=s", VAR=>\$id, DEFAULT =>"1" ,DESC=>"Identity (clustering) (0-1)"},
		{OPT=>"primer-error-rate=s", VAR=>\$primererrorrate, DEFAULT =>"0.15" ,DESC=>"cutadapt primer error rate"},
		{OPT=>"tmpdir=s", VAR=>\$tmpdir, DEFAULT =>"." ,DESC=>"Determines tmp dir where data will be stored with an unique id for every user using the app"}
		
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

if ($single_end eq "FALSE"){

	if ($Adapter_R1 eq 'Empty') {

		system("echo -Paired-end Adapter trimming");
		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $r1 -trim_left $trimA1 -out_good $tmpdir/goodA_R1");
		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $r2 -trim_left $trimA2 -out_good $tmpdir/goodA_R2");
		print("\nDone\n");

	} else {
		system("echo -Paired-end Adapter filtering");

		system("cutadapt -j 0 -e $primererrorrate -a $Adapter_R1 -A $Adapter_R2 -n 2 -e 0.2 -o $tmpdir/R1_filteredA.fastq -p $tmpdir/R2_filteredA.fastq $r1 $r2");

		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/R1_filteredA.fastq -trim_left $trimA1 -out_good $tmpdir/goodA_R1");
		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/R2_filteredA.fastq -trim_left $trimA2 -out_good $tmpdir/goodA_R2");
		print("\nDone\n");

	}

	if ($F_Primer1 eq 'Empty') {

		system("echo -Paired-end primer trimming");
		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/goodA_R1.fastq -trim_left $trimP1 -trim_right $trimP2 -out_good $tmpdir/good_R1");
		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/goodA_R2.fastq -trim_left $trimP2 -trim_right $trimP1 -out_good $tmpdir/good_R2");
		print("\nDone\n");

	} else {
		system("echo -Paired-end primer filtering");

		system("cutadapt -j 0 -e $primererrorrate -m 10 -g $F_Primer1 -a $F_Primer2 -G $R_Primer1 -A $R_Primer2 -n 2 -o $tmpdir/R1_filteredP.fastq -p $tmpdir/R2_filteredP.fastq $tmpdir/goodA_R1.fastq $tmpdir/goodA_R2.fastq");

		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/R1_filteredP.fastq -trim_left $trimP1 -trim_right $trimP2 -out_good $tmpdir/good_R1");
		system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/R2_filteredP.fastq -trim_left $trimP2 -trim_right $trimP1 -out_good $tmpdir/good_R2");
		print("\nDone\n");

	}

system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/good_R1.fastq -fastq2 $tmpdir/good_R2.fastq -min_len $minLen -min_qual_mean 30 -out_good $tmpdir/good_filtered -out_bad $tmpdir/bad");

system("$PATH/bin/ea-utils/clipper/fastq-join -p $Np -m $Nm -o $tmpdir/joined $tmpdir/good_filtered_1.fastq $tmpdir/good_filtered_2.fastq ");

system("$PATH/bin/seqtk/seqtk seq -A $tmpdir/joinedjoin > $tmpdir/final.fasta");

}else{#Single-End Section:

if ($Adapter_R1 eq 'Empty') {

	system("echo -Single-end Adapter trimming");
	system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $r1 -trim_left $trimA1 -out_good $tmpdir/goodA_R1");
	print("\nDone\n");

} else {
	system("echo -Single-end Adapter filtering");

	system("$PATH/bin/cutadapt/cutadapt -j 0 -e $primererrorrate -a $Adapter_R1 -n 2 -o $tmpdir/R1_filteredA.fastq $r1");

	system("$PATH/bin/prinseq/prinseq-lite.pl -fastq R1_filteredA.fastq -trim_left $trimA1 -out_good $tmpdir/goodA_R1");
	print("\nDone\n");


}

if ($F_Primer1 eq 'Empty') {

	system("echo -Single-end primer trimming");
	system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/goodA_R1.fastq -trim_left $trimP1 -trim_right $trimP2 -out_good $tmpdir/good_R1");
	print("\nDone\n");

	} else {
	system("echo -Single-end primer filtering");
	system("$PATH/bin/cutadapt/cutadapt -j 0 -e $primererrorrate -m 10 -g $F_Primer1 -a $F_Primer2 -n 2 -o $tmpdir/R1_filteredP.fastq $tmpdir/goodA_R1.fastq");
	system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/R1_filteredP.fastq -trim_left $trimP1 -trim_right $trimP2 -out_good $tmpdir/good_R1");
	print("\nDone\n");

	}
	system("$PATH/bin/prinseq/prinseq-lite.pl -fastq $tmpdir/good_R1.fastq -min_len $minLen -min_qual_mean 30 -out_good $tmpdir/final -out_bad $tmpdir/bad");
	system("$PATH/bin/seqtk/seqtk seq -A $tmpdir/final.fastq > $tmpdir/final.fasta");
}

system("$PATH/bin/cd-hit/cd-hit-est -bak 1 -c $id -s $cov -i $tmpdir/final.fasta -d 100 -o $tmpdir/cluster -T 0 -M 1000 -g 1");

system("sed 's/at +\\///' $tmpdir/cluster.bak.clstr | sed 's/\\.\\.\\./\\t/' | sed 's/nt, >/\\t/' > $tmpdir/cluster.bak.tsv");


open(CLS,"$tmpdir/cluster");
@fasta = <CLS>;
close CLS;

open(CLSout, ">$tmpdir/cluster.tsv");

foreach $l (@fasta)
{
	if($l =~ />/)
	{
		$l =~ s/>//;
		@c = split(/ /,$l);
		$c[0] =~ s/\s//;
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
