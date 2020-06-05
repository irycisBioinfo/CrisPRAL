#!/usr/bin/perl
use Getopt::Long;
@Options = (
		
		{OPT=>"File=s",	VAR=>\$file,	DESC=>"File introduced"}

);

#Check options and set variables
(@ARGV < 1) && (usage()); 
GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

system("../bin/cd-hit/cd-hit-est -bak 1 -c 1 -s 1 -i $file -d 100 -o ./cluster -T 0 -M 1000 -g 1");

system("sed 's/at +\\///' ./cluster.bak.clstr | sed 's/\\.\\.\\./\\t/' | sed 's/nt, >/\\t/' > ./cluster.bak.tsv");


open(CLS,"./cluster");
@fasta = <CLS>;
close CLS;

open(CLSout, ">./cluster.tsv");

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
