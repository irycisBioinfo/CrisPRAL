#!/usr/bin/perl

use Getopt::Long;
use strict;


my @Options;
my $fasta;
my $tmpdir;
my $PATH;
my $prefix;

$PATH = abs_path($0);
$PATH =~ s/\/clustering_clusters.pl//;

#Example options:
@Options = (
		{OPT=>"prefixr=s",	VAR=>\$prefix,	DESC=>"prefix for output files"},
		{OPT=>"fasta=s",	VAR=>\$fasta,	DESC=>"Fasta file"},
		{OPT=>"tmpdir=s", VAR=>\$tmpdir, DEFAULT =>"." ,DESC=>"Determines tmp dir where data will be stored with an unique id for every user using the app"
		
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

system("$PATH/bin/cd-hit/cd-hit-est -bak 1 -c 1 -s 1 -i $tmpdir/$fasta -d 100 -o $tmpdir/$prefix -T 0 -M 1000 -g 1");

system("sed 's/at +\\///' $tmpdir/$prefix.bak.clstr | sed 's/\\.\\.\\./\\t/' | sed 's/nt, >/\\t/' > $tmpdir/$prefix.bak.tsv");

open(CLS,"$tmpdir/$prefix");
@fasta = <CLS>;
close CLS;

open(CLSout, ">$tmpdir/$prefix.tsv");

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
