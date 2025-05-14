#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;

sub usage{
	print STDERR << "USAGE";
Name

	blast_parser.pl -- parse the BLAST result and convert to tabular format.

Description

	The BLAST result file including many useful information but not intuitionistic and hard to process by program.
	So this program is written, to get the information and list them in lines on screen or be saved in a file by using ">".
	In the process, this program only keep the result of one query in the memory once. So the memory consume is very small.
	The same as other programs, it also gives some parameters, so you can filter the dissatisfactory alignments easily.
    
	The output format is universal for all the blast formats (include blastn, blastp, tblastn, blastx, and tblastx).
	The fields are seperated by "\t" in each line, If the value of a field is empty, we represent it with "--".
	The order number and description tag are listed below, the meanings of these tags are the same from raw blast result.
	1:Query_id  2:Query_length  3:Query_start  4:Query_end  5:Subject_id  6:Subject_length  7:Subject_start  
	8:Subject_end  9:Identity  10:Positive  11:Gap  12:Align_length  13:Score  14:E_value  15:Query_annotation  16:Subject_annotation
	

Version
	
	Author: fanwei    (fanwei\@genomics.org.cn)
	Author: sunjuan	  (sunjuan\@genomics.org.cn)
	Author: zhushilin (zhushilin\@pubbiotech.com)
	Version: 4.0	Date: 2008-12-10
	
Usage

  	perl blast_parser.pl [options] input_file
	-nohead     do not show the first instruction line.
	-tophit     integer, to set how many subjects for a query to be displayed. 
	-topmatch   integer, to set suits(results of one subject match one query) to be displayed. 
	-eval       float or exponent,to filter the alignments which worse than the E-value cutoff.
	-verbose    output verbose information to screen.
	-help       output help information to screen.

Exmple

	1. Run with the default parameters, this will output all the alignments
	perl blast_parser.pl test_chr_123.seq.bgf.pep.1000.10.blast > test_chr_123.seq.bgf.pep.1000.10.blast.tab
	
	2. Run with user specified Parameters:
	perl blast_parser.pl -tophit 2 test_chr_123.seq.bgf.pep.1000.10.blast > test_chr_123.seq.bgf.pep.1000.10.blast.tab	
	perl blast_parser.pl -topmatch 2 test_chr_123.seq.bgf.pep.1000.10.blast > test_chr_123.seq.bgf.pep.1000.10.blast.tab	
	perl blast_parser.pl -tophit 3 -topmatch 2 -eval 1e-5 test_chr_123.seq.bgf.pep.1000.10.blast > test_chr_123.seq.bgf.pep.1000.10.blast.tab
	perl blast_parser.pl -nohead -tophit 3 -topmatch 2 -eval 1e-5 test_chr_123.seq.bgf.pep.1000.10.blast > test_chr_123.seq.bgf.pep.1000.10.blast.tab

USAGE
	exit 0;
}

my ($Nohead,$Tophit,$Topmatch,$Eval);
my ($Verbose,$Help);
GetOptions(
	"nohead"=>\$Nohead,
	"tophit:i"=>\$Tophit,
	"topmatch:i"=>\$Topmatch,
	"eval:f"=>\$Eval,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
usage() if (@ARGV==0 || $Help);

my $blast_file = shift;


##convert blast raw result to tabular format
&parse_blast($blast_file,$Tophit,$Topmatch,$Eval,$Nohead);


####################################################
################### Sub Routines ###################
####################################################


##parse the BLAST files, and output in tabular formats 
####################################################
sub parse_blast{
	my ($file,$tophit,$topmatch,$eval,$nohead) = @_;
	open (BLAST,"$file") || die ("Could not open the blast file.");

	print "Query_id\tQuery_length\tQuery_start\tQuery_end\tSubject_id\tSubject_length\tSubject_start\tSubject_end\t",
		"Identity\tPositive\tGap\tAlign_length\tScore\tE_value\tQuery_annotation\tSubject_annotation\n" unless(defined $nohead);

	$/="Query= ";
	<BLAST>;

	my $database;
	while (<BLAST>) {
		chomp;
		next if(/\* No hits found \*/);
		my @cycle=split (/\n>/,$_);

		my ($pointer,$query,$query_len,$subject,$subject_len,$query_annotation,$subject_annotation);
		if ($cycle[0]=~/(\S+)\s+(.*?)\n\s*Length=([\d\,\.]+)/s) {
			$query=$1;
			$query_annotation=$2;
			$query_len=$3;
			$query_len =~ s/,//g;
			$query_annotation =~ s/\s+/ /g;
			$query_annotation="--" if(!$query_annotation || $query_annotation eq " ");
		}


		shift @cycle;
		for (my $i=0; $i<@cycle; $i++) {
			last if((defined $tophit) && $i>$tophit-1);
			if ($cycle[$i]=~/(\S+)\s+(.*?)\s*Length=([\d\,\.]+)/s) {
				$subject=$1;
				$subject_annotation=$2;
				$subject_len=$3;
				$subject_len =~ s/,//g;
				$subject_annotation=~ s/\s+/ /g;
				$subject_annotation="--" if(!$subject_annotation || $subject_annotation eq " ");
			}

			my @cycle_inner=split (/Score =/,$cycle[$i]);
			shift @cycle_inner;
			for (my $j=0; $j<@cycle_inner; $j++) {
				last if((defined $topmatch) && $j>$topmatch-1);
				$pointer->[$i][$j]{score}=$1 if($cycle_inner[$j]=~/([\d\,\.]+)\s*bits/);
				$pointer->[$i][$j]{e_value}=$1 if($cycle_inner[$j]=~/Expect[^=]*=\s*(\S+)/); 
				$pointer->[$i][$j]{e_value} =~ s/,$//;
				$pointer->[$i][$j]{e_value}=~s/^e/1e/;
				if ($cycle_inner[$j]=~/Identities\s*=\s*[\d\,\.]+\/([\d\,\.]+)\s*\((\S+)\%\)/s) {
					$pointer->[$i][$j]{align_len}=$1;
					$pointer->[$i][$j]{identity}=$2/100;
				}
				last if((defined $eval) && $pointer->[$i][$j]{e_value}>$eval);

				$pointer->[$i][$j]{positive}=($cycle_inner[$j]=~/Positives\s*=\s*\S+\s*\((\S+)\%\)/s)? $1/100 : "--";

				$pointer->[$i][$j]{gap}=($cycle_inner[$j]=~/Gaps\s*=\s*\S+\s*\((\S+)\%\)/)? $1/100 : 0;

				$pointer->[$i][$j]{q_start}=$1 if($cycle_inner[$j]=~/Query\s*([\d\,\.]+)\s*/);
				$pointer->[$i][$j]{s_start}=$1 if($cycle_inner[$j]=~/Sbjct\s*([\d\,\.]+)\s*/);
				$pointer->[$i][$j]{q_start} =~ s/,//g; 
				$pointer->[$i][$j]{s_start} =~ s/,//g; 

				$cycle_inner[$j]=~s/\n\n\n\nLambda\s*.+$//s;
				$cycle_inner[$j]=~s/\s+$//s;
				if ($cycle_inner[$j]=~/Query\s*[\d\,\.]+\s*\D+([\d\,\.]+)\D+?Sbjct\s*[\d\,\.]+\s*\D+([\d\,\.]+)$/s) {
					$pointer->[$i][$j]{q_end}=$1;
					$pointer->[$i][$j]{s_end}=$2;
					$pointer->[$i][$j]{q_end} =~ s/,//g; 
					$pointer->[$i][$j]{s_end} =~ s/,//g; 
				}

				print "$query\t$query_len\t$pointer->[$i][$j]{q_start}\t$pointer->[$i][$j]{q_end}\t",
					"$subject\t$subject_len\t$pointer->[$i][$j]{s_start}\t$pointer->[$i][$j]{s_end}\t",
					"$pointer->[$i][$j]{identity}\t$pointer->[$i][$j]{positive}\t$pointer->[$i][$j]{gap}\t$pointer->[$i][$j]{align_len}\t",
					"$pointer->[$i][$j]{score}\t$pointer->[$i][$j]{e_value}\t$query_annotation\t$subject_annotation\n";
			}
		}
	}
	$/="\n";
	close(BLAST);
}
