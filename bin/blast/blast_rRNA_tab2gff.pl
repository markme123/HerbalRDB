#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use Cwd qw/cwd abs_path/;

sub usage{
	print STDERR << "USAGE";
Description

change blast table to gff

Version

  Author: Dad
  Note:

Usage
  
  blast_rRNA_tab2gff.pl [options] <rRNA.blast.tab> <prefix>
  --rRNA_tab  <str>    rRNA.blast.tab
  --prefix    <str>    set a prefix name for the gene ID in gff3 result file  
  --help               output help information to screen  

Exmple

blast_rRNA_tab2gff.pl --rRNA_tab chr10_11.fa.1.rRNA.blast.tab --prefix chr10_11.fa.1

USAGE
	exit 0;
}

##get options from command line into variables and set default values
my ( $RRNA_tab, $Prefix, $Help );
GetOptions(
	"rRNA_tab:s" => \$RRNA_tab,
	"prefix:s"   => \$Prefix,
	"help"       => \$Help
);

usage() if ( $Help );

####################################################
################### Sub Routines ###################
####################################################

###1:Query_id  2:Query_length  3:Query_start  4:Query_end   5:Subject_id  6:Subject_length  7:Subject_start  8:Subject_end
###9:Identity  10:Positive  11:Gap  12:Align_length  13:Score  14:E_value  15:Query_annotation  16:Subject_annotation
####################################################
sub tab_to_gff3 {
	my $file    = shift;
	my $pre_tag = shift;
	my $pre_tagf = shift;
	my $output;

	#$pre_tag .= "_" if ($pre_tag);
	$pre_tag = ($pre_tag) ? "$pre_tag\_" : "";
	my $mark = "0001";
	open IN, $file || die "fail $file";
	while (<IN>) {
		my @t = split /\t/;
		next if ( ( $t[3] - $t[2] ) / $t[1] < 0.2 );    #### Add by Zhu Shilin 2015-03-12
		my ( $target_id, $rRNA_type ) = ( $1, $2 ) if ( $t[0] =~ /(\w+)\#([\w\.]+)/ );
		my $gene_id = $pre_tag . $rRNA_type . "_" . $mark;
		my $seq_id  = $t[4];
		my $score   = $t[13];
		my $strand  = ( $t[6] <= $t[7] ) ? "+" : "-";
		my ( $gene_start, $gene_end ) = ( $t[6] <= $t[7] ) ? ( $t[6], $t[7] ) : ( $t[7], $t[6] );
		$output .=
		  "$seq_id\tblastn\trRNA\t$gene_start\t$gene_end\t$score\t$strand\t.\tID=$gene_id;Target=$target_id $t[2] $t[3];Annotation=\"$t[14]\"\n";
		$mark++;
	}
	close IN;

	open OUT, ">$pre_tagf.rRNA_blast.gff" || die "fail creat $file";
	print OUT "##gff-version 3\n$output";
	close OUT;

}

tab_to_gff3( "$RRNA_tab", "", "$Prefix" );