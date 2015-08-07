#!/usr/local/bin/perl -w
use strict;
use warnings;
use feature "switch";
#no warnings "experimental::smartmatch";
use File::Basename;

system("clear");

if (scalar @ARGV<4) {
	formatting();
}

my ($species, $mirna, $dbdir, $resultdir, $mode) = @ARGV;
# $dbdir = "where the database locates";
my ($resultfile, $logfile, $temp, $command, $mirnaid, $job, $flag, $len);
$job = fileparse($resultdir);

$resultfile = "$resultdir/results.txt";
$logfile = "$resultdir/log.txt";

checkArguments();
if ($mode == 0) {
	alignment();
}
else {
	predict();
}


sub formatting {
	print "	Description: Predict for non-canonical miRNA targets using the ncMirTar algorithm.

	USAGE:
		perl $0 species miRNA_seq\n";
	exit 1;
}

sub checkArguments {
	open(LOG, '>', $logfile) || die "'$logfile' $!";
	open(OUTPUT, '>', $resultfile) or do {
		print LOG "Could not open file '$resultfile'\n";
		exit 1;
	};

	$species =~ tr/A-Z/a-z/;	# capital to small
	$temp = "$dbdir/$species";
	if(! -d $temp) {
		print OUTPUT "$species : species does not exist\n";
		exit 1;
	}
	$mirna =~ tr/a-z/A-Z/;
	if (! $mirna =~ /[ACGTU]/) {
		print OUTPUT "$mirna has invalid nucleotide\n";
		exit 1;
	}
	if (length($mirna)<15 || length($mirna)>30) {
		print OUTPUT "$mirna miRNA size should be between 15 and 30\n";
		exit 1;
	}
	$mirna =~ tr/T/U/;
}

sub alignment {
	my ($score,$evalue,$size,$percent);
	$flag = 0;
	$len = length($mirna);
	$temp = "$resultdir/blast_querry";
	open(FILE, '>', $temp) or do {
		print LOG "$temp fails to write\n";
		exit 1;
	};
	if ($species eq "human") {
		print FILE ">hsa-";
		$command = "$dbdir/script/./blastn -task blastn-short -evalue 0.001 -strand plus -db $dbdir/script/blastdb/human_miR -query $temp -out $resultdir/blast_result";
	}
	else {
		print FILE ">mmu-";
		$command = "$dbdir/script/./blastn -task blastn-short -evalue 0.001 -strand plus -db $dbdir/script/blastdb/mouse_miR -query $temp -out $resultdir/blast_result";
	}
	print FILE "$job\n$mirna";
	close FILE;

	system($command) ==0 or do {
		print LOG "$command does not complete\n";
		exit 1;
	};
	
	$temp = "$resultdir/blast_result";
	open(FILE, '<', $temp) or do {
		print LOG "$temp fails to open\n";
		exit 1;
	};
	while ($temp=<FILE>) {
		$temp =~ s/\r|\n//g;
		given ($flag) {
			when(0) {
				if (substr($temp,0,1) eq '>') {
					print OUTPUT "$temp, ";
					$mirnaid = substr($temp,2);
					$flag = 1;
				}
			}
			when(1) {
				print OUTPUT "$temp, ";
				$flag = 2;
			}
			when(2) {
				$flag = 3;
			}
			when(3) {
				my @number1 = $temp =~ /(\d+\.?[eE]?-?\d+)/g;	# of matches
				$size = $number1[1];
				$temp = <FILE>;
				my @number2 = $temp =~ /(\d+%)/g;	# of identities
				$percent = $number2[0];
				if (($size==$len)&&($percent eq "100%")) {last;}
				$score = $number1[0];
				$evalue = $number1[2];
				print OUTPUT "$score, $evalue\n";
				$flag = 4;
			}
			when(4) {
				$flag = 5;	# Strand
			}
			when(5) {
				$flag = 6;
			}
			when(6) {
				if ($temp eq '') {$flag=7;}
				else {print OUTPUT "$temp\n";}
			}
			when(7) {
				if (substr($temp,0,1) eq '>') {
					print OUTPUT "$temp, ";
					$mirnaid = substr($temp,2);
					$flag = 1;
				}
			}
		}
	}
	close FILE;
	close OUTPUT;
	close LOG;
	print "$flag\n";
	if ($flag==0) {
		predict();
	}
	elsif ($flag==3) {
		system("perl $dbdir/script/search_ncMirTar.pl $mirnaid NA $dbdir $resultdir");
	}
	system("rm -rf $resultdir/blast_*");
}

sub predict {
	$temp = "$resultdir/blast_querry";
	open(FILE, '>', $temp) or do {
		print LOG "$temp fails to write\n";
		exit 1;
	};
	if ($species eq "human") {
		print FILE ">hsa-";
	}
	else {
		print FILE ">mmu-";
	}
	print FILE "$job\n$mirna";
	close FILE;
	
	$command = "$dbdir/script/./screen_ncMirTar.o $resultdir/blast_querry $dbdir";
	system($command);
	$command = "$dbdir/script/./miR_seed.o $resultdir/blast_querry";
	system($command);
	system("rm -rf $resultdir/screening.fasta");
	$command = "$dbdir/script/./model6_ncMirTar.o $resultdir/blast_querry $dbdir";
	system($command);
	system("rm -rf $resultdir/miRNA_mir");
	$command = "$dbdir/script/./model5_ncMirTar.o $resultdir/blast_querry $dbdir";
	system($command);
	system("rm -rf $resultdir/libsvm*");
	$command = "$dbdir/script/./rank_ncMirTar.o $resultdir/blast_querry $dbdir";
	system($command);
	system("rm -rf $resultdir/blast_*");
	# $temp = "$resultdir/blast_result";
	# do {print ".";sleep(1);}
	# while (-f $temp);
}