#!/usr/local/bin/perl -w
use strict;
use warnings;

system("clear");

if (scalar @ARGV!=4) {
	formatting();
}

my ($mirna, $mrna, $dbdir, $resultdir) = @ARGV;
# $dbdir = "where the database locates";
my ($TSfile, $resultfile, $logfile, $temp, $path, $path_TS, $species, $mirnaname, $refseq, $flag);
$resultfile = "$resultdir/results.txt";
$TSfile = "$resultdir/results_TS.txt";
$logfile = "$resultdir/log.txt";
	
checkArguments();

sub formatting {
	print "	Description: Search for predicted microRNA targets using the ncMirTar algorithm. 

	USAGE:
		perl $0 microRNA_ID mRNA_ID Directory_Dabase Directory_Result\n";
	exit 1;
}

sub checkArguments {
	open(LOG, '>', $logfile) || die "Could not open file '$logfile' $!";
	open(OUTPUT, '>', $resultfile) or do {
		print LOG "Could not open file '$resultfile'\n";
		exit 1;
	};
	open(OUTPUT_TS, '>', $TSfile) or do {
		print LOG "Could not open file '$TSfile'\n";
		exit 1;
	};
	
	if ($mirna ne "NA") {	# if user inputs a miRNA
		substr($mirna,0,3) =~ tr/A-Z/a-z/;	# capital to small
		if (index($mirna,"hsa",0)==0) {
			$species = "human";
		}
		elsif (index($mirna,"mmu",0)==0) {
			$species = "mouse";
		}
		else {
			print OUTPUT "$mirna : species does not exist\n";
			exit 1;
		}
		
		substr($mirna,4,3) =~ tr/LETMIr/letmiR/;
		$mirnaname = substr($mirna,4);
		$path = "$dbdir/$species/$mirnaname";
		$path_TS = "$dbdir/TargetScan/$species/$mirnaname";
		
		if(! -d $path) {
			print OUTPUT "$mirna : microRNA does not exist\n";
			exit 1;
		}
	}

	if ($mrna ne "NA") {
		if (index($mrna,"_",2)==2) {	# if user inputs a RefSeq Accession ID
			substr($mrna,0,2) =~ tr/nm/NM/;
			if (index($mrna,"NM",0)==0) {
				$refseq = $mrna;
				if ($mirna ne "NA") {
					$temp = join('',"$dbdir/gene/RefSeq_", substr($mirna,0,3), '/',$refseq);
				}
				else {
					$temp = "$dbdir/gene/RefSeq_hsa/$refseq";
					$species = "human";
					if (! -f $temp) {
						$temp = "$dbdir/gene/RefSeq_mmu/$refseq";
						$species = "mouse";
					}
				}
				open(FILE, '<', $temp) || do {
					print OUTPUT "$refseq transcript does not exist\n";
					exit 1;
				};
				if ($mrna=<FILE>) {
					chomp $mrna;
					my @strings = $mrna =~ /\((\w+)\)/x;	# extract gene in the bracket
					$mrna = $strings[0];
				}
				else {
					print LOG "$refseq file is broken\n";
					exit 1;
				}
				close FILE;
			}
		}
		else {	# if user inputs a gene name
			if ($mirna ne "NA") {
				$flag = 3;
				$temp = join('',"$dbdir/gene/gene_", substr($mirna,0,3), '/',substr($mrna,0,1));
			}
			else {
				$flag = 1;
				$temp = join('',"$dbdir/gene/gene_hsa/",substr($mrna,0,1));
				$species = "human";
			}
			if (open(FILE, '<', $temp)) {
				while ($temp=<FILE>) {
					chomp $temp;
					$temp =~ s/\r//g;
					if ($temp eq $mrna) {
						$flag = $flag-1;
						last;
					}				
				}
				close FILE;
			}
			if ($flag==3) {
				print OUTPUT "$mrna gene does not exist\n";
				exit 1;
			}
			else {
				if ($flag==1) {
					$temp = join('',"$dbdir/gene/gene_mmu/",substr($mrna,0,1));
					$species = "mouse";
					open(FILE, '<', $temp) || do {
						print OUTPUT "$mrna gene does not exist\n";
						exit 1;
					};
					while ($temp=<FILE>) {
						chomp $temp;
						if ($temp eq $mrna) {
							$flag = $flag-1;
							last;
						}
					}
					close FILE;
					if ($flag==1) {
						print OUTPUT "$mrna gene does not exist\n";
						exit 1;
					}
				}
			}
		}
	}
	
	if ($mirna ne "NA") {
		if ($mrna ne "NA") {
			$temp = "$path/$mrna";
			if(! -e $temp) {
				print OUTPUT "No non-canonical target is found\n";
				if ($refseq) {
					miRNA_refseq_TS();
				}
				else {
					miRNA_gene_TS();
				}
			}
			else {
				if ($refseq) {
					miRNA_refseq();
				}
				else {
					miRNA_gene();
					print OUTPUT_TS "No canonical target is found\n";
				}
			}
			
		}
		else {
			miRNA_all();
			miRNA_all_TS();
		}
	}
	else {
		if ($refseq) {
			refseq_all();
			refseq_all_TS();
		}
		else {
			gene_all();
			gene_all_TS();
		}
	}
	close OUTPUT;
	close OUTPUT_TS;
	close LOG;
}

sub miRNA_refseq {
#	print "One miRNA and one RefSeq:\n";
	$temp = "$path/$mrna";
	open(FILE, '<', $temp) || do {
		print OUTPUT "No non-canonical target is found\n";
		exit 0;
	};
	my $count = 0;
	while ($temp=<FILE>) {
		chomp $temp;
		if (index($temp,"@")==0) {
			$flag = index($temp,',');
			if (substr($temp,1,$flag-1) eq $refseq) {
				$temp = substr($temp,$flag);
				substr($temp,2,1) = '-';	# target type
				substr($temp,5,2) = '';	# count
				print OUTPUT "@";
				print OUTPUT $mirna;
				print OUTPUT "::$mrna,$refseq$temp\n";
				$flag = 1;
				$count = 1;
			}
			else {
				$flag = 0;
			}
		}
		else {
			if ($flag>0) {
				print OUTPUT "$temp\n";
			}
		}
	}
	close FILE;
	if ($count==0) {
		print OUTPUT "No non-canonical target is found\n";
		miRNA_refseq_TS();
	}
	else {
		print OUTPUT_TS "No canonical target is found\n";
	}
}

sub miRNA_refseq_TS {
#	print "One miRNA and one RefSeq from TargetScan:\n";
	open(FILE, '<', $path_TS) || do {
		print OUTPUT_TS "No canonical target is found\n";
		exit 0;
	};
	my @array;
	$flag = 0;
	while ($temp=<FILE>) {
		chomp $temp;	# avoid \n at the end
		@array = split(',', $temp);
		if ($array[1] eq $refseq) {
			print OUTPUT_TS "$mirna\t$mrna\t$refseq\n";
			$flag = 1;
			last;
		}
	}
	close FILE;

	if ($flag==0) {
		print OUTPUT_TS "No canonical target is found\n";
	}
}

sub miRNA_gene {
#	print "One miRNA and one gene:\n";
	$temp = "$path/$mrna";
	open(FILE, '<', $temp) || do {
		print OUTPUT "No non-canonical target is found\n";
		exit 0;
	};
	while ($temp=<FILE>) {
		chomp $temp;
		if (index($temp,"@",0)==0) {
			$temp = substr($temp,1);
			$flag = index($temp,',');
			substr($temp,$flag+2,1) = '-';	# target type
			substr($temp,$flag+5,2) = '';	# count
			print OUTPUT "@";
			print OUTPUT $mirna;
			print OUTPUT "::$mrna,$temp\n";
		}
		else {
			print OUTPUT "$temp\n";
		}
	}
	close FILE;
}

sub miRNA_gene_TS {
#	print "One miRNA and one gene from TargetScan:\n";
	open(FILE, '<', $path_TS) || do {
		print OUTPUT_TS "No canonical target is found\n";
		exit 0;
	};
	my @array;
	$flag = 0;
	while ($temp=<FILE>) {
		chomp $temp;	# avoid \n at the end
		@array = split(',', $temp);
		if ($array[0] eq $mrna) {
			print OUTPUT_TS "$mirna\t$mrna\n";
			$flag = 1;
			last;
		}
	}
	close FILE;

	if ($flag==0) {
		print OUTPUT_TS "No canonical target is found\n";
	}
}

sub miRNA_all {
#	print "One miRNA and all gene:\n";
	foreach my $mrnainitial ('A' .. 'Z') {
		$temp = join('',"$dbdir/gene/gene_", substr($mirna,0,3), '/', $mrnainitial);
		if (open(GENELIST, '<', $temp)) {
			while ($mrna=<GENELIST>) {
				chomp $mrna;
				$mrna =~ s/\r//g;
				$temp = "$path/$mrna";
				open(FILE, '<', $temp) || next;
				while ($temp=<FILE>) {
					chomp $temp;
					if (index($temp,"@",0)==0) {
						$temp = substr($temp,1);
						$flag = index($temp,',');
						substr($temp,$flag+2,1) = '-';	# target type
						substr($temp,$flag+5,2) = '';	# count
						print OUTPUT "@";
						print OUTPUT $mirna;
						print OUTPUT "::$mrna,$temp\n";
					}
					else {
						print OUTPUT "$temp\n";
					}
				}
				close FILE;				
			}
			close GENELIST;
		}
		else {
			print LOG "$temp is not open\n";
		}
	}
}

sub miRNA_all_TS {
#	print "One miRNA and all genes from TargetScan:\n";
	open(FILE, '<', $path_TS) || do {
		print OUTPUT_TS "No canonical target is found\n";
		exit 0;
	};
	$mrna = "";
	my @array;
	while ($temp=<FILE>) {
		chomp $temp;	# avoid \n at the end
		@array = split(',', $temp);
		
		if ($array[0] ne $mrna) {
			$mrna = $array[0];
			print OUTPUT_TS "$mirna\t$mrna\n";
		}
	}
	close FILE;
}
			
sub refseq_all {
#	print "All miRNA and one RefSeq:\n";
	my $count=0;
	$temp = join('',"$dbdir/gene/",$species,"_miRNA.fasta");
	open(MIRLIST, '<', $temp) || do {
		print LOG "$temp is not open\n";
		exit 1;
	};
	while ($mirna=<MIRLIST>) {	# miRNA ID
		$temp=<MIRLIST>;	# sequence
		$mirna =~ s/\r|\n//g;
		$mirna = substr($mirna,1);
		$mirnaname = substr($mirna,4);
		$temp = "$dbdir/$species/$mirnaname/$mrna";
		open(FILE, '<', $temp) || next;
		while ($temp=<FILE>) {
			chomp $temp;
			if (index($temp,"@",0)==0) {
				$flag = index($temp,',');
				if (substr($temp,1,$flag-1) eq $refseq) {
					$temp = substr($temp,$flag);
					substr($temp,2,1) = '-';	# target type
					substr($temp,5,2) = '';	# count
					print OUTPUT "@";
					print OUTPUT $mirna;
					print OUTPUT "::$mrna,$refseq$temp\n";
					$flag = 1;
					$count = 1;
				}
				else {
					$flag = 0;
				}
			}
			else {
				if ($flag>0) {
					print OUTPUT "$temp\n";
				}
			}
		}
		close FILE;
	}
	close MIRLIST;
	if ($count==0) {
		print OUTPUT "No target is found\n";
	}
}

sub refseq_all_TS {
#	print "All miRNA and one RefSeq from TargetScan:\n";
	my $count=0;
	$temp = join('',"$dbdir/gene/",$species,"_miRNA.fasta");
	open(MIRLIST, '<', $temp) || do {
		print LOG "$temp is not open\n";
		exit 1;
	};
	my @array;
	my $mrnabackup = $mrna;
	substr($mrnabackup,0,1) =~ tr/a-z/A-Z/;	# capital to small
	while ($mirna=<MIRLIST>) {	# miRNA ID
		$temp=<MIRLIST>;	# sequence
		$mirna =~ s/\r|\n//g;
		$mirna = substr($mirna,1);
		$mirnaname = substr($mirna,4);
		$temp = "$dbdir/TargetScan/$species/$mirnaname";
		open(FILE, '<', $temp) || next;
		while ($temp=<FILE>) {
			chomp $temp;
			@array = split(',', $temp);
			substr($array[0],0,1) =~ tr/a-z/A-Z/;	# capital to small
			if ($array[0] lt $mrnabackup) {
				next;
			}
			else {
				if ($array[0] gt $mrnabackup) {
					last;
				}
				else {
					if ($array[1] eq $refseq) {
						print OUTPUT_TS "$mirna\t$mrna\t$refseq\n";
						$count = 1;
						last;
					}
				}
			}
		}
		close FILE;
	}
	close MIRLIST;
	if ($count==0) {
		print OUTPUT_TS "No canonical target is found\n";
	}
}

sub gene_all {
#	print "All miRNA and one gene:\n";
	$temp = join('',"$dbdir/gene/",$species,"_miRNA.fasta");
	open(MIRLIST, '<', $temp) || do {
		print LOG "$temp is not open\n";
		exit 1;
	};
	while ($mirna=<MIRLIST>) {
		$temp=<MIRLIST>;	# sequence
		$mirna =~ s/\r|\n//g;
		$mirna = substr($mirna,1);
		$mirnaname = substr($mirna,4);
		$temp = "$dbdir/$species/$mirnaname/$mrna";
		open(FILE, '<', $temp) || next;
		while ($temp=<FILE>) {
			chomp $temp;
			if (index($temp,"@",0)==0) {
				$temp = substr($temp,1);
				$flag = index($temp,',');
				substr($temp,$flag+2,1) = '-';	# target type
				substr($temp,$flag+5,2) = '';	# count
				print OUTPUT "@";
				print OUTPUT $mirna;
				print OUTPUT "::$mrna,$temp\n";
			}
			else {
				print OUTPUT "$temp\n";
			}
		}
		close FILE;
	}
	close MIRLIST;
}

sub gene_all_TS {
#	print "All miRNA and one gene from TargetScan:\n";
	my $count=0;
	$temp = join('',"$dbdir/gene/",$species,"_miRNA.fasta");
	open(MIRLIST, '<', $temp) || do {
		print LOG "$temp is not open\n";
		exit 1;
	};
	my @array;
	my $mrnabackup = $mrna;
	substr($mrnabackup,0,1) =~ tr/a-z/A-Z/;	# capital to small
	while ($mirna=<MIRLIST>) {	# miRNA ID
		$temp=<MIRLIST>;	# sequence
		$mirna =~ s/\r|\n//g;
		$mirna = substr($mirna,1);
		$mirnaname = substr($mirna,4);
		$temp = "$dbdir/TargetScan/$species/$mirnaname";
		open(FILE, '<', $temp) || next;
		while ($temp=<FILE>) {
			chomp $temp;
			@array = split(',', $temp);
			substr($array[0],0,1) =~ tr/a-z/A-Z/;	# capital to small
			if ($array[0] lt $mrnabackup) {
				next;
			}
			else {
				if ($array[0] gt $mrnabackup) {
					last;
				}
				else {
					print OUTPUT_TS "$mirna\t$mrna\n";
					$count = 1;
					last;
				}
			}
		}
		close FILE;
	}
	close MIRLIST;
	if ($count==0) {
		print OUTPUT_TS "No canonical target is found\n";
	}
}
