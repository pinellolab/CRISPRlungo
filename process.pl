use strict;
use warnings;

my $settingsFile = shift;

my $usage = "usage: $0 {settingsFile}";
die $usage unless -e $settingsFile;

#read settings file
my %settings = ();
open IN, $settingsFile or die $!;
while (my $line = <IN>)
{
	chomp $line;
	$line =~ s/#.*$//;
	my @lineEls = split "\t", $line;
	my $key = $lineEls[0];
	$key =~ s/^\s*//;
	$key =~ s/\s*$//;
	my $val = $lineEls[1];
	$val =~ s/^\s*//;
	$val =~ s/\s*$//;

	$settings{$key} = $val;
}
close IN;

my $fragmentSize = $settings{'fragment_size'} || 20; 
my $fragmentStepSize = $settings{'fragment_step_size'} || 10; 
my $alignmentExtension = $settings{'alignment_extension'} || 50; #number of bp to extend beyond av read length around cut site for custom index
my @cuts = ();
@cuts = (split " ", $settings{'cut_sites'}) if exists $settings{'cut_sites'};

#may need to merge paired fastqs here
my $fastq = $settings{'fastq'};
open IN, $fastq or die $!;
my $sumLen = 0;
my $numberReadsToCheckLen = 50; 
for (my $i = 0; $i < $numberReadsToCheckLen; $i++)
{
	my $id = <IN>;
	my $seq = <IN>;
	my $readLength = length ($seq);
	$sumLen += $readLength;
	my $plus = <IN>;
	my $qual = <IN>;
}
my $averageReadLength = int($sumLen/$numberReadsToCheckLen);

my $root = "$settingsFile.CRISPRlungo";

my $inputForGenomeMapping = $settings{'fastq'};
my $customMappedBamFile = "$root.customMapped.bam"; #file only generated if cut sites provided
my $customIndexFasta = "$root.customIndex.fa";
my $customAlignedCount = 0;

my $primerInfo = $settings{'primer_site'};
my $primerChr = "";
my $primerLoc = -1;
if ($primerInfo)
{
	($primerChr, $primerLoc) = split "_", $primerInfo;
}

my $printedCustomAmpliconsCount = 0;

if (@cuts > 0)
{
	#first create an index with predicted outcomes
	open OUT, ">$customIndexFasta" or die $!;

	for (my $i = 0; $i < @cuts; $i++)
	{
		#first, add wildtype AA
		my ($chrA, $siteA) = split "_", $cuts[$i];
		my $cutStartA = $siteA - ($averageReadLength + $alignmentExtension);
		my $cutStartAStop = $siteA - 1;
		my $cutEndA = $siteA + ($averageReadLength + $alignmentExtension);

		my $primerIsInLeftBitA = ($primerChr != "" && $primerChr == $chrA && $primerLoc >= $cutStartA && $primerLoc <= $cutStartAStop);
		my $primerIsInRightBitA = ($primerChr != "" && $primerChr == $chrA && $primerLoc >= $siteA && $primerLoc <= $cutEndA);

		my $leftBitA = `samtools faidx -n 10000 $settings{'genome'} $chrA:$cutStartA-$cutStartAStop | tail -n 1`;
		chomp $leftBitA;
		my $leftBitARC = reverse($leftBitA); #reverse complement
		$leftBitARC =~ tr/ACGTacgt/TGCAtgca/;

		my $rightBitA = `samtools faidx -n 10000 $settings{'genome'} $chrA:$siteA-$cutEndA | tail -n 1`;
		chomp $rightBitA;
		my $rightBitARC = reverse($rightBitA); #reverse complement
		$rightBitARC =~ tr/ACGTacgt/TGCAtgca/;

		my $wtSeq = $leftBitA.$rightBitA;
		chomp $wtSeq;
		if ($primerChr eq "" || $primerIsInLeftBitA || $primerIsInRightBitA)
		{
			print OUT ">CRISPRlungo_wt$i\n$wtSeq\n";
			$printedCustomAmpliconsCount++;
		}

		#next add splices with other cuts AB and BA
		for (my $j = $i+1; $j < @cuts; $j++)
		{
			my ($chrB, $siteB) = split "_", $cuts[$j];
			my $cutStartB = $siteB - ($averageReadLength + 50);
			my $cutStartBStop = $siteB - 1;
			my $cutEndB = $siteB + ($averageReadLength + 50);


			my $primerIsInLeftBitB = ($primerChr != "" && $primerChr == $chrB && $primerLoc >= $cutStartB && $primerLoc <= $cutStartBStop);
			my $primerIsInRightBitB = ($primerChr != "" && $primerChr == $chrB && $primerLoc >= $siteB && $primerLoc <= $cutEndB);

			my $leftBitB = `samtools faidx -n 10000 $settings{'genome'} $chrB:$cutStartB-$cutStartBStop | tail -n 1`;
			chomp $leftBitB;
			my $leftBitBRC = reverse($leftBitB); #reverse complement
			$leftBitBRC =~ tr/ACGTacgt/TGCAtgca/;

			my $rightBitB = `samtools faidx -n 10000 $settings{'genome'} $chrB:$siteB-$cutEndB | tail -n 1`;
			chomp $rightBitB;
			my $rightBitBRC = reverse($rightBitB); #reverse complement
			$rightBitBRC =~ tr/ACGTacgt/TGCAtgca/;

			if ($primerChr eq "" || $primerIsInLeftBitA || $primerIsInRightBitB)
			{
				my $LARB = $leftBitA.$rightBitB;
				print OUT ">CRISPRlungo_L$i.R$j\n$LARB\n";
				$printedCustomAmpliconsCount++;
			}

			if ($primerChr eq "" || $primerIsInLeftBitB || $primerIsInRightBitA)
			{
				my $LBRA = $leftBitB.$rightBitA;
				print OUT ">CRISPRlungo_L$j.R$i\n$LBRA\n";
				$printedCustomAmpliconsCount++;
			}

			if ($primerChr eq "" || $primerIsInLeftBitA || $primerIsInLeftBitB)
			{
				my $LALB = $leftBitA.$leftBitBRC;
				print OUT ">CRISPRlungo_L$i.L$j\n$LALB\n";
				$printedCustomAmpliconsCount++;
			}
			if ($primerChr eq "" || $primerIsInRightBitA || $primerIsInRightBitB)
			{
				my $RARB = $rightBitA.$rightBitBRC;
				print OUT ">CRISPRlungo_R$i.R$j\n$RARB\n";
				$printedCustomAmpliconsCount++;
			}
		}
	}
	close OUT;
	print "Created $printedCustomAmpliconsCount custom amplicons\n";
	print "Creating custom index\n";
	unless (-e "$customIndexFasta.1.bt2")
	{
		print `bowtie2-build $customIndexFasta $customIndexFasta`;
	}

	print "Aligning reads to custom index\n";
	#next, align to the entire genome
	my $customUnmappedFastqFile = "$root.customUnmapped.fastq";
	my $customAlnCommand = "bowtie2 -x $customIndexFasta -U $fastq --un-gz $customUnmappedFastqFile 2> $root.customBowtie2Log | samtools view -Shu - | samtools sort -o $customMappedBamFile - && samtools index $customMappedBamFile";

	if (!-e $customMappedBamFile)
	{
		print `$customAlnCommand`;
	}
	$inputForGenomeMapping = $customUnmappedFastqFile;
	$customAlignedCount = `samtools view -F 4 -c $customMappedBamFile`;
	chomp $customAlignedCount;
}

print "Aligning remainder of reads to genome\n";
#next, align to the entire genome
my $genomeMappedBamFile = "$root.genomeMapped.bam";
my $alnCommand = "bowtie2 -x $settings{'bowtie2Genome'} -U $inputForGenomeMapping 2> $root.genomeBowtie2Log | samtools view -Shu - | samtools sort -o $genomeMappedBamFile - && samtools index $genomeMappedBamFile";

if (!-e $genomeMappedBamFile)
{
	#print `$alnCommand 2> $root.bowtie2Log`;
	print `$alnCommand`;
}

my $genomeAlignedCount = `samtools view -F 4 -c $genomeMappedBamFile`;
chomp $genomeAlignedCount;

print "Creating chopped reads\n";
#next, run through all reads, and chop unaligned reads
open IN, "samtools view $genomeMappedBamFile |" or die $!;
my $unmappedFile = "$root.unmapped.fastq";
open UNMAPPED, ">$unmappedFile" or die $!;
my %mappedChrs = ();
my $unmappedID = 0;
my %alignedLocs = (); #aligned reads will not be chopped, but keep track of where they aligned for CRISPResso output
my %alignedChrCounts = ();
my $maxFrags = 0; 
my %fragsPerRead = ();
while (my $line = <IN>)
{
	my @lineEls = split "\t", $line;
	my $lineChr = $lineEls[2];
	$mappedChrs{$lineChr}++;

	my $mapq = $lineEls[5];
	my $unmapped = $lineEls[1] & 0x4;#unmapped
	my $start = $lineEls[3]-1;#position


	#chop unmapped reads
	if ($unmapped)
	{
		my $seq = $lineEls[9];#sequence
		my $qual = $lineEls[10];
		my $id = $lineEls[0];
		my $fragNum = 0;
		my $offset = 0;
		while ($offset + $fragmentSize < length($seq))
		{
			my $newID = "\@$id.CLID$unmappedID.CLO$offset.CLF$fragNum";
			my $newSeq = substr($seq,$offset,$fragmentSize);
			my $newQual = substr($qual,$offset,$fragmentSize);
			print UNMAPPED "$newID\n$newSeq\n+\n$newQual\n";

			$fragNum++;
			$offset += $fragmentStepSize;
		}

		$offset = length($seq)-$fragmentSize;
		my $newID = "\@$id.CLID$unmappedID.CLO$offset.CLF$fragNum";
		my $newSeq = substr($seq,$offset,$fragmentSize);
		my $newQual = substr($qual,$offset,$fragmentSize);
		print UNMAPPED "$newID\n$newSeq\n+\n$newQual\n";
		$maxFrags = $fragNum if $fragNum > $maxFrags;

		$fragsPerRead{$fragNum}++;

		$unmappedID++;
	}
	else
	{
		$alignedLocs{$lineChr}{$start}++;
		$alignedChrCounts{$lineChr}++;
	}
}
close UNMAPPED;
my $numberUnmappedReadsChopped = $unmappedID;
print "Created chopped reads for $numberUnmappedReadsChopped reads\n";
open OUT, ">$root.fragsPerUnalignedRead" or die $!;
print OUT "numFragments\tnumReads\n";
foreach my $numFrags (sort {$a<=>$b} keys %fragsPerRead)
{
	print OUT "$numFrags\t$fragsPerRead{$numFrags}\n";
}
close OUT;
open OUT, ">$root.globalAlignedChrs" or die $!;
print OUT "chr\tnumReads\n";
foreach my $chr (sort keys %alignedChrCounts)
{
	print OUT "$chr\t$alignedChrCounts{$chr}\n";
}
close OUT;


print "Aligning chopped reads\n";
#align the chopped reads
my $mappedChoppedSamFile = "$root.mappedChopped.sam";
my $alnCommand2 = "bowtie2 -x $settings{'bowtie2Genome'} -U $unmappedFile -S $mappedChoppedSamFile --end-to-end";

if (!-e $mappedChoppedSamFile)
{
	print `$alnCommand2 2>$root.bowtie2ChoppedLog`;
}


print "Processing chopped reads\n";
my $fragMetaFile = "$root.fragMeta";
open OUT, ">$fragMetaFile" or die $!;
my $head = "ID";
for (my $i = 0; $i < $maxFrags; $i++)
{
	$head .= "\t$i";
}
print OUT "$head\n";

open IN, $mappedChoppedSamFile or die $!;
my $line1 = <IN>;
while ($line1 and $line1 !~ /^\@PG/)
{
	$line1 = <IN>;
}
my $currID = "";
my %currIDChrs = (); #keep track of which chrs frags from this id map to
my %currIDVals = (); #keep track of where each frag maps to (this is a hash but it's like an array, so we can check for uninitialized values)
my %translocations = (); #hash of all translocation locations
my $translocationCount = 0;
my $unidentifiedCount = 0;
my %chromsPerFragReadCount = (); #keep track for all reads how many chrs the fragments mapped to
my $fragsMappedCount = 0;
while (my $line = <IN>)
{
	my @lineEls = split "\t", $line;
	my $lineChr = $lineEls[2];
	$mappedChrs{$lineChr}++;
	$fragsMappedCount++;

	my $mapq = $lineEls[5];
	my $unmapped = $lineEls[1] & 0x4;#unmapped
	my $start = $lineEls[3]-1;#position
	
	my $orig_id;
	my $lungoID;
	my $lungoOffset;
	my $lungoFrag;

	my $id = $lineEls[0];
	if ($id =~ /(.*)\.CLID(\d+)\.CLO(\d+)\.CLF(\d+)$/)
	{
		$orig_id = $1;
		$lungoID = $2;
		$lungoOffset = $3;
		$lungoFrag = $4;
	}
	else
	{
		die "Can't parse id: $id from line $line. Perhaps the line was trimmed?\n";
	}

	#if the current line is the next id, print the prevous line 
	if ($currID ne $orig_id and $currID ne "")
	{
		my $outline = $currID;
		for (my $i = 0; $i < $maxFrags; $i++)
		{
			my $val = "";
			if (exists ($currIDVals{$i}))
			{
				$val = $currIDVals{$i};
				my ($valChr, $valPos) = split " ", $currIDVals{$i};
				$currIDChrs{$valChr}++ unless $valChr eq "*";
			}
			$outline .= "\t$val";
		}
		print OUT "$outline\n";

		if ($currIDVals{"0"} =~ /\*/)
		{
			$unidentifiedCount++;
		}
		elsif ($currIDVals{"lastChr"} =~ /\*/)
		{
			$unidentifiedCount++;
		}
		else
		{
			my ($startChr, $startPos) = split " ", $currIDVals{"0"};
			$translocations{"$startChr $currIDVals{'lastChr'}"}++;
			$translocationCount++;
		}

		my $currChrCount = keys %currIDChrs;
		$chromsPerFragReadCount{$currChrCount}++;
		%currIDChrs = ();
		%currIDVals = ();
	}
	$currID = $orig_id;

	my $inferredStart = $start - $lungoOffset;
	$currIDVals{$lungoFrag} = "$lineChr $inferredStart";
	$currIDVals{'last'} = "$lineChr $inferredStart";
	$currIDVals{'lastChr'} = $lineChr;
}

if ($fragsMappedCount > 0)
{
	#print last one
	my $outline = $currID;
	for (my $i = 0; $i < $maxFrags; $i++)
	{
		my $val = "";
		if (exists ($currIDVals{$i}))
		{
			$val = $currIDVals{$i};
			my ($valChr, $valPos) = split " ", $currIDVals{$i};
			$currIDChrs{$valChr}++;
		}
		$outline .= "\t$val";
	}
	print OUT "$outline\n";
	if ($currIDVals{"0"} =~ /\*/)
	{
		$unidentifiedCount++;
	}
	else
	{
		my ($startChr, $startPos) = split " ", $currIDVals{"0"};
		$translocations{"$startChr $currIDVals{'lastChr'}"}++;
		$translocationCount++;
	}
	my $currChrCount = keys %currIDChrs;
	$chromsPerFragReadCount{$currChrCount}++;
	#done with last one
}

close OUT;
print "Found $translocationCount translocations and $unidentifiedCount unidentified reads\n";

open OUT, ">$root.fragsAlignedChrs" or die $!;
print OUT "chr\tnumReads\n";
foreach my $chr (sort keys %mappedChrs)
{
	print OUT "$chr\t$mappedChrs{$chr}\n";
}
close OUT;

open OUT, ">$root.fragChromsPerRead" or die $!;
print OUT "numChroms\tnumReads\n";
foreach my $numChroms (sort {$a<=>$b} keys %chromsPerFragReadCount)
{
	print OUT "$numChroms\t$chromsPerFragReadCount{$numChroms}\n";
}
close OUT;

print "Reporting\n";
#report
my %translocationTable = ();
open OUT, ">$root.translocationReport" or die $!;
print OUT "from\tto\tcount\n";
foreach my $key (sort keys %translocations)
{
	my ($fromChr, $toChr) = split " ", $key;
	print OUT "$fromChr\t$toChr\t$translocations{$key}\n";
	$translocationTable{$fromChr}{$toChr} += $translocations{$key};
}
close OUT;
#
#report
open OUT, ">$root.translocationReport.table" or die $!;
my @keys = sort keys %translocationTable;
my $translocationHead = "data\t".join("\t", @keys);
print OUT "$translocationHead\n";
foreach my $key (@keys)
{
	my $line = "$key";
	foreach my $key2 (@keys)
	{
		my $val = 0;
		$val = $translocationTable{$key}{$key2} if exists $translocationTable{$key} and exists $translocationTable{$key}{$key2};
		$line .= "\t$val";
	}
	print OUT "$line\n";
}
close OUT;

print "Preparing for CRISPResso\n";
my $lenFile = $settings{'genome'};
$lenFile =~ s/.fa$/.dict/;

open IN, $lenFile or die $!;
my %chrLens = ();
my @chrs = ();
while (my $line = <IN>)
{
	my @lineEls = split "\t", $line;
	my $chr = $lineEls[1];
	next if $chr !~ /SN:chr/;
	$chr =~ s/SN://;
	my $len = $lineEls[2];
	$len =~ s/LN://;
	$chrLens{$chr} = $len;
	push @chrs, $chr;
}
close IN;

#merge overlapping regions and create crispreso runs for reads aligned to genome
my $crispressoCutoff = $settings{'crispresso2_min_count'};
open INFO, ">$root.CRISPResso.info" or die $!;
print INFO "name\tchr\tstart\tend\treadCount\tamplicon\n";
my @crispressoNames = ();
my @crispressoCommands = ();
my $crispressoCount = 0;
foreach my $chr (@chrs)
{
	my @startLocs = (0);
	my @endLocs = (0);
	my $lastInd = 0;
	my @keys = sort keys %{$alignedLocs{$chr}};
	foreach my $start (sort{$a<=>$b} keys %{$alignedLocs{$chr}})
	{
		my $count = $alignedLocs{$chr}{$start} || 0;
		if ($count > $crispressoCutoff)
		{
			if ($endLocs[$lastInd] > $start)
			{
				$endLocs[$lastInd] = $start + $averageReadLength;
			}
			else
			{
				$lastInd++;
				$startLocs[$lastInd] = $start;
				$endLocs[$lastInd] = $start + $averageReadLength;
			}
		}
	}

	for (my $i = 1; $i < @startLocs; $i++)
	{
		my $name = $chr."_".$startLocs[$i];
		push @crispressoNames, $name;
		my $ampSeq = `samtools faidx -n 10000 $settings{'genome'} $chr:$startLocs[$i]-$endLocs[$i] | tail -n 1`;
		chomp $ampSeq;
		
		my $readCount = 0;
		open READS, "samtools view -F 4 $genomeMappedBamFile $chr:$startLocs[$i]-$endLocs[$i] | " or die $!;
		my $readsFile = "$root.crispresso.$name.fastq";
		open READSOUT, ">$readsFile" or die $!;
		while (my $line = <READS>)
		{
			my @lineEls = split "\t", $line;
			my $seq = $lineEls[9];#sequence
			my $qual = $lineEls[10];
			my $id = $lineEls[0];
			print READSOUT "\@$id\n$seq\n+\n$qual\n";
			$readCount++;
		}
		close READSOUT;
		close READS;
		print INFO "$name\t$chr\t$startLocs[$i]\t$endLocs[$i]\t$readCount\t$ampSeq\n";

		my $command = "CRISPResso -n $name -a $ampSeq -r1 $readsFile &> $readsFile.log";
		push @crispressoCommands, $command;
		$crispressoCount++;

	}
}

#create crispresso runs for reads aligned to custom index
my $customTranslocationCount = 0;
my $customWildtypeCount = 0;
if (@cuts > 0)
{
	my $queryBpAroundCut = 50; # pull reads out that are within 50bp of cut site
	my @customChrNames = ();
	my @customNames1 = ();
	my @customNames2 = ();
	for (my $i = 0; $i < @cuts; $i++)
	{
		push @customChrNames, "CRISPRlungo_wt$i";
		push @customNames1, $cuts[$i];
		push @customNames2, $cuts[$i];
		for (my $j = $i+1; $j < @cuts; $j++)
		{
			push @customChrNames, ("CRISPRlungo_L$i.R$j", "CRISPRlungo_L$j.R$i", "CRISPRlungo_L$i.L$j", "CRISPRlungo_R$i.R$j");
			push @customNames1, ($cuts[$i])x4;
			push @customNames2, ($cuts[$j])x4;
		}
	}

	my $startInd = ($averageReadLength + $alignmentExtension) - $queryBpAroundCut;
	my $endInd = ($averageReadLength + $alignmentExtension) + $queryBpAroundCut;

	for (my $i = 0; $i < @customChrNames; $i++)
	{
		my $customChrName = $customChrNames[$i];

		my $readCount = `samtools view -F 4 -c $customMappedBamFile $customChrName:$startInd-$endInd`;
		chomp $readCount;
		if ($customChrName =~ /CRISPRlungo_wt/)
		{
			$customWildtypeCount+= $readCount;
		}
		else
		{
			$customTranslocationCount+= $readCount;
		}
		if ($readCount > $settings{'crispresso2_min_count'})
		{
			open READS, "samtools view -F 4 $customMappedBamFile $customChrName:$startInd-$endInd | " or die $!;
			my $readsFile = "$root.crispresso.$customChrName.fastq";
			open READSOUT, ">$readsFile" or die $!;
			my $printedReadCount = 0;
			while (my $line = <READS>)
			{
				my @lineEls = split "\t", $line;
				my $seq = $lineEls[9];#sequence
				my $qual = $lineEls[10];
				my $id = $lineEls[0];
				print READSOUT "\@$id\n$seq\n+\n$qual\n";
				$printedReadCount++;
			}
			close READSOUT;
			close READS;
			my $ampSeq = `samtools faidx -n 10000 $customIndexFasta $customChrName:$startInd-$endInd | tail -n 1`;
			chomp $ampSeq;

			print INFO "$customChrName\t$customChrName\t$customNames1[$i]\t$customNames2[$i]\t$readCount\t$ampSeq\n";

			my $command = "CRISPResso -n $customChrName -a $ampSeq -r1 $readsFile &> $readsFile.log";
			push @crispressoCommands, $command;
			$crispressoCount++;
		}
	}
		
}

open OUT, ">$root.log" or die $!;
print OUT "aligned_custom\twildtype_custom\ttranslocations_custom\taligned_global\tunaligned_chopped\ttranslocations\tunidentified\n";
print OUT "$customAlignedCount\t$customWildtypeCount\t$customTranslocationCount\t$genomeAlignedCount\t$numberUnmappedReadsChopped\t$translocationCount\t$unidentifiedCount\n";
close OUT;

open OUT, ">$root.CRISPResso.commands" or die $!;
print OUT join("\n", @crispressoCommands);
close OUT;
open OUT, ">$root.CRISPResso.commands" or die $!;
print "Created $crispressoCount crispresso runs\n";
for (my $i = 0; $i < @crispressoCommands; $i++)
{
	print "Running command $i\n";
	my $output = `$crispressoCommands[$i]`;
	print OUT $output;
}
close OUT;

print "Finished\n";
