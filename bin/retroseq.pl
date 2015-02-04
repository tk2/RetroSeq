#!/usr/bin/env perl
=pod
This file is part of RetroSeq.

    RetroSeq is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RetroSeq is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RetroSeq.  If not, see <http://www.gnu.org/licenses/>.
=cut

use Carp;
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;

use lib dirname(__FILE__).'/..';
use RetroSeq::Vcf;
use RetroSeq::Utilities;

my $DEFAULT_ID = 80;
my $DEFAULT_LENGTH = 36;
my $DEFAULT_ANCHORQ = 20;
my $DEFAULT_MAX_DEPTH = 200;
my $DEFAULT_READS = 10;
my $MAX_READ_GAP_IN_REGION = 120;
my $DEFAULT_MIN_CLUSTER_READS = 2;
my $DEFAULT_MAX_CLUSTER_DIST = 4000;

my $HEADER = qq[#retroseq v:].$RetroSeq::Utilities::VERSION;
my $FOOTER = qq[#END_CANDIDATES];

my $BAMFLAGS = 
{
    'paired_tech'    => 0x0001,
    'read_paired'    => 0x0002,
    'unmapped'       => 0x0004,
    'mate_unmapped'  => 0x0008,
    'reverse_strand' => 0x0010,
    'mate_reverse'   => 0x0020,
    '1st_in_pair'    => 0x0040,
    '2nd_in_pair'    => 0x0080,
    'not_primary'    => 0x0100,
    'failed_qc'      => 0x0200,
    'duplicate'      => 0x0400,
};

my ($discover, $call, $bam, $bams, $ref, $eRefFofn, $length, $id, $output, $anchorQ, $region, $input, $reads, $depth, $noclean, $tmpdir, $readgroups, $filterFile, $orientate, $ignoreRGsFofn, $callNovel, $refTEs, $excludeRegionsDis, $doAlign, $help);

GetOptions
(
    #actions
    'discover'      =>  \$discover,
    'call'          =>  \$call,
    
    #parameters
    'bam=s'         =>  \$bam,
    'bams=s'        =>  \$bams,
    'eref=s'        =>  \$eRefFofn,
    'ref=s'         =>  \$ref,
    'len=s'         =>  \$length,
    'id=s'          =>  \$id,
    'q=s'           =>  \$anchorQ,
    'output=s'      =>  \$output,
    'input=s'       =>  \$input,
    'reads=s'       =>  \$reads,
    'depth=s'       =>  \$depth,
    'noclean'       =>  \$noclean,
    'region=s'      =>  \$region,
    'tmp=s'         =>  \$tmpdir,
    'rgs=s'         =>  \$readgroups,
    'filter=s'      =>  \$filterFile,
    'orientate=s'   =>  \$orientate,
    'ignoreRGs=s'   =>  \$ignoreRGsFofn,
    'novel'         =>  \$callNovel,
    'refTEs=s'      =>  \$refTEs,
    'exd=s'         =>  \$excludeRegionsDis,
    'align'         =>  \$doAlign,
    'h|help'        =>  \$help,
);

print <<MESSAGE;

RetroSeq: A tool for discovery of transposable elements from short read alignments

Version: $RetroSeq::Utilities::VERSION
Author: Thomas Keane (thomas.keane\@sanger.ac.uk)

MESSAGE

my $USAGE = <<USAGE;
Usage: $0 -<command> options

            -discover       Takes a BAM and a set of reference TE (fasta) and calls candidate supporting read pairs (BED output)
            -call           Takes multiple output of discovery stage and a BAM and outputs a VCF of TE calls
            
NOTE: $0 requires samtools, bcftools, exonerate, unix sort, bedtools to be in the default path

READ NAMES: This software assumes that reads that make up a pair have identical read names in the BAM file

USAGE

( $discover || $call || $help) or die $USAGE;

if( $discover )
{
    ( $bam && $output ) or die <<USAGE;
Usage: $0 -discover -bam <string> -eref <string> -output <string> [-q <int>] [-id <int>] [-len <int> -noclean]
    
    -bam        BAM file of paired reads mapped to reference genome
    -output     Output file to store candidate supporting reads (required for calling step)
    [-refTEs    Tab file with TE type and BED file of reference elements. These will be used to quickly assign discordant reads the TE types and avoid alignment. Using this will speed up discovery dramatically.]
    [-noclean   Do not remove intermediate output files. Default is to cleanup.]
    [-q         Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-id        Minimum percent ID for a match of a read to the transposon references. Default is 90.]
    [-rgs       Comma separated list of readgroups to operate on. Default is all.]
    [-exd       Fofn of BED files of regions where discordant mates falling into will be excluded e.g. simple repeats, centromeres, telomeres]
    [-align     Do the computational expensive exonerate PE discordant mate alignment step]
    [-eref      Tab file with list of transposon types and the corresponding fasta file of reference sequences (e.g. SINE   /home/me/refs/SINE.fasta). Required when the -align option is used.]
    [-len       Minimum length of a hit to the transposon references when using the -align option. Default is 36bp.]
    
USAGE
    
    croak qq[Cant find BAM file: $bam] unless -f $bam || -l $bam;
    croak qq[Cant find TE probes sequence file (-eref parameter)] if( $doAlign && ( ! $eRefFofn || ! -f $eRefFofn ) );
    
    my $erefs;
    if( $doAlign || $eRefFofn )
    {
        print qq[Reading -eref file: $eRefFofn\n];
        $erefs = _tab2Hash( $eRefFofn );
        foreach my $type ( keys( %{$erefs} ) )
        {
            if( ! -f $$erefs{$type} ){croak qq[Cant find transposon reference file: ].$$erefs{ $type };}
        }
        if( ! $doAlign ){print qq[Forcing -align option to be switched on as transposon sequence files were provided\n];$doAlign = 1;}
    }
    
    my $refTEsF;
    if( $refTEs )
    {
        print qq[Reading -refTEs file: $refTEs\n];
        $refTEsF = _tab2Hash( $refTEs );
        foreach my $type ( keys( %{$refTEsF} ) )
        {
            if( ! -f $$refTEsF{$type} ){croak qq[Cant find $type reference TEs file: ].$$refTEsF{ $type };}
        }
    }
    
    $anchorQ = defined( $anchorQ ) && $anchorQ > -1 ? $anchorQ : $DEFAULT_ANCHORQ;
    $id = defined( $id ) && $id < 101 && $id > 0 ? $id : $DEFAULT_ID;
    $length = defined( $length ) && $length > 25 ? $length : $DEFAULT_LENGTH;
    my $clean = defined( $noclean ) ? 0 : 1;
        
    if( $readgroups && length( $readgroups ) > 0 )
    {
        my @s = split( /,/, $readgroups );foreach my $rg ( @s ){if( $rg !~ /[A-Za-z0-9]|\.|-|_+/ ){croak qq[Invalid readgroup: $rg\n];}}
        print qq[Restricting discovery phase to read groups: $readgroups\n];
    }
    
    print qq[\nMin anchor quality: $anchorQ\nMin percent identity: $id\nMin length for hit: $length\n\n];
    
    #test for samtools
    RetroSeq::Utilities::checkBinary( q[samtools], qq[0.1.16], qq[0.1.19] );
    RetroSeq::Utilities::checkBinary( q[exonerate], qq[2.2.0] ) if( $doAlign );
    RetroSeq::Utilities::checkBinary( q[bedtools] );
    
    _findCandidates( $bam, $erefs, $id, $length, $anchorQ, $output, $readgroups, $refTEsF, $excludeRegionsDis, $doAlign, $clean );
}
elsif( $call )
{
    ( $bam && $input && $ref && $output ) or die <<USAGE;
Usage: $0 -call -bam <string> -input <string> -ref <string> -output <string> [ -filter <BED file> -cleanup -reads <int> -depth <int>]
    
    -bam            BAM file OR BAM fofn
    -input          Either a single output file from the PE discover stage OR a prefix of a set of files from discovery to be combined for calling OR a fofn of discovery stage output files
    -ref            Fasta of reference genome
    -output         Output file name (VCF)
    [-filter        Tab file with TE type and BED file of reference elements. These will be filtered out from the calling.]
    [-region        Call a particular chromosome only (chr) OR region (chr:start-end) only]
    [-depth         Max average depth of a region to be considered for calling. Default is 200.]
    [-reads         It is the minimum number of reads required to make a call. Default is 5.]
    [-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30.]
    [-ignoreRGs     Single read group name OR a file of readgroups that should be ignored. Default is none.]
    [-noclean       Do not remove intermediate output files. Default is to cleanup.]
    
USAGE
    
    croak qq[Cant find BAM or BAM fofn: $bam] unless -f $bam || -l $bam;
    
    if( ! -f $input )
    {
        my @files = glob( qq[$input*] );
        if( ! @files || @files == 0 )
        {
            croak qq[Cant find input or files with prefix: $input];
        }
    }
    croak qq[Cant find reference genome fasta: $ref] unless -f $ref;
    croak qq[Cant find reference genome index - please index your reference file with samtools faidx] unless -f qq[$ref.fai];
    
    if( $reads ){die qq[Invalid reads parameter: $reads] unless ( $reads =~ /^\d+$/ && $reads > -1 ); print qq[Min reads: $reads\n];}else{$reads = $DEFAULT_READS;}
    $reads = defined( $reads ) && $reads =~ /^\d+$/ && $reads > -1 ? $reads : $DEFAULT_READS;
    $depth = defined( $depth ) && $depth =~ /^\d+$/ && $depth > -1 ? $depth : $DEFAULT_MAX_DEPTH;
    $anchorQ = defined( $anchorQ ) && $anchorQ > -1 ? $anchorQ : $DEFAULT_ANCHORQ;
    my $clean = defined( $noclean ) ? 0 : 1;
    my %filterBEDs;
    if( $filterFile )
    {
        croak qq[Cant find filter file: $filterFile] unless -f $filterFile;
        
        open( my $tfh, $filterFile ) or die $!;
        while(my $line = <$tfh> ){chomp( $line );my @s = split( /\t/, $line );die qq[Cant find $s[0] filter file: $s[ 1 ]\n] unless -f $s[ 1 ];$filterBEDs{ $s[ 0 ] } = $s[ 1 ];}
    }
    
    #test for samtools
    RetroSeq::Utilities::checkBinary( q[samtools], qq[0.1.16] );
    RetroSeq::Utilities::checkBinary( q[bcftools] );
    RetroSeq::Utilities::checkBinary( q[bedtools] );
    
    my @bams;
    my $first = `head -1 $bam`;chomp( $first );
    if( -f $first ) #is this a fofn
    {
        open( my $ifh, $bam ) or die qq[failed to BAM fofn: $input\n];
        while(my $file = <$ifh> )
        {
            chomp( $file );
            if( ( -l $file || -f $file ) && ( -l $file.qq[.bai] || -f $file.qq[.bai] ) ){push(@bams, $file);}else{die qq[Cant find BAM input file or BAM index file: $file\n];}
        }
        print qq[Found ].scalar(@bams).qq[ BAM files\n\n];
        close( $ifh );
    }
    else{if( ( -l $bam || -f $bam ) && ( -l $bam.qq[.bai] || -f $bam.qq[.bai] ) ){push( @bams, $bam );}else{die qq[Cant find BAM input file or BAM index file: $bam\n];}}
    
    my $sampleName = RetroSeq::Utilities::getBAMSampleName( \@bams );
    print qq[Calling sample $sampleName\n];
    

    print qq[Beginning paired-end calling...\n];
    _findInsertions( \@bams, $sampleName, $input, $ref, $output, $reads, $depth, $anchorQ, $region, $clean, \%filterBEDs, 1, $orientate, $ignoreRGsFofn, $callNovel );
	exit;
}
else
{
    print qq[You did not specify an action!\n\n$USAGE];
}

print qq[RetroSeq finished successfully\n\n];
exit;

sub _findCandidates
{
    my $bam = shift;
    my $erefs = shift;
    my $id = shift;
    my $length = shift;
    my $minAnchor = shift;
    my $output = shift;
    my $readgroups = shift;
    my $refTEsF = shift;
    my $excludeRegionsDis = shift;
    my $doAlign = shift;
    my $clean = shift;
    
    my $readgroupsFile = qq[$$.readgroups];
    if( $readgroups )
    {
        open( my $tfh, qq[>$readgroupsFile] ) or die $!;
        foreach my $rg (split(/,/, $readgroups )){print $tfh qq[$rg\n];}
        close( $tfh );
    }
    
    my $fastaCounter = 0;
    my $candidatesFasta = $doAlign ? qq[$$.candidates.$fastaCounter.fasta] : undef;
    my $candidatesBed = qq[$$.candidate_anchors.bed];
    my $discordantMatesBed = qq[$$.discordant_mates.bed];
    my $clipFasta = qq[$$.clip_candidates.fasta];
    my %candidates = %{ _getCandidateTEReadNames($bam, $readgroups, $minAnchor, $candidatesFasta, $candidatesBed, $discordantMatesBed ) };
    
    print scalar( keys( %candidates ) ).qq[ candidate reads remain to be found after first pass....\n];
    
    if( keys( %candidates ) > 0 )
    {
        my $ffh;
        if( $doAlign ){open( $ffh, qq[>>$candidatesFasta] ) or die qq[ERROR: Failed to create fasta file: $!\n];}
        
        #now go and get the reads from the bam (annoying have to traverse through the bam a second time - but required for reads where the mate is aligned distantly)
        #also dump out their mates as will need these later as anchors
        open( my $bfh, qq[samtools view ].( defined( $readgroups ) ? qq[-R $$.readgroups ] : qq[ ] ).qq[$bam |] ) or die $!;
        open( my $dfh, qq[>>$discordantMatesBed] ) or die $!;
        my $currentChr = '';
        my $readsFound = 0;
        while( my $sam = <$bfh> )
        {
            chomp( $sam );
            my @s = split( /\t/, $sam );
            my $name = $s[ 0 ];
            my $flag = $s[ 1 ];
            my $ref = $s[ 2 ];
            if( $candidates{ $name } )
            {
                #       read is 1st in pair                looking for 1st in pair      read in 2nd in pair                     looking for second
                if( ($flag & $$BAMFLAGS{'1st_in_pair'}) && $candidates{ $name } == 1 || ($flag & $$BAMFLAGS{'2nd_in_pair'}) && $candidates{ $name } == 2 )
                {
                    my $seq = $s[ 9 ];
                    if( _sequenceQualityOK( $seq ) )
                    {
                        if( $doAlign ){print $ffh qq[>$name\n$seq\n];}
                        
                        #record the read position details in the BED file of discordant mates
                        my $pos = $s[ 3 ];
                        my $qual = $s[ 4 ];
                        my $dir = ($flag & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
                        my $endPos = RetroSeq::Utilities::getSAMendpos($pos,$s[5]);
                        print $dfh qq[$ref\t$pos\t$endPos\t$name\t$dir\t$qual\n];
                        $readsFound ++;
                    }
                    delete( $candidates{ $name } );
                }
            }
            if( $currentChr ne $ref && $ref ne '*' ){print qq[Reading chromosome: $ref\n];$currentChr = $ref;}
        }
        if( $doAlign ){close( $ffh );}close( $dfh );
        
        if( $readsFound == 0 )
        {
            die qq[ERROR: Failed to recover discordant anchor reads - are your mate read names identical? RetroSeq assumes that mates have identical read names.];
        }
    }
    undef %candidates; #finished with this hash
    
    #user provided a file of regions where discordant mates should be excluded
    if( $excludeRegionsDis && -f $excludeRegionsDis )
    {
        #cat the files together
        open( my $tfh, $excludeRegionsDis ) or die $!;
        while( my $f = <$tfh> )
        {
            chomp( $f );
            if( ! -f $f ){die qq[ERROR: Cant find exclude file: $f\n];}
            system(qq[cat $f >> $$.regions_exclude.bed]) == 0 or die qq[Failed to cat files together\n];
        }
        close( $tfh );
        
        #bedtools intersect with the mates BED file - collect the readnames
        system( qq[bedtools intersect -a $discordantMatesBed -b $$.regions_exclude.bed -f 0.5 -u > $$.discordant_mates.remove.bed] ) == 0 or die qq[Failed to run bedtools intersect to filter];
        my %reads;
        open( my $rfh, qq[bedtools intersect -a $discordantMatesBed -b $$.regions_exclude.bed -f 0.5 -u | awk -F"\t" '{print \$4}' | ] ) or die $!;
        while( my $l = <$rfh> ){chomp( $l );$reads{ $l } = 1;}
        
        #filter the fasta file on these readnames
        if( $doAlign )
        {
            open( $tfh, $candidatesFasta ) or die $!;
            $fastaCounter ++; my $newFasta = qq[$$.candidates.$fastaCounter.fasta];
            open( my $tofh, qq[>$newFasta] ) or die $!;
            while( my $h = <$tfh> )
            {
                chomp( $h );
                my $seq = <$tfh>;
                chomp( $seq );
                if( ! $reads{ substr( $h, 1 ) } )
                {
                    print $tofh qq[$h\n$seq\n];
                }
            }
            close( $tfh );close( $tofh );
            $candidatesFasta = $newFasta;
        }
    }
    
    #if the user provided a tab file with mappings of TE type to BED file of locations
    #then use this info to assign discordant mates to these types
    if( $refTEsF )
    {
        print qq[Using reference TE locations to assign discordant mates...\n];
        foreach my $type(keys(%$refTEsF))
        {
            my $file = $$refTEsF{ $type };
            print qq[Screening for hits to: $type\n];
            system( qq[bedtools intersect -a $discordantMatesBed -b $file -u | awk -F"\t" '{print \$4,\$5}' > $$.$type.mates.bed] ) == 0 or die qq[Failed to run bedtools intersect];
            
            #print the mates (i.e. the anchors) of these reads into the discovery output file
            #first load up the readnames
            my %reads;
            open( my $t1fh, qq[$$.$type.mates.bed] ) or die $!;
            while( my $l = <$t1fh> )
            {
                chomp( $l );my @s = split( /\s/, $l );
                $reads{ $s[0] } = $s[1]; #stick the read names in a hash
            }
            close( $t1fh );
            
            #now get the reads from the anchors BED file and print them
            open( my $bfh, qq[$candidatesBed] ) or die $!;
            open( my $tofh, qq[>>$output] ) or die $!;
            while( my $l = <$bfh> )
            {
                chomp( $l );
                my @s = split( /\t/, $l );
                if( $reads{ $s[ 3 ] } )
                {
                    print $tofh qq[$s[0]\t$s[1]\t$s[2]\t$type\t$s[3]\t$s[4]\t].$reads{ $s[ 3 ] }.qq[\t90\n];
                }
            }
            close( $bfh );close( $tofh );
            
            #now filter the fasta file to remove these reads if doing alignment step subsequently
            if( $doAlign )
            {
                open( my $tfh, $candidatesFasta ) or die $!;
                $fastaCounter ++; my $newFasta = qq[$$.candidates.$fastaCounter.fasta];
                open( $tofh, qq[>$newFasta] ) or die $!;
                while( my $h = <$tfh> )
                {
                    chomp( $h );
                    my $seq = <$tfh>;
                    chomp( $seq );
                    if( ! $reads{ substr( $h, 1 ) } )
                    {
                        print $tofh qq[$h\n$seq\n];
                    }
                }
                close( $tfh );close( $tofh );
                $candidatesFasta = $newFasta;
            }
        }
    }
    
    if( $doAlign )
    {
        #create a single fasta file with all the TE ref seqs
        my $refsFasta = qq[$$.allrefs.fasta];
        open( my $tfh, qq[>$refsFasta] ) or die $!;
        foreach my $type ( keys( %{ $erefs } ) )
        {
            open( my $sfh, $$erefs{ $type } ) or die $!;
            my $seqCount = 1;
            while( my $line = <$sfh>)
            {
                chomp( $line );
                if( $line =~ /^>/ )
                {
                    print $tfh qq[>$type\n]; #call them all by the type label
                }
                else
                {
                    print $tfh qq[$line\n];
                }
            }
            close( $sfh );
        }
        close( $tfh );
        
        #if there arent any candidates, then we are done
        if( -s $candidatesFasta == 0 )
        {
            if( $clean )
            {
                #delete the intermediate files
                unlink( glob( qq[$$.*] ) ) or die qq[Failed to remove intermediate files: $!];
            }
            print qq[Failed to find any candidate reads - exiting];
            exit;
        }
        
        #run exonerate and parse the output from the stream (dump out hits for different refs to diff temp files)
        #output format:                                                    read hitlen %id +/- refName
        open( my $efh, qq[exonerate -m affine:local --bestn 5 --ryo "INFO: %qi %qal %pi %tS %ti\n"].qq[ $candidatesFasta $refsFasta | egrep "^INFO|completed" | ] ) or die qq[Exonerate failed to run: $!];
        print qq[Parsing PE alignments....\n];
        my $lastLine;
        my %anchors;
        while( my $hit = <$efh> )
        {
            chomp( $hit );
            $lastLine = $hit;
            last if ( $hit =~ /^-- completed/ );
            
            my @s = split( /\s+/, $hit );
            #    check min identity	  check min length
            if( $s[ 3 ] >= $id && $s[ 2 ] >= $length )
            {
                $anchors{ $s[ 1 ] } = [ $s[ 5 ], $s[ 4 ], $s[ 3 ] ]; #TE type, mate orientation, percent ID. This could be a memory issue (possibly dump out to a file and then use unix join to intersect with the bed)
            }
        }
        close( $efh );
        
        if( $lastLine ne qq[-- completed exonerate analysis] ){die qq[Alignment did not complete correctly\n];}
        
        #now put all the anchors together into a single file per type
        open( my $afh, qq[>>$output] ) or die $!;
        open( my $cfh, $candidatesBed ) or die $!;
        while( my $anchor = <$cfh> )
        {
            chomp( $anchor );
            my @s = split( /\t/, $anchor );
            if( defined( $anchors{ $s[ 3 ] } ) )
            {
                my $mate_orientation = $anchors{ $s[ 3 ] }[ 1 ];
                my $type = $anchors{ $s[ 3 ] }[ 0 ];
                #note the anchor orientation is contained is 3rd last, then mate orientation 2nd entry, and then percent id is last
                print $afh qq[$s[0]\t$s[1]\t$s[2]\t$type\t$s[3]\t$s[4]\t$mate_orientation\t].$anchors{ $s[ 3 ] }[ 2 ].qq[\n];
            }
        }
        close( $afh );close( $cfh );
        undef( %anchors );
    }
    
    if( $clean )
    {
        #delete the intermediate files
        unlink( glob( qq[$$.*] ) ) or die qq[Failed to remove intermediate files: $!];
    }
}

sub _findInsertions
{
    my $bamsRef = shift;my @bams = @{$bamsRef};
    my $sampleName = shift;
    my $input = shift;
    my $ref = shift;
    my $output = shift;
    my $minReads = shift;
    my $depth = shift;
    my $minQ = shift;
    my $region = shift;
    my $clean = shift;
    my $tempR  = shift;my %filterBEDs; %filterBEDs = %{$tempR} if( defined( $tempR ) );
    my $hets = shift;
    my $orientate = shift;
    my $ignoreRGs = shift;
    my $novel = shift;
    
    my @files;
    if( -f $input )
    {
        my $first = `head -1 $input`;chomp( $first );
        if( -f $first ) #is this a fofn
        {
            open( my $ifh, $input ) or die qq[failed to open fofn of discovery output files: $input\n];
            while(my $file = <$ifh> )
            {
                chomp( $file );
                if( -f $file ){push(@files, $file);}else{die qq[Cant find discovery output file: $file\n];}
            }
            print qq[Found ].scalar(@files).qq[ PE discovery stage input files\n\n];
            close( $ifh );
        }else{push( @files, $input );}#looks like a discovery output file
    }
    else
    {
        @files = glob( qq[$input*] ) or die qq[Failed to glob files: $input*\n];
        die qq[Cant find any inputs files with prefix $input*\n] if( ! @files || @files == 0 );
    }
    
    #parse the region
    my ($chr, $start, $end) = (undef, undef, undef);
    if( defined( $region ) )
    {
        if( $region =~ /^([A-Za-z0-9]+):([0-9]+)-([0-9]+)$/ )
        {
            $chr = $1;$start=$2;$end=$3;
            print qq[Restricting calling to region: $region\n];
        }
        else{$chr = $region;print qq[Restricting calling to Chr: $chr\n];}
    }
    
    #merge the discovery files together and sort them by type, chr, position cols
    my $merged = qq[$$.merge.PE.tab];
    foreach my $file( @files )
    {
        system(qq[cat $file >> $$.merge.PE.tab]) == 0 or die qq[Failed to merge PE input files\n];
    }
    
    if( $chr )
    {
        if( $start && $end )
        {
            system(qq[echo -e "$chr\t$start\t$end" > $$.region.bed ]) == 0 or die qq[failed to defined region];
            system(qq/awk '\$1=="$chr"&&\$2>$start&&\$3<$end' $merged > $$.merge.PE.region.tab/) == 0 or die qq[failed to grep chr out from reads file];
            $merged = qq[$$.merge.PE.region.tab];
        }
        else #just the chr is defined
        {
            system(qq/awk '\$1=="$chr"' $merged > $$.merge.PE.region.tab/) == 0 or die qq[failed to grep chr out from reads file];
            $merged = qq[$$.merge.PE.region.tab];
        }
    }
    
    my $sortedCandidates = qq[$$.merge.PE.sorted.tab];
    system( qq[sort -k4,4d -k1,1d -k2,2n $merged | uniq > $sortedCandidates] ) == 0 or die qq[Failed to sort merged input files\n];
    
    my $ignoreRGsFormatted = undef;
    if( $ignoreRGs )
    {
        $ignoreRGsFormatted = qq[$$.ignore.rgs];
        open( my $tfh, qq[>$ignoreRGsFormatted] ) or die qq[Failed to create tmp file: $ignoreRGsFormatted\n];
        if( -f $ignoreRGs )
        {
            open( my $ifh, $ignoreRGs )  or die qq[Failed to open readgroups file: $ignoreRGs $!\n];
            while( my $rg = <$ifh> ){chomp( $rg );print $tfh qq[RG:Z:$rg\t\n];print qq[Ignoring RG: $rg\n];}
        }
        else
        {
            print $tfh qq[RG:$ignoreRGs\t\n];
            print qq[Ignoring RG: $ignoreRGs\n];
        }
        print qq[\n];
        close( $tfh );
    }
    
    #pull out the reads categorised as 'unknown' (if any exist in the file) into a separate file
    my $unknownCandidatesFile = qq[$$.merge.PE.sorted.unknown.tab];
    system(qq[awk -F"\t" '{if(\$4=="unknown"){print \$1"\t"\$2"\t"\$3"\tunknown\t"\$6}}' $sortedCandidates > $$.merge.PE.sorted.unknown.tab]) == 0 or die qq[Failed to extract unknown candidate reads\n];
    
    open( my $tfh, $sortedCandidates ) or die qq[Failed to open merged tab file: $!\n];
    my %typeBEDFiles;
    my $currentType = '';
    my $cfh;
    my $readCount = 0;
    my $raw_candidates = qq[$output.candidates];
    open( my $dfh, qq[>$raw_candidates] ) or die $!;
    print $dfh qq[FILTER: chr\tstart\tend\ttype\_sample\tL_Fwd_both\tL_Rev_both\tL_Fwd_single\tL_Rev_single\tR_Fwd_both\tR_Rev_both\tR_Fwd_single\tR_Rev_single\tL_Last_both\tR_First_both\tDist\n];
    close( $dfh );
    my $currentTypeAnchorsFile;
    while( my $line = <$tfh> )
    {
        chomp( $line );
        my @s = split( /\t/, $line );
        if( $currentType eq '' )
        {
            $currentType = $s[ 3 ];
            die qq[Failed to detect the TE type to be called from line: $line\n] unless $currentType;
            
            print qq[PE Call: $currentType\n];
            $currentTypeAnchorsFile = qq[$$.$currentType.pe_anchors.bed];
            open( $cfh, qq[>$currentTypeAnchorsFile] ) or die $!;
        }
        elsif( $currentType ne $s[ 3 ] || eof( $tfh ) )
        {
            close( $cfh );
            
            if( $currentType ne 'unknown' || ($novel && $currentType eq 'unknown') )
            {
                #convert to a region BED (removing any candidates with very low numbers of reads)
                print qq[$currentType Calling initial rough boundaries of insertions....\n];
                my $rawTECalls1 = qq[$$.raw_calls.1.$currentType.tab];
                
                #add in the unknown mates (if any were present)
                if( -s $unknownCandidatesFile )
                {
                    my $newfile = qq[$$.$currentType.pe_anchors.inc_unknown.bed];
                    system(qq[cat $$.$currentType.pe_anchors.bed $unknownCandidatesFile | sort -k1,1d -k2,2n > $newfile]) == 0 or die qq[Failed to add in unknown candidate reads to $currentType file\n];
                    $currentTypeAnchorsFile = $newfile;
                }
                
                my $numRegions = RetroSeq::Utilities::convertToRegionBedPairsWindowBED( $currentTypeAnchorsFile, $DEFAULT_MIN_CLUSTER_READS, $currentType, $MAX_READ_GAP_IN_REGION, $DEFAULT_MAX_CLUSTER_DIST, 1, 0, $rawTECalls1 );
                
                if( $numRegions > 0 )
                {
                    #if a bed filter file of ref types was provided - then apply te filter now to reduce number of calls to test
                    if( %filterBEDs && $filterBEDs{ $currentType } )
                    {
                        print qq[Filtering reference elements for: $currentType\n];
                        
                        #remove the regions specified in the exclusion BED file
                        my $rawTECalls1Filtered = qq[$$.raw_calls.1.$currentType.filtered.tab];
                        RetroSeq::Utilities::filterOutRegions( $rawTECalls1, $filterBEDs{ $currentType }, $rawTECalls1Filtered );
                        $rawTECalls1 = $rawTECalls1Filtered;
                    }
                    
                    #remove extreme depth calls
                    print qq[Removing calls with extremely high depth (>$depth)....\n];
                    my $rawTECalls2 = qq[$$.raw_calls.2.$currentType.tab];
                    _removeExtremeDepthCalls( $rawTECalls1, \@bams, $depth, $rawTECalls2, $raw_candidates );
                    
                    #new calling filtering code
                    print qq[Filtering and refining candidate regions into calls....\n];
                    my $homCalls = qq[$$.raw_calls.3.$currentType.hom.bed];
                    my $hetCalls = qq[$$.raw_calls.3.$currentType.het.bed];
                    
                    _filterCallsBedMinima( $rawTECalls2, \@bams, 10, $minQ, $ref, $raw_candidates, $hets, $homCalls, $hetCalls, $ignoreRGsFormatted, $minReads );
                    
                    #remove close duplicated calls
                    my $rmdupHomCalls = qq[$$.raw_calls.3.$currentType.hom.rmdup.bed];
                    RetroSeq::Utilities::_removeDups( $homCalls, $rmdupHomCalls );
                    $typeBEDFiles{ $currentType }{hom} = $rmdupHomCalls if( -f $rmdupHomCalls && -s $rmdupHomCalls > 0 );
                    
                    my $rmdupHetCalls;
                    if( -s $hetCalls > 0 )
                    {
                        $rmdupHetCalls = qq[$$.raw_calls.3.$currentType.het.rmdup.bed];
                        RetroSeq::Utilities::_removeDups( $hetCalls, $rmdupHetCalls );
                        $typeBEDFiles{ $currentType }{het} = $rmdupHetCalls;
                    }
                    
                    #remove the temporary files for this 
                    if( $clean )
                    {
                        unlink( qq[$$.raw_reads.0.$currentType.tab], qq[$$.raw_calls.1.$currentType.tab], qq[$$.raw_calls.2.$currentType.tab] ) or die qq[Failed to delete temp files\n];
                        unlink( $homCalls ) if( -f $homCalls && -s $homCalls == 0 );
                        unlink( $hetCalls ) if ( -f $hetCalls && -s $hetCalls == 0 );
                    }
                }
            }
            
            last if( eof( $tfh ) );
            
            $readCount = 0;
            $currentType = $s[ 3 ];
            print qq[PE Call: $currentType\n];
            $currentTypeAnchorsFile = qq[$$.$currentType.pe_anchors.bed];
            open( $cfh, qq[>$currentTypeAnchorsFile] ) or die $!;
        }
        else
        {
            #             chr    start  stop   TE_type   a_orien  m_orien %id
            print $cfh qq[$s[0]\t$s[1]\t$s[2]\t$currentType\t$s[5]\n];
            $readCount ++;
        }
    }
    close( $tfh );
    close( $cfh ) if( $cfh );
    
    #identify hybrid insertions from unused clusters (in $$.hybrid.[pos|neg].out files)
    if( -s qq[$$.hybrid.pos.out] && -s qq[$$.hybrid.neg.out] )
    {
        print qq[PE Call: Hybrid elements\n];
        my $rawTECalls1 = qq[$$.raw_calls.1.hybrid.tab];
        open(my $ofh, qq[>$rawTECalls1] ) or die $!;
        my $window = $DEFAULT_MAX_CLUSTER_DIST/2;
        open(my $bfh, qq[bedtools window -a $$.hybrid.pos.out -b $$.hybrid.neg.out -l $window -r $window | ] ) or die qq[bedtools failed for hybrids];
        my $numRegions = 0;
        while( my $line=<$bfh> )
        {
            chomp($line);
            my @s = split( /\t/, $line );
            
            #do some sanity checks??
            my $totalReads = $s[ 5 ] + $s[ 11 ];
            my $end = $s[2]>$s[8]?$s[2]:$s[8];
            my $start = $s[1]<$s[7]?$s[1]:$s[7];
            
            print $ofh qq[$s[0]\t$start\t$end\t$s[3]-$s[9]-hybrid\t$totalReads\t].$RetroSeq::Utilities::PASS.qq[\n];$numRegions ++;
        }
        
        close($bfh);

        #now get the clusters where there is only support on either 5 prime or 3 prime (so these can be tested for unknown reads on other side)
        open( $bfh, qq[bedtools window -a $$.hybrid.pos.out -b $$.raw_calls.1.hybrid.tab -v | ] ) or die qq[bedtools failed for hybrids-unknown];
        while( my $line=<$bfh> )
        {
            chomp($line);
            my @s = split( /\t/, $line );
            
            my $totalReads = $s[ 5 ];
            if( $totalReads >= $minReads )
            {
                print $ofh qq[$s[0]\t$s[1]\t].($s[2]+500).qq[\t$s[3]-unknown\t$totalReads].$RetroSeq::Utilities::PASS.qq[\n];$numRegions ++;
            }
        }
        
        open( $bfh, qq[bedtools window -a $$.hybrid.neg.out -b $$.raw_calls.1.hybrid.tab -v | ] ) or die qq[bedtools failed for hybrids-unknown];
        while( my $line=<$bfh> )
        {
            chomp($line);
            my @s = split( /\t/, $line );

            my $totalReads = $s[ 5 ];
            if( $totalReads >= $minReads )
            {
                print $ofh qq[$s[0]\t].($s[1]-500).qq[\t$s[2]\tunknown-$s[3]\t$totalReads].$RetroSeq::Utilities::PASS.qq[\n];$numRegions ++;
            }
        }
        close($ofh);
        
        if( $numRegions > 0 )
        {
            #remove extreme depth calls
            print qq[Removing calls with extremely high depth (>$depth)....\n];
            my $rawTECalls2 = qq[$$.raw_calls.2.hybrid.tab];
            _removeExtremeDepthCalls( $rawTECalls1, \@bams, $depth, $rawTECalls2, $raw_candidates );
            
            #new calling filtering code
            print qq[Filtering and refining candidate regions into calls....\n];
            my $homCalls = qq[$$.raw_calls.3.hybrid.hom.bed];
            my $hetCalls = qq[$$.raw_calls.3.hybrid.het.bed];
          
            _filterCallsBedMinima( $rawTECalls2, \@bams, 10, $minQ, $ref, $raw_candidates, $hets, $homCalls, $hetCalls, $ignoreRGsFormatted, $minReads );
            
            #remove close duplicated calls
            my $rmdupHomCalls = qq[$$.raw_calls.3.hybrid.hom.rmdup.bed];
            RetroSeq::Utilities::_removeDups( $homCalls, $rmdupHomCalls );
            $typeBEDFiles{ hybrid }{hom} = $rmdupHomCalls if( -f $rmdupHomCalls && -s $rmdupHomCalls > 0 );
         
            my $rmdupHetCalls;
            if( -s $hetCalls > 0 )
            {
                $rmdupHetCalls = qq[$$.raw_calls.3.hybrid.het.rmdup.bed];
                RetroSeq::Utilities::_removeDups( $hetCalls, $rmdupHetCalls );
                $typeBEDFiles{ hybrid }{het} = $rmdupHetCalls;
            }
            
            #remove the temporary files for this 
            if( $clean )
            {
                unlink( qq[$$.raw_reads.0.hybrid.tab], qq[$$.raw_calls.1.hybrid.tab], qq[$$.raw_calls.2.hybrid.tab] ) or die qq[Failed to delete temp files\n];
                unlink( $homCalls ) if( -f $homCalls && -s $homCalls == 0 );
                unlink( $hetCalls ) if ( -f $hetCalls && -s $hetCalls == 0 );
            }
        }
    }
    
    if( $novel )
    {
        #remove all the unknown calls that overlap with calls that are assigned to a type
        my $unknownHoms = $typeBEDFiles{unknown}{hom};
        my $unknownHets = $typeBEDFiles{unknown}{het};
        foreach my $type( keys( %typeBEDFiles ) )
        {
            next if( $type eq 'unknown' );
            my $f = $typeBEDFiles{$type}{hom};
            system( qq[cat $f >> $$.PE.allcalls.bed] ) == 0 or die qq[Failed to cat file: $f] if( $f && -s $f );
            
            $f = $typeBEDFiles{$type}{het}; 
            system( qq[cat $f >> $$.PE.allcalls.bed] ) == 0 or die qq[Failed to cat file: $f] if( $f && -s $f );
        }
        
        if( -f $unknownHoms )
        {
            my $removed = RetroSeq::Utilities::filterOutRegions( $unknownHoms, qq[$$.PE.allcalls.bed], qq[$$.unknownHoms.filtered] );
            if( $removed > 0 )
            {
                print qq[Removed $removed hom unknown calls\n];
                $typeBEDFiles{unknown}{hom} = qq[$$.unknownHoms.filtered];
            }
        }
        
        if( defined( $unknownHets ) && -s $unknownHets )
        {
            my $removed = RetroSeq::Utilities::filterOutRegions( $unknownHets, qq[$$.PE.allcalls.bed], qq[$$.unknownHets.filtered] );
            if( $removed > 0 )
            {
                print qq[Removed $removed het unknown calls\n];
                $typeBEDFiles{unknown}{het} = qq[$$.unknownHets.filtered];
            }
        }
    }
    
    #output calls in VCF and BED format
    print qq[Creating VCF file of calls....\n];
    _outputCalls( \%typeBEDFiles, $sampleName, $ref, $output );
    
	if( $clean )
	{
	    #delete the intermediate files
	    unlink( glob( qq[$$.*] ) ) or die qq[Failed to remove intermediate files: $!];
	}
}

=pod
the new and improved calling code from mouse paper
takes a BED of rough call regions and refines them into breakpoints and does checking of the supporting
reads either side
=cut
sub _filterCallsBedMinima
{
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 11;
    
    my $bedin = shift;
	my @bams = @{ $_[ 0 ] };shift;
	my $minDepth = shift;
	my $minMapQ = shift;
	my $ref = shift;
	my $thrown_out_file = shift;
	my $hets = shift;
	my $bedoutHoms = shift;
	my $bedoutHets = shift;
	my $ignoreRGsFormatted = shift;
	my $minReads = shift;
	
	open( my $ifh, $bedin ) or die $!;
	open( my $homsfh, qq[>$bedoutHoms] ) or die $!;
	my $hetsfh = undef;
	open( $hetsfh, qq[>$bedoutHets] ) or die $!;
	
	open( my $dfh, qq[>>$thrown_out_file] ) or die $!;
	my %breakpointsTested;
	while( my $originalCall = <$ifh> )
	{
	    chomp( $originalCall );
	    print qq[Considering call $originalCall\n];
	    
	    my @originalCallA = split( /\t/, $originalCall );
	    my $strain = (split( /\./, $originalCallA[3] ))[ 0 ];
	    
	    my $start = $originalCallA[ 1 ];my $end = $originalCallA[ 2 ];my $chr = $originalCallA[ 0 ];
	    
	    my $t = RetroSeq::Utilities::getCandidateBreakPointsDirVote($originalCallA[0], $originalCallA[1], $originalCallA[2],\@bams,20 );
	    if( ! $t ){warn qq[Failed to get candidate breakpoints for call $originalCall\n];next;}
	    my $breakpoint = $t->[ 0 ];
	    
	    if( $breakpointsTested{ $breakpoint } ){print qq[Duplicate breakpoint found: $breakpoint. Skipping.\n];next;}else{$breakpointsTested{$breakpoint} = 1;}
	    
        print qq[Testing breakpoint $breakpoint\n];
        
        my $result = RetroSeq::Utilities::testBreakPoint( $originalCallA[ 0 ], $breakpoint, \@bams, $minMapQ, $originalCall, $dfh, $ignoreRGsFormatted, $minReads, 0 );
        
        my $flag = $result->[0];
        my $call = $result->[1];
        my $ratio = $result->[2];
        my $spanningRPs = $result->[4];
        
        if( $flag == $RetroSeq::Utilities::INV_BREAKPOINT )
        {
            print qq[Failed - at inversion breakpoint\n];
        }
        elsif( $flag > $RetroSeq::Utilities::NOT_ENOUGH_READS_CLUSTER )
        {
            #note to self - disabled the hom/het decision code (needs to be rewritten)
            print qq[Breakpoint score: $ratio Flag: $flag\n];
            print $homsfh $call.qq[\t$flag\t$spanningRPs\n];
        }
        else{print qq[Discarding: $call : $flag\n];}
	}
	close( $ifh );
	close( $homsfh );
	close( $hetsfh );
	close( $dfh );
}

#***************************INTERNAL HELPER FUNCTIONS********************

sub _getCandidateTEReadNames
{
    my $bam = shift;
    my $readgroups = shift;
    my $minAnchor = shift;
    my $candidatesFasta = shift;
    my $candidatesBed = shift;
    my $discordantMatesBed = shift;
    
    my %candidates;
    my $ffh;
    if( $candidatesFasta ){open( $ffh, qq[>>$candidatesFasta] );}
    open( my $afh, qq[>$candidatesBed] ) or die qq[ERROR: Failed to create anchors file: $!\n];
    my $sebfh;
    
    my $cfh;
    open( my $dfh, qq[>>$discordantMatesBed] ) or die qq[ERROR: Failed to create discordant mates BED: $!\n];
    
    print qq[Opening BAM ($bam) and getting initial set of candidate mates....\n];
    
    open( my $bfh, qq[samtools view ].( defined( $readgroups ) ? qq[-R $$.readgroups ] : qq[ ] ).qq[$bam |] ) or die $!;
    my $currentChr = '';
    while ( my $samLine = <$bfh> )
    {
        chomp( $samLine );
        my @sam = split( /\t/, $samLine );
        my $flag = $sam[ 1 ];
        my $qual = $sam[ 4 ];
        my $name = $sam[ 0 ];
        my $ref = $sam[ 2 ];
        my $mref = $sam[ 6 ];
        my $cigar = $sam[ 5 ];
        
        if( $candidates{ $name } )
        {
            if( ($flag & $$BAMFLAGS{'1st_in_pair'}) && $candidates{ $name } == 1 || ($flag & $$BAMFLAGS{'2nd_in_pair'}) && $candidates{ $name } == 2 )
            {
                my $seq = $sam[ 9 ];
                
                if( _sequenceQualityOK( $seq ) )
                {
                    if( defined($candidatesFasta) ){print $ffh qq[>$name\n$seq\n];}

                    #record the read position details in the BED file of discordant mates
                    my $pos = $sam[ 3 ];
                    my $dir = ($flag & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
                    my $endPos = RetroSeq::Utilities::getSAMendpos($pos,$cigar);
                    print $dfh qq[$ref\t$pos\t$endPos\t$name\t$dir\t$qual\n];
                }
                delete( $candidates{ $name } );
            }
        }
        
        my $supporting = RetroSeq::Utilities::isSupportingClusterRead( $flag, $sam[ 8 ], $qual, $minAnchor, $cigar, undef, undef, undef );
        
        if( $supporting > 0 )
        {
            if( $supporting == 1 ) #single-ended mapped
            {
               $candidates{ $name } = $flag & $$BAMFLAGS{'1st_in_pair'} ? 2 : 1; #mate is recorded
               
               my $pos = $sam[ 3 ];
               my $dir = ($flag & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
               my $endPos = RetroSeq::Utilities::getSAMendpos($pos,$sam[5]);
               print $afh qq[$ref\t$pos\t$endPos\t$name\t$dir\n];
            }
            elsif( $supporting == 2 ) #discordantly mapped
            {
               $candidates{ $name } = $flag & $$BAMFLAGS{'1st_in_pair'} ? 2 : 1; #mate is recorded
                   
               my $pos = $sam[ 3 ];
               my $dir = ($flag & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
               my $endPos = RetroSeq::Utilities::getSAMendpos($pos,$sam[5]);
               print $afh qq[$ref\t$pos\t$endPos\t$name\t$dir\n];
            }
        }
        if( $currentChr ne $ref && $ref ne '*' ){print qq[Reading chromosome: $ref\n];$currentChr = $ref;}
    }
    close( $bfh );
    if( defined($candidatesFasta) ){close( $ffh );}
    close( $afh );close( $dfh );
    
    return \%candidates;
}

sub _checkAvgQuality
{
    my $quals = shift;
    my $minAvg = 15;
    
    my @chars = split( //, $quals );
    my $sum = 0;
    foreach my $q (@chars)
    {
        $sum += ( unpack( 'C', $q ) - 33 );
    }
    
    return ($sum/length($quals)) >= $minAvg ? 1 : 0;
}

sub _checkDiscoveryOutput
{
    my $file = shift;
    
    open( my $tfh, $file ) or die $!;
    my $line = <$tfh>;chomp( $line );
    if( $line !~ /^$HEADER/ ){die qq[Malformed header of input file: $file\n$line\n$HEADER\n];}
    my $lastLine;
    while(<$tfh>){chomp;$lastLine=$_;}
    close( $tfh );
    
    if( $lastLine ne $FOOTER ){die qq[Malformed footer of input file: $file\n];}
    return 1;
}

#input is a list of samples and a bed file, output is a VCF + BED file of calls
sub _outputCalls
{
    my $t = shift; my %typeBedFiles = %{$t}; #BED/tab format
    my $sample = shift;
    my $reference = shift;
    my $output = shift;
    
    open( my $vfh, qq[>$output] ) or die $!;
    
    my $vcf_out = RetroSeq::Vcf->new();
    $vcf_out->add_columns($sample);
    
    my $header = RetroSeq::Utilities::getVcfHeader($vcf_out);
    print $vfh $header;
    
    open( my $bfh, qq[>$output.bed] ) or die $!; 
    foreach my $type ( keys( %typeBedFiles ) )
    {
        #check for hom calls
        if( $typeBedFiles{ $type }{hom} )
        {
            die qq[Cant find BED file for type: $type\n] if( ! -f $typeBedFiles{ $type }{hom} );
            open( my $cfh, $typeBedFiles{ $type }{hom} ) or die qq[Failed to open TE calls file: ].$typeBedFiles{ $type }.qq[\n];
            while( <$cfh> )
            {
                chomp;
                my @s = split( /\t/, $_ );
                my $pos = int( ( $s[ 1 ] + $s[ 2 ] ) / 2 );
                my $refbase = _getReferenceBase( $reference, $s[ 0 ], $pos );
                next if( ! $refbase );
                
                my $ci1 = $s[ 1 ] - $pos;
                my $ci2 = $s[ 2 ] - $pos;
                my $flag = $s[ 5 ];
                my $spanning = $s[ 6 ];
                
                print $bfh qq[$s[0]\t$pos\t].($pos+1).qq[\t$type==$sample\t$s[4]\tNA\n];
                
                my %out;
                $out{CHROM}  = $s[ 0 ];
                $out{POS}    = $pos;
                $out{ID}     = '.';
                $out{ALT}    = ['<INS:ME>'];
                $out{REF}    = $refbase;
                $out{QUAL}   = $s[ 4 ];
                $out{INFO} = { SVTYPE=>'INS', NOT_VALIDATED=>undef, MEINFO=>qq[$s[3],$s[1],$s[2],NA] };
                $out{FORMAT} = ['GT', 'GQ', 'FL', 'SP'];
                
                $out{gtypes}{$sample}{GT} = qq[<INS:ME>/<INS:ME>];
                $out{gtypes}{$sample}{GQ} = qq[$s[4]];
                $out{gtypes}{$sample}{FL} = $flag;
                $out{gtypes}{$sample}{SP} = $spanning;
                
                $vcf_out->format_genotype_strings(\%out);
                print $vfh $vcf_out->format_line(\%out);
            }
            close( $cfh );
        }
        
        if( $typeBedFiles{ $type }{het} )
        {
            die qq[Cant find BED file for type: $type\n] if( ! -f $typeBedFiles{ $type }{het} );
            open( my $cfh, $typeBedFiles{ $type }{het} ) or die qq[Failed to open TE calls file: ].$typeBedFiles{ $type }{het}.qq[\n];
            while( <$cfh> )
            {
                chomp;
                my @s = split( /\t/, $_ );
                my $pos = int( ( $s[ 1 ] + $s[ 2 ] ) / 2 );
                my $refbase = _getReferenceBase( $reference, $s[ 0 ], $pos );
                next if( ! $refbase );
                
                my $ci1 = $s[ 1 ] - $pos;
                my $ci2 = $s[ 2 ] - $pos;
                my $flag = $s[ 5 ];
                my $spanning = $s[ 6 ];
                
                print $bfh qq[$s[0]\t$pos\t].($pos+1).qq[\t$type==$sample\t$s[4]\tNA\n];
                
                my %out;
                $out{CHROM}  = $s[ 0 ];
                $out{POS}    = $pos;
                $out{ID}     = '.';
                $out{ALT}    = [$refbase,'<INS:ME>'];
                $out{REF}    = $refbase;
                $out{QUAL}   = $s[ 4 ];
                $out{INFO} = { SVTYPE=>'INS', NOT_VALIDATED=>undef, MEINFO=>qq[$s[3],$s[1],$s[2],NA] };
                $out{FORMAT} = ['GT','GQ', 'FL', 'SP'];
                
                $out{gtypes}{$sample}{GT} = qq[$refbase/<INS:ME>];
                $out{gtypes}{$sample}{GQ} = qq[$s[4]];
                $out{gtypes}{$sample}{FL} = $flag;
                $out{gtypes}{$sample}{SP} = $spanning;
                
                $vcf_out->format_genotype_strings(\%out);
                print $vfh $vcf_out->format_line(\%out);
            }
            close( $cfh );
        }
    }
    close( $vfh );close( $bfh );
}

#remove calls where there is very high depth (measured from pileup) in the region
sub _removeExtremeDepthCalls
{
	my $calls = shift;
	my @bams = @{ $_[ 0 ] };shift;
	my $maxDepth = shift;
	my $outputbed = shift;
	my $thrown_out_candidates_file = shift;
	
	my $bamStr = '';
	foreach my $bam( @bams ){$bamStr.=qq[ $bam];}
	
	open( my $cfh, $calls ) or die $!;
	open( my $ofh, ">$outputbed" ) or die $!;
	open( my $rfh, qq[>>$thrown_out_candidates_file] ) or die $!;
	while( <$cfh> )
	{
		chomp;
		
		my @s = split( /\t/, $_ );
		
		my $start = $s[ 1 ] - 100 > 0 ? $s[ 1 ] - 50 : 0;
		my $end = $s[ 2 ] + 100;
		my $size = $end - $start;
		my $chr = $s[ 0 ];
		my $label = $s[ 3 ];
		my $score = $s[ 4 ];
		
		#get the avg depth over the region
		my $totalDepth = `samtools mpileup -r $chr:$start-$end $bamStr | awk -F"\t" '{SUM += \$4} END {print SUM}'`;chomp( $totalDepth );
		my $avgDepth = $totalDepth / $size;
		if( $avgDepth < $maxDepth )
		{
			print $ofh qq[$_\n];
		}
		else
		{
		    print "Excluding call due to high depth: $_ AvgDepth: $avgDepth\n";
		    print $rfh qq[$chr\t$start\t$end\t$label\_depth\t$avgDepth\t$score\n];
		}
	}
	close( $cfh );
	close( $ofh );
	close( $rfh );
	
	return 1;
}

sub _revCompDNA
{
	croak "Usage: revCompDNA string\n" unless @_ == 1;
	my $seq = shift;
	$seq= uc( $seq );
	$seq=reverse( $seq );
	
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

#from an indexed fasta - get the reference base at a specific position
sub _getReferenceBase
{
    my $fasta = shift;
    my $chr = shift;
    my $pos = shift;
    
    my $base = `samtools faidx $fasta $chr:$pos-$pos | tail -1`;
    chomp( $base );
    
    if( length( $base ) != 1 || $base !~ /[acgtnACGTN]/ ){warn qq[Failed to get reference base at $chr:$pos-$pos: $base\n];return undef;}
    return $base;
}

sub _sequenceQualityOK
{
    my $seq = shift;
    
    #dont include the read if it is >80% N's
    $seq = uc( $seq );
    if( ($seq =~ tr/N//) < length( $seq ) * 0.8 && $seq =~ /^[ACGTN]+$/ )
    {
        return 1;
    }
    return 0;
}

sub _tab2Hash
{
    my $file = shift;
    
    my %hash;
    open( my $tfh, $file ) or die $!;
    while( my $entry = <$tfh> )
    {
        chomp( $entry );
        next if( $entry =~ /^#/ );
        die qq[Each line of tab delimited file should have two entries separated by single tab: $file\n] unless $entry =~ /^(.+)(\t)(.+)$/;
        $hash{ $1 } = $3;
    }
    close( $tfh );
    
    return \%hash;
}

