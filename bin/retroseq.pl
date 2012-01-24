#!/usr/bin/env perl
# 
# Author:       tk2
# Maintainer:   tk2
# Created:      Fri Sep 10 14:07:53 BST 2010 @588 /Internet Time/
# Updated:      Fri Sep 10 14:08:03 BST 2010 @588 /Internet Time/

=pod
Status: 
Should work and compile
Test: validation and genotyping steps

Issues/improvements:
Genotyping
    - test it
    - should dump out genotype likelihoods also
=cut

use Carp;
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Path qw(make_path);

use lib dirname(__FILE__).'/../lib/';
use Vcf;
use Utilities;

my $DEFAULT_ID = 80;
my $DEFAULT_LENGTH = 36;
my $DEFAULT_ANCHORQ = 20;
my $DEFAULT_MAX_DEPTH = 200;
my $DEFAULT_READS = 5;
my $DEFAULT_MIN_GENOTYPE_READS = 3;
my $MAX_READ_GAP_IN_REGION = 120;
my $DEFAULT_MIN_CLUSTER_READS = 2;
my $DEFAULT_SR_MIN_CLUSTER_READS = 3;
my $DEFAULT_SR_MIN_UNKNOWN_CLUSTER_READS = 10;
my $DEFAULT_MAX_CLUSTER_DIST = 4000;
my $DEFAULT_MAX_SR_CLUSTER_DIST = 30;
my $DEFAULT_MIN_SOFT_CLIP = 30;

my $HEADER = qq[#retroseq v:].substr($Utilities::VERSION,0,1);
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

my ($discover, $call, $genotype, $bam, $bams, $ref, $eRefFofn, $length, $id, $output, $anchorQ, $region, $input, $reads, $depth, $noclean, $tmpdir, $readgroups, $filterFile, $heterozygous, $orientate, $ignoreRGsFofn, $srmode, $minSoftClip, $srOutputFile, $srInputFile, $callNovel, $help);

GetOptions
(
    #actions
    'discover'      =>  \$discover,
    'call'          =>  \$call,
    'genotype'      =>  \$genotype,
    
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
    'hets'          =>  \$heterozygous,
    'filter=s'      =>  \$filterFile,
    'orientate=s'   =>  \$orientate,
    'ignoreRGs=s'   =>  \$ignoreRGsFofn,
    'srmode'        =>  \$srmode,
    'minclip=i'     =>  \$minSoftClip,
    'srcands=s'     =>  \$srOutputFile,
    'srInput=s'     =>  \$srInputFile,
    'novel'         =>  \$callNovel,
    'h|help'        =>  \$help,
);

print <<MESSAGE;

RetroSeq: A tool for discovery and genotyping of transposable elements from short read alignments
Version: $Utilities::VERSION
Author: Thomas Keane (thomas.keane\@sanger.ac.uk)

MESSAGE

my $USAGE = <<USAGE;
Usage: $0 -<command> options

            -discover       Takes a BAM and a set of reference TE (fasta) and calls candidate supporting read pairs (BED output)
            -call           Takes multiple output of discovery stage and a BAM and outputs a VCF of TE calls
            -genotype       Input is a VCF of TE calls and a set of sample BAMs, output is a new VCF with genotype calls for new samples
            
NOTE: $0 requires samtools, bcftools, exonerate, unix sort to be in the default path

USAGE

( $discover || $call || $genotype || $help) or die $USAGE;

if( $discover )
{
    ( $bam && $eRefFofn && $output ) or die <<USAGE;
Usage: $0 -discover -bam <string> -eref <string> -output <string> [-srmode] [-q <int>] [-id <int>] [-len <int> -noclean]
    
    -bam        BAM file of paired reads mapped to reference genome
    -eref       Tab file with list of transposon types and the corresponding fasta file of reference sequences (e.g. SINE   /home/me/refs/SINE.fasta)
    -output     Output file to store candidate supporting reads (required for calling step)
    [-srmode]   Search for split reads in the BAM file
    [-minclip]  Minimum length of soft clippped portion of read to be considered for split-read analysis. Default is 30bp.
    [-noclean   Do not remove intermediate output files. Default is to cleanup.]
    [-q         Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-id        Minimum percent ID for a match of a read to the transposon references. Default is 90.]
    [-len       Minimum length of a hit to the transposon references. Default is 36bp.]
    [-rgs       Comma separated list of readgroups to operate on. Default is all.]
    
USAGE
    
    croak qq[Cant find BAM file: $bam] unless -f $bam;
    croak qq[Cant find TE tab file: $eRefFofn] unless -f $eRefFofn;
    
    my $erefs = _tab2Hash( $eRefFofn );
    foreach my $type ( keys( %{$erefs} ) )
    {
        if( ! -f $$erefs{$type} ){croak qq[Cant find transposon reference file: ].$$erefs{ $type };}
    }
    
    $anchorQ = defined( $anchorQ ) && $anchorQ > -1 ? $anchorQ : $DEFAULT_ANCHORQ;
    $id = defined( $id ) && $id < 101 && $id > 0 ? $id : $DEFAULT_ID;
    $length = defined( $length ) && $length > 25 ? $length : $DEFAULT_LENGTH;
    my $clean = defined( $noclean ) ? 0 : 1;
    
    if( $srmode )
    {
        if( ! $srOutputFile ){die qq[You must specify the -srcands output file parameter in split-read mode\n];}
        if( ! $minSoftClip ){$minSoftClip = $DEFAULT_MIN_SOFT_CLIP;}
        print qq[Running split-read discovery mode\n];
    }else{undef($minSoftClip);}
    
    if( $readgroups && length( $readgroups ) > 0 )
    {
        my @s = split( /,/, $readgroups );foreach my $rg ( @s ){if( $rg !~ /[A-Za-z0-9]|\.|-|_+/ ){croak qq[Invalid readgroup: $rg\n];}}
        print qq[Restricting discovery phase to read groups: $readgroups\n];
    }
    
    print qq[\nMin anchor quality: $anchorQ\nMin percent identity: $id\nMin length for hit: $length\n\n];
    
    #test for samtools
    Utilities::checkBinary( q[samtools], qq[0.1.16] );
    Utilities::checkBinary( q[exonerate], qq[2.2.0] );
    
    _findCandidates( $bam, $erefs, $id, $length, $anchorQ, $output, $readgroups, $minSoftClip, $srOutputFile, $clean );
}
elsif( $call )
{
    ( $bam && $input && $ref && $output ) or die <<USAGE;
Usage: $0 -call -bam <string> -input <string> -ref <string> -output <string> [-srinput <SR candidates file> -filter <BED file> -cleanup -reads <int> -depth <int> -hets]
    
    -bam            BAM file OR BAM fofn
    -input          Either a single output file from the PE discover stage OR a prefix of a set of files from discovery to be combined for calling OR a fofn of discovery stage output files
    -ref            Fasta of reference genome
    -output         Output file name (VCF)
    [-srinput       Either a single output from split-read discovery stage OR a prefix of a set of files to be combined for calling (will trigger split-read calling to run)]
    [-hets          Call heterozygous insertions. Default is homozygous.]
    [-orientate     Attempt to predicate the orientation of the calls. Default is no.]
    [-filter        BED file of regions to exclude from final calls e.g. reference elements, low complexity etc.]
    [-region        Call a particular chromosome only (chr) OR region (chr:start-end) only]
    [-depth         Max average depth of a region to be considered for calling. Default is 200.]
    [-reads         It is the minimum number of reads required to make a call. Default is 5. Parameter is optional.]
    [-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-ignoreRGs     Single read group name OR a file of readgroups that should be ignored. Default is none.]
    [-novel]    Call novel sequence insertions also (from non-TE mate pairs)
    [-noclean       Do not remove intermediate output files. Default is to cleanup.]
    
USAGE
    
    croak qq[Cant find BAM or BAM fofn: $bam] unless -f $bam;
    
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
    if( $filterFile )
    {
        croak qq[Cant find filter file: $filterFile] unless -f $filterFile;
    }
    
    #test for samtools
    Utilities::checkBinary( q[samtools], qq[0.1.16] );
    Utilities::checkBinary( q[bcftools] );
    
    my @bams;
    my $first = `head -1 $bam`;chomp( $first );
    if( -f $first ) #is this a fofn
    {
        open( my $ifh, $bam ) or die qq[failed to BAM fofn: $input\n];
        while(my $file = <$ifh> )
        {
            chomp( $file );
            if( -f $file && -f $file.qq[.bai] ){push(@bams, $file);}else{die qq[Cant find BAM input file or BAM index file: $file\n];}
        }
        print qq[Found ].scalar(@bams).qq[ BAM files\n\n];
        close( $ifh );
    }
    else{if( -f $bam && -f $bam.qq[.bai] ){push( @bams, $bam );}else{die qq[Cant find BAM input file or BAM index file: $bam\n];}}
    
    my $sampleName = Utilities::getBAMSampleName( \@bams );
    print qq[Calling sample $sampleName\n];
    
    #what sort of calling to run - PE and/or SR
    if( $srInputFile )
    {
        print qq[Beginning split-read calling...\n];
        _findInsertionsSR( \@bams, $sampleName, $srInputFile, $ref, $output.qq[.SR], $reads, $depth, $region, $clean, $heterozygous, $ignoreRGsFofn, $callNovel );
        
        print qq[Beginning paired-end calling...\n];
        _findInsertions( \@bams, $sampleName, $input, $ref, $output.qq[.PE], $reads, $depth, $anchorQ, $region, $clean, $filterFile, $heterozygous, $orientate, $ignoreRGsFofn, $callNovel );
        
        #now merge the VCF files into a single VCF
        #implement later..
    }
    else
    {
        print qq[Beginning paired-end calling...\n];
        _findInsertions( \@bams, $sampleName, $input, $ref, $output.qq[.PE], $reads, $depth, $anchorQ, $region, $clean, $filterFile, $heterozygous, $orientate, $ignoreRGsFofn, $callNovel );
    }
	exit;
}
elsif( $genotype )
{
    ( $ref && $bams && $input && $output ) or die <<USAGE;
Usage: $0 -genotype -bams <string> -input <string> -ref <string> -output <string> [-cleanup -reads <int> -region <chr:[start-end]>
    
    -bams           File of BAM file names (one per sample to be genotyped)
    -input          VCF file of TE calls
    -ref            Fasta of reference genome
    -output         Output VCF file (will be annotated with new genotype calls)
    [-hets          Call heterozygous insertions. Default is homozygous.]
    [-orientate     Attempt to predict the orientation of the calls. Default is no.]
    [-region        Call a particular chromosome only (chr) OR region (chr:start-end) only]
    [-depth         Max average depth of a region to be considered for calling. Default is 200.]
    [-reads         It is the minimum number of reads required to make a call. Default is $DEFAULT_MIN_GENOTYPE_READS. Parameter is optional.]
    [-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-noclean       Do not remove intermediate output files. Default is to cleanup.]
USAGE
    
    croak qq[Cant find BAM fofn: $bams] unless -f $bams;
    croak qq[Cant find input: $input] unless -f $input;
    
    croak qq[Cant find reference genome fasta: $ref] unless -f $ref;
    croak qq[Cant find reference genome index - please index your reference file with samtools faidx] unless -f qq[$ref.fai];
    
    my $clean = defined( $noclean ) ? 0 : 1;
    $reads = defined( $reads ) && $reads =~ /^\d+$/ && $reads > -1 ? $reads : $DEFAULT_MIN_GENOTYPE_READS;
    $anchorQ = defined( $anchorQ ) && $anchorQ > -1 ? $anchorQ : $DEFAULT_ANCHORQ;
    
    #test for samtools
    Utilities::checkBinary( q[samtools], qq[0.1.16] );
    Utilities::checkBinary( q[bcftools] );
    
    _genotype( $bams, $input, $ref, $region, $anchorQ, $output, $clean, $heterozygous, $orientate, $reads );
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
    my $minSoftClip = shift;
    my $srOutput = shift;
    my $clean = shift;
    
    #test for exonerate
    Utilities::checkBinary( q[exonerate], qq[2.2.0] );
    
    my $readgroupsFile = qq[$$.readgroups];
    if( $readgroups )
    {
        open( my $tfh, qq[>$readgroupsFile] ) or die $!;
        foreach my $rg (split(/,/, $readgroups )){print $tfh qq[$rg\n];}
        close( $tfh );
    }
    
    my $candidatesFasta = qq[$$.candidates.fasta];
    my $candidatesBed = qq[$$.candidate_anchors.bed];
    my $clipFasta = qq[$$.clip_candidates.fasta];
    my %candidates = %{ _getCandidateTEReadNames($bam, $readgroups, $minAnchor, $minSoftClip, $candidatesFasta, $candidatesBed, $clipFasta ) };
    
    print scalar( keys( %candidates ) ).qq[ candidate reads remain to be found after first pass....\n];
    
    if( keys( %candidates ) > 0 )
    {
        open( my $ffh, qq[>>$candidatesFasta] ) or die qq[ERROR: Failed to create fasta file: $!\n];
        
        #now go and get the reads from the bam (annoying have to traverse through the bam a second time - but required for reads where the mate is aligned distantly)
        #also dump out their mates as will need these later as anchors
        open( my $bfh, qq[samtools view ].( defined( $readgroups ) ? qq[-R $$.readgroups ] : qq[ ] ).qq[$bam |] ) or die $!;
        my $currentChr = '';
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
                        print $ffh qq[>$name\n$seq\n];
                    }
                    delete( $candidates{ $name } );
                }
            }
            if( $currentChr ne $ref && $ref ne '*' ){print qq[Reading chromosome: $ref\n];$currentChr = $ref;}
        }
        close( $ffh );
    }
    undef %candidates; #finished with this hash
    
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
    open( my $efh, qq[exonerate -m affine:local --bestn 5 --ryo "INFO: %qi %qal %pi %tS %ti\n"].qq[ $$.candidates.fasta $refsFasta | egrep "^INFO|completed" | ] ) or die qq[Exonerate failed to run: $!];
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
    open( my $afh, qq[>$output] ) or die $!;
    open( my $cfh, qq[$$.candidate_anchors.bed] ) or die $!;
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
        else
        {
            print $afh qq[$s[0]\t$s[1]\t$s[2]\tunknown\t$s[3]\t$s[4]\n];
        }
    }
    close( $afh );close( $cfh );
    undef( %anchors );
    
    #if running split-read mode, then run exonerate on these sequences too
    if( $srmode && -s $clipFasta )
    {
        open( my $efh, qq[exonerate -m affine:local --bestn 5 --ryo "INFO: %qi %qal %pi %tS %ti\n"].qq[ $clipFasta $refsFasta | egrep "^INFO|completed" | uniq | ] ) or die qq[Exonerate failed to run: $!];
        print qq[Parsing split-read alignments....\n];
        open( my $cofh, qq[>$srOutput] ) or die qq[Failed to create clipped candidates output file: $!\n];
        $lastLine = '';
        my %readsAligned;
        while( my $hit = <$efh> )
        {
            chomp( $hit );
            $lastLine = $hit;
            last if ( $hit =~ /^-- completed/ );
            
            my @s = split( /\s+/, $hit );
            #    check min identity	  check min length
            if( $s[ 3 ] >= $id && $s[ 2 ] >= $length )
            {
                #get the readname
                my @readDetails = split( /###/, $s[ 1 ] );
                
                #               chr                pos              pos          TE     readname            flag             cigar        align_strand  insertSize
                my $print = qq[$readDetails[1]\t$readDetails[2]\t$readDetails[2]\t$s[5]\t$readDetails[0]\t$readDetails[3]\t$readDetails[4]\t$s[4]\t$s[3]\t$readDetails[5]\n];
                if( $print ne $lastLine )
                {
                    print $cofh $print;
                }
                $readsAligned{ $readDetails[ 0 ] } = 1;
                $lastLine = $print;
            }
        }
        close( $efh );
        
        #now get all the reads that didnt align to a TE and annotate them as unknown type
        open( my $tfh, $clipFasta )  or die $!;
        while( my $l = <$tfh> )
        {
            next unless $l =~ /^>/;
            $l = substr( $l, 1 );
            my @readDetails = split( /###/, $l );
            if( ! $readsAligned{ $readDetails[ 0 ] } )
            {
                my $print = qq[$readDetails[1]\t$readDetails[2]\t$readDetails[2]\tunknown\t$readDetails[0]\t$readDetails[3]\t$readDetails[4]\n];
                print $cofh $print;
            }
        }
        close( $cofh );
        
        if( $lastLine ne qq[-- completed exonerate analysis] ){die qq[SR alignment did not complete correctly\n];}
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
    my $filterBED  = shift;
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
    
    #merge the SR discovery files together and sort them by type, chr, position cols
    foreach my $file( @files )
    {
        system(qq[cat $file >> $$.merge.PE.tab]) == 0 or die qq[Failed to merge SR input files\n];
    }
    
    system( qq[sort -k4,4d -k1,1d -k2,2n $$.merge.PE.tab | uniq > $$.merge.PE.sorted.tab] ) == 0 or die qq[Failed to sort merged SR input files\n];
    
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
    
    open( my $tfh, qq[$$.merge.PE.sorted.tab] ) or die qq[Failed to open merged tab file: $!\n];
    my %typeBEDFiles;
    my $currentType = '';
    my $cfh;
    my $readCount = 0;
    my $raw_candidates = qq[$output.candidates];
    open( my $dfh, qq[>$raw_candidates] ) or die $!;
    print $dfh qq[FILTER: chr\tstart\tend\ttype\_sample\tL_Fwd_both\tL_Rev_both\tL_Fwd_single\tL_Rev_single\tR_Fwd_both\tR_Rev_both\tR_Fwd_single\tR_Rev_single\tL_Last_both\tR_First_both\tDist\n];
    close( $dfh );
    while( my $line = <$tfh> )
    {
        chomp( $line );
        my @s = split( /\t/, $line );
        if( $currentType eq '' )
        {
            $currentType = $s[ 3 ];
            print qq[PE Call: $currentType\n];
            open( $cfh, qq[>$$.$currentType.pe_anchors.bed] ) or die $!;
        }
        elsif( $currentType ne $s[ 3 ] || eof( $tfh ) )
        {
            close( $cfh );
            
            if( $currentType ne 'unknown' || ($novel && $currentType eq 'unknown') )
            {
                #convert to a region BED (removing any candidates with very low numbers of reads)
                print qq[$currentType Calling initial rough boundaries of insertions....\n];
                my $rawTECalls1 = qq[$$.raw_calls.1.$currentType.tab];
                my $numRegions = Utilities::convertToRegionBedPairs( qq[$$.$currentType.pe_anchors.bed], $DEFAULT_MIN_CLUSTER_READS, $currentType, $MAX_READ_GAP_IN_REGION, $DEFAULT_MAX_CLUSTER_DIST, 1, 0, $rawTECalls1 );
                
                if( $numRegions == 0 )
                {
                    unlink( qq[$$.$currentType.pe_anchors.bed], $rawTECalls1 ) or die qq[Failed to delete temp files\n] if( $clean );
                    
                    $readCount = 0;
                    $currentType = $s[ 3 ];
                    print qq[PE Call: $currentType\n];
                    open( $cfh, qq[>$$.$currentType.pe_anchors.bed] ) or die $!;
                    next;
                }
                
                if( defined( $filterBED ) )
                {
                    #remove the regions specified in the exclusion BED file
                    my $filtered = qq[$$.raw_calls.1.filtered.$currentType.tab];
                    Utilities::filterOutRegions( $rawTECalls1, $filterBED, $filtered );
                    $rawTECalls1 = $filtered;
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
                
                $typeBEDFiles{ $currentType }{hom} = $homCalls if( -f $homCalls && -s $homCalls > 0 );
                if( $hets && -f $hetCalls && -s $hetCalls > 0 ){$typeBEDFiles{ $currentType }{het} = $hetCalls;}
                
                #remove the temporary files for this 
                if( $clean )
                {
                    unlink( qq[$$.raw_reads.0.$currentType.tab], qq[$$.raw_calls.1.$currentType.tab], qq[$$.raw_calls.2.$currentType.tab] ) or die qq[Failed to delete temp files\n];
                    unlink( $homCalls ) if( -f $homCalls && -s $homCalls == 0 );
                    unlink( $hetCalls ) if ( -f $hetCalls && -s $hetCalls == 0 );
                }
            }
            
            last if( eof( $tfh ) );
            
            $readCount = 0;
            $currentType = $s[ 3 ];
            print qq[PE Call: $currentType\n];
            open( $cfh, qq[>$$.$currentType.pe_anchors.bed] ) or die $!;
        }
        else
        {
            #             chr    start  stop   TE_type   a_orien  m_orien %id
            print $cfh qq[$s[0]\t$s[1]\t$s[2]\t$currentType\t$s[5]\n];
            $readCount ++;
        }
    }
    close( $tfh );
    close( $cfh );
    
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
            
            if( $hets )
            {
                $f = $typeBEDFiles{$type}{het}; 
                system( qq[cat $f >> $$.PE.allcalls.bed] ) == 0 or die qq[Failed to cat file: $f] if( $f && -s $f );
            }
        }
        
        if( -f $unknownHoms )
        {
            my $removed = Utilities::filterOutRegions( $unknownHoms, qq[$$.PE.allcalls.bed], qq[$$.unknownHoms.filtered] );
            if( $removed > 0 )
            {
                print qq[Removed $removed hom unknown calls\n];
                $typeBEDFiles{unknown}{hom} = qq[$$.unknownHoms.filtered];
            }
        }
        
        if( $hets && defined( $unknownHets ) && -s $unknownHets )
        {
            my $removed = Utilities::filterOutRegions( $unknownHets, qq[$$.PE.allcalls.bed], qq[$$.unknownHets.filtered] );
            if( $removed > 0 )
            {
                print qq[Removed $removed het unknown calls\n];
                $typeBEDFiles{unknown}{het} = qq[$$.unknownHets.filtered];
            }
        }
    }
    
    #make a VCF file with the calls
    _outputCalls( \%typeBEDFiles, $sampleName, $ref, $output );
	
    #output calls in VCF and BED format
    print qq[Creating VCF file of calls....\n];
    _outputCalls( \%typeBEDFiles, $sampleName, $ref, $output );
    
	if( $clean )
	{
	    #delete the intermediate files
	    unlink( glob( qq[$$.*] ) ) or die qq[Failed to remove intermediate files: $!];
	}
}

sub _findInsertionsSR
{
    my $bamsRef = shift;my @bams = @{$bamsRef};
    my $sample = shift;
    my $input = shift;
    my $ref = shift;
    my $output = shift;
    my $minReads = shift;
    my $depth = shift;
    my $region = shift;
    my $clean = shift;
    my $hets = shift;
    my $ignoreRGsFormatted = shift;
    my $novel = shift;
    
    my @files;
    if( -f $input )
    {
        my $first = `head -1 $input`;chomp( $first );
        if( -f $first ) #is this a fofn
        {
            open( my $ifh, $input ) or die qq[failed to open fofn of SR discovery output files: $input\n];
            while(my $file = <$ifh> )
            {
                chomp( $file );
                if( -f $file ){push(@files, $file);}else{die qq[Cant find discovery output file: $file\n];}
            }
            print qq[Found ].scalar(@files).qq[ SR discovery stage input files\n\n];
            close( $ifh );
        }else{push( @files, $input );}#looks like a discovery output file
    }
    else
    {
        @files = glob( qq[$input*] ) or die qq[Failed to glob files: $input*\n];
        die qq[Cant find any inputs files with prefix $input*\n] if( ! @files || @files == 0 );
    }
    
    #merge the SR discovery files together and sort them by type, chr, position cols
    foreach my $file( @files )
    {
        system(qq[cat $file >> $$.merge.SR.tab]) == 0 or die qq[Failed to merge SR input files\n];
    }
    
    system( qq[sort -k4,4d -k1,1d -k2,2n $$.merge.SR.tab | uniq > $$.merge.SR.sorted.tab] ) == 0 or die qq[Failed to sort merged SR input files\n];
    
    #grab the reads that didnt match any TE into a separate file (to use as evidence when calling individual elements) - NEED to adjust the coords to the actual breakpoints first
    open( my $tfh, qq[$$.merge.SR.sorted.tab] ) or die qq[failed to open candidates file: $!\n];
    open( my $ufh, qq[>$$.merge.SR.unknown.tab] );
    while( my $line = <$tfh> )
    {
        chomp( $line );
        my @s = split( /\t/, $line );
        next unless $s[ 4 ] eq 'unknown';
        my $breakpoint = Utilities::calculateCigarBreakpoint( $s[ 1 ], $s[ 6 ] );
        my $orientation = ($s[ 5 ] & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
        print $ufh qq[$s[0]\t$breakpoint\t$breakpoint\tunknown\t$orientation\n];
    }
    close( $tfh );
    close( $ufh );
    
    open( $tfh, qq[$$.merge.SR.sorted.tab] ) or die qq[Failed to open merged tab file: $!\n];
    my %typeBEDFiles;
    my $currentType = '';
    my $cfh;
    my $readCount = 0;
    my %currentReads;
    while( my $line = <$tfh> )
    {
        my @s = split( /\t/, $line );#print qq[$line\n];
        if( $currentType eq '' )
        {
            $currentType = $s[ 3 ];
            print qq[SR Call: $currentType\n];
            open( $cfh, qq[>$$.$currentType.sr_anchors.bed] ) or die $!;
        }
        elsif( $currentType ne $s[ 3 ] || eof( $tfh ) )
        {
            close( $cfh );
            
            if( $currentType ne 'unknown' || ($novel && $currentType eq 'unknown') )
            {
                #resort the bed file by pos
                system( qq[sort -k1,1d -k2,2n $$.$currentType.sr_anchors.bed $$.merge.SR.unknown.tab > $$.$currentType.sr_anchors.sorted.bed] ) == 0 or die qq[Sort failed on $currentType\n];
                
                
                #call the regions
                my $minReads = $currentType eq 'unknown' ? $DEFAULT_SR_MIN_UNKNOWN_CLUSTER_READS : $DEFAULT_SR_MIN_CLUSTER_READS;
                if( $readCount > 0 && Utilities::convertToRegionBedPairs( qq[$$.$currentType.sr_anchors.sorted.bed], $minReads, $currentType, $MAX_READ_GAP_IN_REGION, $DEFAULT_MAX_SR_CLUSTER_DIST, 0, 1, qq[$$.$currentType.sr_calls.bed] ) > 0 )
                {
                    $typeBEDFiles{ $currentType }{hom} = qq[$$.$currentType.sr_calls.bed];
                    #unlink( qq[$$.$currentType.sr_anchors.bed], qq[$$.$currentType.sr_anchors.sorted.bed] );
                }else{unlink( qq[$$.$currentType.sr_anchors.bed], qq[$$.$currentType.sr_anchors.sorted.bed] );}
            }
            last if( eof( $tfh ) );
            
            %currentReads = ();
            $readCount = 0;
            $currentType = $s[ 3 ];
            open( $cfh, qq[>$$.$currentType.sr_anchors.bed] ) or die $!;
        }
        elsif( ! $currentReads{ $s[4] } )
        {
            #check it is paired correctly   AND there is only 1 breakpoint on the read
            if( ( $s[ 5 ] & $$BAMFLAGS{'read_paired'}) && ($s[6]=~tr/S/S/) == 1 )
            {
                my $orientation = ($s[ 5 ] & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
                
                #get the actual breakpoint from the cigar string, flag, and position field
                my $breakpoint = Utilities::calculateCigarBreakpoint( $s[ 1 ], $s[ 6 ] );
                
                print $cfh qq[$s[0]\t$breakpoint\t$breakpoint\t$currentType\t$orientation\n];
                $readCount ++;
                $currentReads{ $s[4] } = 1;
            }
        }
    }
    close( $tfh );
    close( $cfh );
    
    if( $novel )
    {
        #remove all the unknown calls that overlap with calls that are assigned to a type
        my $unknownHoms = $typeBEDFiles{unknown}{hom};
        foreach my $type( keys( %typeBEDFiles ) )
        {
            next if( $type eq 'unknown' );
            my $f = $typeBEDFiles{$type}{hom};
            system( qq[cat $f >> $$.allcalls.bed] ) == 0 or die qq[Failed to cat file: $f] if( $f && -s $f );
        }
        
        my $removed = Utilities::filterOutRegions( $unknownHoms, qq[$$.allcalls.bed], qq[$$.unknownHoms.filtered] );
        if( $removed > 0 )
        {
            print qq[Removed $removed hom unknown calls\n];
            $typeBEDFiles{unknown}{hom} = qq[$$.unknownHoms.filtered];
        }
        
        #test if the breakpoint looks like an inversion breakpoint
        foreach my $type( keys( %typeBEDFiles ) )
        {
            my $homs = $typeBEDFiles{$type}{hom};
            my $removed = Utilities::checkBreakpointsSR($homs, \@bams, ($DEFAULT_SR_MIN_CLUSTER_READS * 2), $ignoreRGsFormatted, $DEFAULT_ANCHORQ, qq[$homs.breakpoints_checked]);
            if( $removed > 0 ){$typeBEDFiles{$type}{hom} = qq[$homs.breakpoints_checked];}
        }
    }
    
    #make a VCF file with the calls
    _outputCalls( \%typeBEDFiles, $sample, $ref, $output );
    
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
	if( $hets )
	{
	    open( $hetsfh, qq[>$bedoutHets] ) or die $!;
	}
	
	open( my $dfh, qq[>>$thrown_out_file] ) or die $!;
	while( my $originalCall = <$ifh> )
	{
	    chomp( $originalCall );
	    print qq[Considering call $originalCall\n];
	    
	    my @originalCallA = split( /\t/, $originalCall );
	    my $strain = (split( /\./, $originalCallA[3] ))[ 0 ];
	    
	    my $start = $originalCallA[ 1 ];my $end = $originalCallA[ 2 ];my $chr = $originalCallA[ 0 ];
	    
	    my $t = Utilities::getCandidateBreakPointsDirVote($originalCallA[0], $originalCallA[1], $originalCallA[2],\@bams,20 );
	    if( ! $t ){warn qq[Failed to get candidate breakpoints for call $originalCall\n];next;}
	    my $breakpoint = $t->[ 0 ];my $depth = $t->[ 1 ];
	    
        print qq[Testing breakpoint $breakpoint\n];
        
        my $result = Utilities::testBreakPoint( $originalCallA[ 0 ], $breakpoint, \@bams, $minMapQ, $originalCall, $dfh, $ignoreRGsFormatted, $minReads, $DEFAULT_MIN_SOFT_CLIP, 0 );
        
        my $flag = $result->[0];
        my $call = $result->[1];
        my $ratio = $result->[2];
        
        if( $flag == $Utilities::INV_BREAKPOINT )
        {
            print qq[Failed - at inversion breakpoint\n];
        }
        elsif( $flag > $Utilities::NOT_ENOUGH_READS_CLUSTER )
        {
            #decide if its a hom or het call by the depth
            if( $depth <= 10 )
            {
                print qq[Hom breakpoint score: $ratio Flag: $flag\n];
                print $homsfh $call.qq[\t$flag\n];
            }
            elsif( $hets )
            {
                print qq[Het breakpoint score: $ratio Flag: $flag\n];
                print $hetsfh $call.qq[\t$flag\n];
            }
        }
        else{print qq[Discarding: $call : $flag\n];}
	}
	close( $ifh );
	close( $homsfh );
	close( $hetsfh ) if( $hets );
	close( $dfh );
}

=pod
Idea is to take a list of bams for new samples (e.g. low cov),
and look in the region around the VCF of calls to see if there is some
support for the call in the new sample and output a new VCF with the
new genotypes called
=cut
sub _genotype
{
    croak qq[Incorrect number of parameters to _gentoype function: ].scalar(@_) unless @_ == 10;
    my $bam_fofn = shift;
	my $input = shift;
	my $ref = shift;
    my $region = shift;
	my $minMapQ = shift;
	my $output = shift;
	my $clean = shift;
	my $hets = shift;
	my $orientate = shift;
	my $minReads = shift;
	
	#get the list of sample names
    my %sampleBAM;
    open( my $tfh, $bams ) or die $!;
    while(<$tfh>)
    {
        chomp;
        my $bam = $_;
        die qq[Cant find BAM file: $_\n] unless -f $bam;
        die qq[Cant find BAM index for BAM: $bam\n] unless -f qq[$bam.bai];
        
        my $bams = [$bam];
        my $s = Utilities::getBAMSampleName( $bams );
        if( $s ){$sampleBAM{ $s } = $bam;}else{die qq[Failed to determine sample name for BAM: $bam\nCheck SM tag in the read group entries.\n];exit;}
    }close( $tfh );
    
    #parse the region
    my ($chr, $start, $end) = (undef, undef, undef);
    if( defined( $region ) )
    {
        if( $region =~ /^([A-Za-z0-9]+):([0-9]+)-([0-9]+)$/ )
        {
            $chr = $1;$start=$2;$end=$3;
            print qq[Restricting genotyping to region: $region\n];
        }
        else{$chr = $region;print qq[Restricting genotyping to Chr: $chr\n];}
    }
    
    my $vcf = Vcf->new(file=>$input);
    $vcf->parse_header();
    my $vcf_out = Vcf->new();
    open( my $out, qq[>$output] ) or die $!;
    foreach my $sample ( keys( %sampleBAM ) )
	{
	    $vcf_out->add_columns( $sample );
	}
	my $header = Utilities::getVcfHeader( $vcf_out );
	print $out qq[$header];
	
	open( my $ofh, qq[>$output.candidates] ) or die $!;
    print $ofh qq[FILTER: chr\tstart\tend\ttype\_sample\tL_Fwd_both\tL_Rev_both\tL_Fwd_single\tL_Rev_single\tR_Fwd_both\tR_Rev_both\tR_Fwd_single\tR_Rev_single\tL_Last_both\tR_First_both\tDist\n];
    while( my $entry = $vcf->next_data_hash() )
    {
        my $chr_ = $$entry{CHROM};
        my $pos = $$entry{POS};
        
        if( defined( $chr ) )
        {
            if( $chr_ ne $chr )
            {
                next;
            }
            elsif( defined( $start ) && defined( $end ) )
            {
                if( $pos < $start ){next;}
                elsif( $pos > $end ){last;}
            }
        }
        
        print qq[Genotyping call: $chr_\t$pos\n];
        
        #get the existing GT call
        my @samples = keys( %{$$entry{gtypes}});
        my $gt = $$entry{gtypes}{$samples[ 0 ]}{GT};
        my $type = $$entry{INFO}{MEINFO};my @typeInfo = split( /,/, $type );
        
        foreach my $sample ( sort( keys( %sampleBAM ) ) )
        {
            my @bams = $sampleBAM{ $sample };
            my $quality = Utilities::testBreakPoint($chr_, $pos, \@bams, $minMapQ, qq[$chr_\t$pos\t].($pos+1).qq[\t].$typeInfo[0].qq[\t].$typeInfo[3].qq[\n], $ofh, undef, $minReads, 1 );
            if( $quality )
	        {
	            $$entry{gtypes}{$sample}{GT} = $gt;
                $$entry{gtypes}{$sample}{GQ} = $quality->[0];
	        }
	        else
	        {
	            $$entry{gtypes}{$sample}{GT} = qq[./.];
                $$entry{gtypes}{$sample}{GQ} = '.';
	        }
        }
        $vcf->format_genotype_strings($entry);
        print $out $vcf_out->format_line($entry);
    }
    close( $ofh );
    $vcf->close();
	close( $out );
}

#***************************INTERNAL HELPER FUNCTIONS********************

sub _getCandidateTEReadNames
{
    my $bam = shift;
    my $readgroups = shift;
    my $minAnchor = shift;
    my $minSoftClip = shift;
    my $candidatesFasta = shift;
    my $candidatesBed = shift;
    my $clippedFasta = shift;
    
    my %candidates;
    open( my $ffh, qq[>$candidatesFasta] ) or die qq[ERROR: Failed to create fasta file: $!\n];
    open( my $afh, qq[>$candidatesBed] ) or die qq[ERROR: Failed to create anchors file: $!\n];
    my $cfh;
    if( $minSoftClip )
    {
        open( $cfh, qq[>$clippedFasta] ) or die qq[ERROR: Failed to create clip fasta file: $!\n];
    }
    
    print qq[Opening BAM ($bam) and getting initial set of candidate mates....\n];
    
    open( my $bfh, qq[samtools view ].( defined( $readgroups ) ? qq[-R $$.readgroups ] : qq[ ] ).qq[$bam |] ) or die $!;
    my $currentChr = '';
    while ( my $samLine = <$bfh> )
    {
        chomp( $samLine );
        my @sam = split( /\t/, $samLine );
        my $flag = $sam[ 1 ];;
        my $qual = $sam[ 4 ];
        my $name = $sam[ 0 ];
        my $ref = $sam[ 2 ];
        my $mref = $sam[ 6 ];
        my $readLen = length( $sam[ 9 ] );
        my $cigar = $sam[ 5 ];
        
        if( $candidates{ $name } )
        {
            if( ($flag & $$BAMFLAGS{'1st_in_pair'}) && $candidates{ $name } == 1 || ($flag & $$BAMFLAGS{'2nd_in_pair'}) && $candidates{ $name } == 2 )
            {
                my $seq = $sam[ 9 ];
                
                if( _sequenceQualityOK( $seq ) )
                {
                    print $ffh qq[>$name\n$seq\n];
                }
                delete( $candidates{ $name } );
            }
        }
        
        my $supporting = Utilities::isSupportingClusterRead( $flag, $sam[ 8 ], $qual, $minAnchor, $minSoftClip, $cigar );
        
        if( $supporting > 0 )
        {
            if( $supporting == 1 )
            {
               $candidates{ $name } = $flag & $$BAMFLAGS{'1st_in_pair'} ? 2 : 1; #mate is recorded
                   
               my $pos = $sam[ 3 ];
               my $dir = ($flag & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
               print $afh qq[$ref\t$pos\t].($pos+$readLen).qq[\t$name\t$dir\n];
            }
            elsif( $supporting == 2 )
            {
               $candidates{ $name } = $flag & $$BAMFLAGS{'1st_in_pair'} ? 2 : 1; #mate is recorded
                   
               my $pos = $sam[ 3 ];
               my $dir = ($flag & $$BAMFLAGS{'reverse_strand'}) ? '-' : '+';
               print $afh qq[$ref\t$pos\t].($pos+$readLen).qq[\t$name\t$dir\n];
            }
            
            #if in SR mode - then see if it is a candidate split read
            #                               read mapped                             mate mapped                            mate on same chr             ins size sensible          has single soft clip in cigar bigger than minimum clip size
            if( $minSoftClip && ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ( $mref eq '=' || $mref eq $ref ) && $sam[ 8 ] < 3000 && ($cigar=~tr/S/S/) == 1 && $cigar =~ /(\d+)(S)/ && $1 > $minSoftClip )        #( $flag & $$BAMFLAGS{'read_paired'} ) 
            {
                my $seq = $sam[ 9 ];
                my $quals = $sam[ 10 ];
                
                #check there arent high numbers of mismatches in the aligned portion (via MD tag if available)
                if( $samLine =~ /\tMD:Z:([0-9ACGT\^]+)/ )
                {
                    my $md = $1;
                    my $mismatches = ( $md =~ tr/ACGT/n/ );
                    if( $mismatches <= 5 )
                    {
                        #clip at start of alignment
                        if( $cigar =~ /^(\d+)S/ )
                        {
                            my $clip = $1;
                            if( ( $flag & $$BAMFLAGS{'reverse_strand'} ) )
                            {
                                if( _checkAvgQuality( substr($quals, length($seq)-$clip) ) ){print $cfh qq[>$name###$sam[2]###$sam[3]###$flag###$cigar###$sam[8]\n];print $cfh substr($seq, length($seq)-$clip).qq[\n];}
                            }
                            else
                            {
                                if( _checkAvgQuality( substr($quals, 0, $clip ) ) ){print $cfh qq[>$name###$sam[2]###$sam[3]###$flag###$cigar###$sam[8]\n];print $cfh substr($seq, 0, $clip ).qq[\n];}
                            }
                        }
                        elsif( $cigar =~ /(\d+)S$/ )
                        {
                            my $clip = $1;
                            if( ( $flag & $$BAMFLAGS{'reverse_strand'} ) )
                            {
                                if( _checkAvgQuality( substr($quals, 0, $clip ) ) ){print $cfh qq[>$name###$sam[2]###$sam[3]###$flag###$cigar###$sam[8]\n];print $cfh substr($seq, 0, $clip ).qq[\n];}
                            }
                            else
                            {
                                if( _checkAvgQuality( substr($quals, length($seq)-$clip) ) ){print $cfh qq[>$name###$sam[2]###$sam[3]###$flag###$cigar###$sam[8]\n];print $cfh substr($seq, length($seq)-$clip).qq[\n];}
                            }
                        }
                    }
                }
            }
        }
        if( $currentChr ne $ref && $ref ne '*' ){print qq[Reading chromosome: $ref\n];$currentChr = $ref;}
    }
    close( $bfh );
    close( $ffh );
    close( $afh );
    close( $cfh );
    
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
    
    my $vcf_out = Vcf->new();
    $vcf_out->add_columns($sample);
    
    my $header = Utilities::getVcfHeader($vcf_out);
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
                my $dir = defined( $s[ 6 ] ) && length( $s[ 6 ] ) > 0 ? $s[ 6 ] : 'NA';
                
                print $bfh qq[$s[0]\t$pos\t].($pos+1).qq[\t$type==$sample\t$s[4]\t$dir\n];
                
                my %out;
                $out{CHROM}  = $s[ 0 ];
                $out{POS}    = $pos;
                $out{ID}     = '.';
                $out{ALT}    = ['<INS:ME>'];
                $out{REF}    = $refbase;
                $out{QUAL}   = $s[ 4 ];
                $out{INFO} = { SVTYPE=>'INS', NOT_VALIDATED=>undef, MEINFO=>qq[$type,$s[1],$s[2],$dir] };
                $out{FORMAT} = ['GT', 'GQ', 'FL'];
                
                if( $flag == $Utilities::PASS )
                {
                    $out{gtypes}{$sample}{GT} = qq[<INS:ME>/<INS:ME>];
                }
                else
                {
                    $out{gtypes}{$sample}{GT} = qq[0/0];
                }
                $out{gtypes}{$sample}{GQ} = qq[$s[4]];
                $out{gtypes}{$sample}{FL} = $flag;
                
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
                my $dir = defined( $s[ 6 ] ) && length( $s[ 6 ] ) > 0 ? $s[ 6 ] : 'NA';
                
                print $bfh qq[$s[0]\t$pos\t].($pos+1).qq[\t$type==$sample\t$s[4]\t$dir\n];
                
                my %out;
                $out{CHROM}  = $s[ 0 ];
                $out{POS}    = $pos;
                $out{ID}     = '.';
                $out{ALT}    = [$refbase,'<INS:ME>'];
                $out{REF}    = $refbase;
                $out{QUAL}   = $s[ 4 ];
                $out{INFO} = { SVTYPE=>'INS', NOT_VALIDATED=>undef, MEINFO=>qq[$type,$s[1],$s[2],$dir] };
                $out{FORMAT} = ['GT','GQ', 'FL'];
                
                if( $flag == $Utilities::PASS )
                {
                    $out{gtypes}{$sample}{GT} = qq[$refbase/<INS:ME>];
                }
                else
                {
                    $out{gtypes}{$sample}{GT} = qq[0/0];
                }
                $out{gtypes}{$sample}{GQ} = qq[$s[4]];
                $out{gtypes}{$sample}{FL} = $flag;
                
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

sub _sortBED
{
    my $input = shift;
    my $output = shift;
    
    croak qq[Cant find intput file for BED sort: $input\n] unless -f $input;
    system( qq[sort -k 1,1d -k 2,2n $input > $output] ) == 0 or die qq[Failed to sort BED file: $input\n];
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
        die qq[Tab file should have entries separated by single tab: $file\n] unless $entry =~ /^(.+)(\t)(.+)$/;
        $hash{ $1 } = $3;
    }
    close( $tfh );
    
    return \%hash;
}

sub _filterBED
{
    my $bed = shift;
    my $chr = shift;
    my $start;my $end;
    if( @_ == 2 )
    {
        $start = shift;
        $end = shift;
    }
    
    open( my $ofh, qq[>$$.filter.temp] ) or die $!;
    open( my $ifh, $bed ) or die $!;
    while( my $line = <$ifh> )
    {
        chomp( $line );
        my @s = split( /\t/, $line );
        if( defined($start) && defined($end) )
        {
            next unless $s[ 0 ] eq $chr && $s[ 1 ] > $start && $s[ 2 ] < $end;
        }
        else{next unless $line =~ /^$chr\t/;}
        print $ofh qq[$line\n];
    }
    close( $ifh );close( $ofh );
    
    unlink( $bed ) or die qq[Failed to remove sorted BED: $!];
    rename(qq[$$.filter.temp], $bed) or die qq[Failed to rename filtered BED: $!];
}

sub _mergeDiscoveryOutputs
{
    my $t = shift;
    my @files = @{ $t };
    my $output = shift;
    
    my %typeFile;
    my $typeCount = 0;
    my $currentFh;
    my $header;
    
    #iterate through the files
    foreach my $file ( @files )
    {
        die qq[Cant find candidates file: $file\n] unless -f $file;
        
        #check its a valid discovery output file - not truncated etc.
        _checkDiscoveryOutput( $file );

        open( my $tfh, $file ) or die $!;
        $header = <$tfh>; chomp( $header );
        while( my $line = <$tfh> )
        {
            chomp( $line );
            next if( $line =~ /^#/ || $line =~ /^TE_TYPE_END/ );
            if( $line =~ /^(TE_TYPE_START)(\s+)(.+)$/ )
            {
                my $currentType = $3;
                                
                if( ! $typeFile{ $currentType } )
                {
                    $typeFile{ $currentType } = qq[$$.merged.$typeCount];
                    $typeCount ++;
                }
                
                if( $currentFh ){close( $currentFh );}
                open( $currentFh, qq[>>].$typeFile{ $currentType } ) or die $!;
            }
            else
            {
                print $currentFh qq[$line\n];
            }
        }
        close( $tfh );
    }
    if( $currentFh ){close( $currentFh );}
    
    open( my $ofh, qq[>$output] ) or die $!;
    #now merge into a single file
    print $ofh $HEADER.qq[\n];
    foreach my $type ( keys( %typeFile ) )
    {
        print $ofh qq[TE_TYPE_START $type\n];
        open( my $ifh, $typeFile{ $type } ) or die $!;
        while( <$ifh> )
        {
            chomp;
            print $ofh qq[$_\n];
        }
        close( $ifh );
        print $ofh qq[TE_TYPE_END\n];
        
        #delete all of the intermediate files
        unlink( $typeFile{ $type } ) or die qq[Failed to delete intermediate file\n];
    }
    print $ofh $FOOTER.qq[\n];
    close( $ofh );
}

