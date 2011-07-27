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
The BAM could be required to be sorted by readname before input
    - this would mean that you could get the candidate reads from a single pass through
    - drawback is that for further stages you need the BAM sorted by coordinate
    - for big projects, it would be extra effort to pre-sort by name
Calling of heterozygotes is not handled right now
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

my $VERSION = 0.2;

my $DEFAULT_ID = 90;
my $DEFAULT_LENGTH = 36;
my $DEFAULT_ANCHORQ = 20;
my $DEFAULT_MAX_DEPTH = 200;
my $DEFAULT_READS = 5;
my $DEFAULT_MIN_GENOTYPE_READS = 3;
my $MAX_READ_GAP_IN_REGION = 50;

my $HEADER = qq[#retroseq v:].substr($VERSION,0,1);
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

my ($discover, $call, $genotype, $bam, $bams, $ref, $eRefFofn, $length, $id, $output, $anchorQ, $region, $input, $reads, $depth, $noclean, $tmpdir, $readgroups, $filterFile, $heterozygous, $orientate, $ignoreRGsFofn, $help);

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
    'ignoreRGs=s'     =>  \$ignoreRGsFofn,
    'h|help'        =>  \$help,
);

print <<MESSAGE;

RetroSeq: A tool for discovery and genotyping of transposable elements from short read alignments
Version: $VERSION
Author: Thomas Keane (thomas.keane\@sanger.ac.uk)

MESSAGE

my $USAGE = <<USAGE;
Usage: $0 -<command> options

            -discover       Takes a BAM and a set of reference TE (fasta) and calls candidate supporting read pairs (BED output)
            -call           Takes multiple output of discovery stage and a BAM and outputs a VCF of TE calls
            -genotype       Input is a VCF of TE calls and a set of sample BAMs, output is a new VCF with genotype calls for new samples
            
NOTE: $0 requires samtools, bcftools, exonerate, unix sort to be in the default path

USAGE

( $discover || $call || $help) or die $USAGE;

if( $discover )
{
    ( $bam && $eRefFofn && $output ) or die <<USAGE;
Usage: $0 -discover -bam <string> -eref <string> -output <string> [-q <int>] [-id <int>] [-len <int> -clean <yes/no>]
    
    -bam        BAM file of paired reads mapped to reference genome
    -eref       Tab file with list of transposon types and the corresponding fasta file of reference sequences (e.g. SINE   /home/me/refs/SINE.fasta)
    -output     Output file to store candidate supporting reads (required for calling step)
    [-noclean   Do not remove intermediate output files. Default is to cleanup.]
    [-q         Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-id        Minmum percent ID for a match of a read to the transposon references. Default is 90.]
    [-len       Miniumum length of a hit to the transposon references. Default is 36bp.]
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
    
    if( $readgroups && length( $readgroups ) > 0 )
    {
        my @s = split( /,/, $readgroups );foreach my $rg ( @s ){if( $rg !~ /[A-Za-z0-9]|\.|-|_+/ ){croak qq[Invalid readgroup: $rg\n];}}
        print qq[Restricting discovery phase to read groups: $readgroups\n];
    }
    
    print qq[\nMin anchor quality: $anchorQ\nMin percent identity: $id\nMin length for hit: $length\n\n];
    
    #test for samtools
    Utilities::checkBinary( q[samtools], qq[0.1.16] );
    Utilities::checkBinary( q[exonerate], qq[2.2.0] );
    
    _findCandidates( $bam, $erefs, $id, $length, $anchorQ, $output, $readgroups, $clean );
}
elsif( $call )
{
    ( $bam && $input && $ref && $output ) or die <<USAGE;
Usage: $0 -call -bam <string> -input <string> -ref <string> -output <string> [-filter <BED file> -cleanup -reads <int> -depth <int> -hets]
    
    -bam            BAM file OR BAM fofn
    -input          Either a single output file from the discover stage OR a prefix of a set of files from discovery to be combined for calling OR a fofn of discovery stage output files
    -ref            Fasta of reference genome
    -output         Output file name (VCF)
    [-hets          Call heterozygous insertions. Default is homozygous.]
    [-orientate     Attempt to predicate the orientation of the calls. Default is no.]
    [-filter        BED file of regions to exclude from final calls e.g. reference elements, low complexity etc.]
    [-region        Call a particular chromosome only (chr) OR region (chr:start-end) only]
    [-depth         Max average depth of a region to be considered for calling. Default is 200.]
    [-reads         It is the minimum number of reads required to make a call. Default is 5. Parameter is optional.]
    [-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-ignoreRGs     Single read group name OR a file of readgroups that should be ignored. Default is none.]
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
    
    _findInsertions( $bam, $input, $ref, $output, $reads, $depth, $anchorQ, $region, $clean, $filterFile, $heterozygous, $orientate, $ignoreRGsFofn );
}
elsif( $genotype )
{
    ( $bams && $input && $reads && $noclean && $tmpdir && $output ) or die <<USAGE;
Usage: $0 -genotype -bams <string> -input <string> -eref <string> -output <string> [-cleanup -reads <int> -region <chr:[start-end]>
    
    -bams           File of BAM file names (one per sample to be genotyped)
    -input          VCF file of TE calls
    -ref            Fasta of TE reference genome
    -output         Output VCF file (will be annotated with new genotype calls)
    [-hets          Call heterozygous insertions. Default is homozygous.]
    [-orientate     Attempt to predicate the orientation of the calls. Default is no.]
    [-region        Call a particular chromosome only (chr) OR region (chr:start-end) only]
    [-depth         Max average depth of a region to be considered for calling. Default is 200.]
    [-reads         It is the minimum number of reads required to make a call. Default is 5. Parameter is optional.]
    [-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-noclean       Do not remove intermediate output files. Default is to cleanup.]
USAGE
    
    croak qq[Cant find BAM fofn: $bams] unless -f $bams;
    croak qq[Cant find input: $input] unless -f $input;
    
    my $clean = defined( $noclean ) ? 0 : 1;
    $reads = defined( $reads ) && $reads =~ /^\d+$/ ? $reads > -1 : $DEFAULT_MIN_GENOTYPE_READS;
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
    my %candidates = %{ _getCandidateTEReadNames($bam, $readgroups, undef, undef, undef, $minAnchor, $candidatesFasta, $candidatesBed ) };
    
    print scalar( keys( %candidates ) ).qq[ candidate reads remain to be found after first pass....\n];
    
    if( keys( %candidates ) > 0 )
    {
        open( my $ffh, qq[>>$candidatesFasta] ) or die qq[ERROR: Failed to create fasta file: $!\n];
        open( my $afh, qq[>>$candidatesBed] ) or die qq[ERROR: Failed to create anchors file: $!\n];
        
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
            if( $currentChr ne $ref ){print qq[Reading chromosome: $ref\n];$currentChr = $ref;}
        }
        close( $ffh );
        close( $afh );
    }
    undef %candidates; #finished with this hash
    
    #create a single fasta file with all the TE ref seqs
    my %refLabels;
    my $counter = 1;
    my $refsFasta = qq[$$.allrefs.fasta];
    open( my $tfh, qq[>$refsFasta] ) or die $!;
    foreach my $type ( keys( %{ $erefs } ) )
    {
        open( my $sfh, $$erefs{ $type } ) or die $!;
        my $seqCount = 1;
        $refLabels{ $type } = $counter;
        while( my $line = <$sfh>)
        {
            chomp( $line );
            if( $line =~ /^>/ )
            {
                print $tfh qq[>$counter\_$seqCount\n];
                $seqCount ++;
            }
            else
            {
                print $tfh qq[$line\n];
            }
        }
        close( $sfh );
        $counter ++;
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
    #system( qq[exonerate --bestn 5 --percent ].($id-10).q[ --ryo "INFO: %qi %qal %pi %tS %ti\n"].qq[ $$.candidates.fasta $refsFasta | egrep "^INFO|completed" > $$.exonerate.out ] ) == 0 or die qq[Exonerate exited incorrectly\n];
    
    #open( my $efh, qq[$$.exonerate.out] ) or die qq[Failed to read exonerate alignment files: $!\n];
    open( my $efh, qq[exonerate --bestn 5 --percent ].($id-10).q[ --ryo "INFO: %qi %qal %pi %tS %ti\n"].qq[ $$.candidates.fasta $refsFasta | egrep "^INFO|completed" | ] ) or die qq[Exonerate failed to run: $!];
    print qq[Parsing alignments....\n];
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
            $s[ 5 ] =~ /(\d+)_(\d+)/;
            $anchors{ $s[ 1 ] } = [ $1, $s[ 4 ] ]; #this could be a memory issue (possibly dump out to a file and then use unix join to intersect with the bed)
        }
    }
    close( $efh );
    
    if( $lastLine ne qq[-- completed exonerate analysis] ){die qq[Alignment did not complete correctly\n];}
    
    #unlink( qq[$$.exonerate.out] );
    
    #setup filehandles for the various TE types
    my %TE_fh;
    foreach my $type ( keys( %{ $erefs } ) )
    {
        my $id = $refLabels{ $type };
        open( my $tfh, qq[>$$.$id.candidates] ) or die $!;
        $TE_fh{ $refLabels{ $type } } = $tfh; #keys are the IDs of the types
    }
    
    #now run through the anchors file and pull out the locations of the anchors
    open( my $cfh, qq[$$.candidate_anchors.bed] ) or die $!;
    while( my $anchor = <$cfh> )
    {
        chomp( $anchor );
        my @s = split( /\t/, $anchor );
        if( defined( $anchors{ $s[ 3 ] } ) )
        {
            my $fh = $TE_fh{ $anchors{ $s[ 3 ] }[ 0 ] }; #get the ID of the type the read was aligned to
            my $mate_orientation = $anchors{ $s[ 3 ] }[ 1 ];
            print $fh qq[$anchor\t$mate_orientation\n];
        }
    }
    close( $cfh );
    
    foreach my $type ( keys( %TE_fh ) )
    {
        close( $TE_fh{ $type } );
    }
    
    #now put all the anchors together into a single file per type
    open( my $afh, qq[>$output] ) or die $!;
	print $afh qq[$HEADER\n];
    foreach my $type ( keys( %{ $erefs } ) )
    {
        my $id = $refLabels{ $type };
        open( my $tfh, qq[$$.$id.candidates] ) or die $!;
        print $afh qq[TE_TYPE_START $type\n];
        while( my $line = <$tfh> )
        {
            print $afh qq[$line];
        }
        print $afh qq[TE_TYPE_END $type\n];
    }
    print $afh qq[$FOOTER\n]; #write an end of file marker
	close( $afh );
	
	if( $clean )
	{
	    #delete the intermediate files
	    unlink( glob( qq[$$.*] ) ) or die qq[Failed to remove intermediate files: $!];
	}
}

sub _findInsertions
{
    my $bam = shift;
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
            print qq[Found ].scalar(@files).qq[ discovery stage input files\n\n];
            close( $ifh );
        }else{push( @files, $input );}#looks like a discovery output file
    }
    else
    {
        @files = glob( qq[$input*] );
        die qq[Cant find any inputs files with prefix $input*\n];
    }
    
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
    
    my $input_;
    if( @files && @files > 1 )
    {
        $input_ = qq[$$.merged.discovery];
        _mergeDiscoveryOutputs( \@files, $input_ );
    }
    else
    {
        #check the eof markers are there from the discovery stage
        _checkDiscoveryOutput( $input );
        $input_ = $input;
    }
    
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
    
    my $sampleName = Utilities::getBAMSampleName( \@bams );
    print qq[Calling sample $sampleName\n];
    
    #parse the region
    my ($chr, $start, $end) = (undef, undef, undef);
    if( defined( $region ) )
    {
        if( $region =~ /^([A-Za-z0-9]+):([0-9]+)-([0-9]+)$/ )
        {
            $chr = $1;$start=$2;$end=$3;
            print qq[Restricting caling to region: $region\n];
        }
        else{$chr = $region;print qq[Restricting calling to Chr: $chr\n];}
    }
    
    #for each type in the file - call the insertions
    open( my $ifh, $input_ ) or die $!;
    my $currentType = '';
    my $tempUnsorted = qq[$$.reads.0.tab]; #a temporary file to dump out the reads for this TE type
    my %typeBEDFiles;
    my $count = 0;
    my $tfh;
    my $raw_candidates = qq[$output.candidates];
    open( my $dfh, qq[>$raw_candidates] ) or die $!;
    print $dfh qq[FILTER: chr\tstart\tend\ttype\_sample\tL_Fwd_both\tL_Rev_both\tL_Fwd_single\tL_Rev_single\tR_Fwd_both\tR_Rev_both\tR_Fwd_single\tR_Rev_single\tL_Last_both\tR_First_both\tDist\n];
    close( $dfh );
    my %sortedCandidatesBED;
    while(1)
    {
        my $line = <$ifh>;
        last unless defined( $line );
        chomp( $line );
        
        next if( $line =~ /^#/ );
        if( $line =~ /^(TE_TYPE_START)(\s+)(.+)$/ )
        {
            $currentType = $3;
            open( $tfh, qq[>$tempUnsorted] ) or die $!;
        }
        elsif( $line =~ /^TE_TYPE_END/ )
        {
            close( $tfh );
                print qq[Calling TE type: $currentType\n];
                if( defined( $chr ) )
                {
                    if( defined( $start ) && defined( $end ) )
                    {
                        #filter the BED file by the chromosome only
                        _filterBED( $tempUnsorted, $chr, $start, $end );
                    }
                    else
                    {
                        #filter the BED file by the chromosome only
                        _filterBED( $tempUnsorted, $chr );
                    }
                }
                
                #call the insertions
                $sortedCandidatesBED{$currentType} = qq[$$.raw_reads.0.$count.tab];
                _sortBED( $tempUnsorted, $sortedCandidatesBED{$currentType} );
                
                #convert to a region BED (removing any candidates with very low numbers of reads)
                print qq[Calling initial rough boundaries of insertions....\n];
                my $rawTECalls1 = qq[$$.raw_calls.1.$count.tab];
                Utilities::convertToRegionBED( $sortedCandidatesBED{$currentType}, $minReads, $sampleName, $MAX_READ_GAP_IN_REGION, $rawTECalls1 );
                
                if( defined( $filterBED ) )
                {
                    #remove the regions specified in the exclusion BED file
                    my $filtered = qq[$$.raw_calls.1.filtered.$count.tab];
                    Utilities::filterOutRegions( $rawTECalls1, $filterBED, $filtered );
                    $rawTECalls1 = $filtered;
                }
                
                #remove extreme depth calls
                print qq[Removing calls with extremely high depth (>$depth)....\n];
                my $rawTECalls2 = qq[$$.raw_calls.2.$count.tab];
                _removeExtremeDepthCalls( $rawTECalls1, \@bams, $depth, $rawTECalls2, $raw_candidates );
                
                #new calling filtering code
                print qq[Filtering and refining candidate regions into calls....\n];
                $typeBEDFiles{ $currentType }{hom} = qq[$$.raw_calls.3.$count.hom.bed]; #homozygous calls
                $typeBEDFiles{ $currentType }{het} = undef; #het calls
                if( $hets ){$typeBEDFiles{ $currentType }{het} = qq[$$.raw_calls.3.$count.het.bed];}
                _filterCallsBedMinima( $rawTECalls2, \@bams, 10, $minQ, $ref, $raw_candidates, $hets, $typeBEDFiles{ $currentType }{hom}, $typeBEDFiles{ $currentType }{het}, $ignoreRGsFormatted );
                $count ++;
        }
        else
        {
            print $tfh qq[$line\n];
        }
    }
    
    if( $orientate )
    {
        #add predicted orientation for the calls
        print qq[Orientating calls...\n];
        foreach my $type ( keys( %typeBEDFiles ) )
        {
            if( $typeBEDFiles{ $type }{hom} && -f $typeBEDFiles{ $type }{hom} )
            {
                Utilities::annotateCallsBED( $typeBEDFiles{ $type }{hom}, $sortedCandidatesBED{$type} );
            }
            
            if( $hets && $typeBEDFiles{ $type }{het} && -f $typeBEDFiles{ $type }{het} )
            {
                Utilities::annotateCallsBED( $typeBEDFiles{ $type }{het}, $sortedCandidatesBED{$type} );
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
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 10;
    
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
	    
	    my $t = Utilities::getCandidateBreakPointsDepth( $originalCallA[0], $originalCallA[1], $originalCallA[2], \@bams, $ref, $dfh );
	    if( ! $t ){warn qq[Failed to get candidate breakpoints for call $originalCall\n];next;}
	    my @t = @{ $t };
	    
	    my @positions = @{ $t[ 0 ] };my %min = %{ $t[ 1 ] };
	    
	    #test each point to see if has the desired signature of fwd / rev pointing reads
	    my $found = 0;my $tested = 0;
	    my $lastRefIndex = -1;
	    my $minRatio = 100000;my $minRatioCall = undef;
	    while( $tested < 5 )
	    {
	        #check the distance to the set tested so far (i.e. dont want to retest with a cluster of local minima)
	        my $newIndex = $lastRefIndex + 1;
	        while( $newIndex < @positions )
	        {
	            my $closeby = 0;
	            for( my $j = 0; $j < $newIndex; $j ++ )
	            {
	                if( abs( $positions[$newIndex] - $positions[ $j ] ) < 50 ){$closeby = 1;last;}
	            }
	            last if( $closeby == 0 );
	            $newIndex ++;
	        }
	        last if $newIndex == @positions;
	        
	        print qq[Testing hom breakpoint $positions[ $newIndex ]\n];
	        $lastRefIndex = $newIndex;
	        my $depth = $min{$positions[$newIndex]};
	        my $refPos = $positions[ $newIndex ];
	        
	        last unless $depth <= $minDepth;
	        
	        my $result = Utilities::testBreakPoint( $originalCallA[ 0 ], $refPos, \@bams, $minMapQ, $originalCall, $dfh, $ignoreRGsFormatted );
	        
	        if( $result && $result->[0] < $minRatio )
	        {
	            $minRatio = $result->[0];
	            $minRatioCall = $result->[1];
	        }
	        
	        $tested ++;
	    }
	    
	    if( $minRatioCall )
	    {
	        print $homsfh $minRatioCall;
	    }
	    elsif( $hets ) #if we are also testing for het calls - then try alternative method to call a het
	    {
	        my $candidateBreaks = Utilities::getCandidateBreakPointsDirVote( $originalCallA[ 0 ], $originalCallA[ 1 ], $originalCallA[ 2 ], \@bams, $minMapQ );
	        next if( !defined( $candidateBreaks ) );
	        
	        my $tested = 0;
	        foreach my $candPos( @{$candidateBreaks} )
	        {
	            last if( ! defined $candPos );
	            print qq[Testing het breakpoint: $candPos\n];
	            my $result = Utilities::testBreakPoint( $originalCallA[ 0 ], $candPos, \@bams, $minMapQ, $originalCall, $dfh, $ignoreRGsFormatted );
	            if( $result )
	            {
	                $minRatio = $result->[0];
	                $minRatioCall = $result->[1];
	                last;
	            }
	            $tested ++;last if $tested > 5;
	        }
	        
	        if( $minRatioCall )
	        {
	            print $hetsfh $minRatioCall;
	        }
	    }
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
    
    my $vcf = Vcf->new(file=>$input);
    $vcf->parse_header();
    my $vcf_out = Vcf->new();
    open( my $out, qq[>$output] ) or die $!;
    foreach my $sample ( keys( %sampleBAM ) )
	{
	    $vcf_out->add_columns( $sample );
	}
	_writeVcfHeader( $vcf_out, $out );
	
    while( my $entry = $vcf->next_data_hash() )
    {
        my $chr_ = $$entry{CHROM};
        my @ci = split( /,/, $$entry{INFO}{CIPOS} );
        my $pos = $$entry{POS};
        my $start = ($$entry{POS} + $ci[ 0 ]) > 1 ? $$entry{POS} + $ci[ 0 ] : 1;
        my $end = $$entry{POS} + $ci[ 1 ];
        
        #get the existing GT call
        my @samples = keys( %{$$entry{gtypes}});
        my $gt = $$entry{gtypes}{$samples[ 0 ]}{GT};
        
        foreach my $sample ( sort( keys( %sampleBAM ) ) )
        {
            my $quality = Utilities::genotypeRegion($chr_, $start, $end, $sampleBAM{ $sample }, $minMapQ, $minReads );
            if( $call )
	        {
	            $$entry{gtypes}{$sample}{GT} = $gt;
                $$entry{gtypes}{$sample}{GQ} = $quality;
	        }
	        else
	        {
	            $$entry{gtypes}{$sample}{GT} = qq[./.];
                $$entry{gtypes}{$sample}{GQ} = '.';
	        }
        }
        $vcf->format_genotype_strings($$entry);
        print $out $vcf_out->format_line($$entry);
    }
    $vcf->close();
	close( $out );
}

#***************************INTERNAL HELPER FUNCTIONS********************

sub _getCandidateTEReadNames
{
    my $bam = shift;
    my $readgroups = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $minAnchor = shift;
    my $candidatesFasta = shift;
    my $candidatesBed = shift;
    
    my %candidates;
    open( my $ffh, qq[>>$candidatesFasta] ) or die qq[ERROR: Failed to create fasta file: $!\n];
    open( my $afh, qq[>>$candidatesBed] ) or die qq[ERROR: Failed to create anchors file: $!\n];
    
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
        
        #            read is not a duplicate        map quality is >= minimum
        if( ! ( $flag & $$BAMFLAGS{'duplicate'} ) && $qual >= $minAnchor )
        {
            #           read is mapped                       mate is unmapped
            if( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ( $flag & $$BAMFLAGS{'mate_unmapped'} ) )
            {
               $candidates{ $name } = $flag & $$BAMFLAGS{'1st_in_pair'} ? 2 : 1; #mate is recorded
               
               my $pos = $sam[ 3 ];
               print $afh qq[$ref\t$pos\t].($pos+$readLen).qq[\t$name\n];
            }
            #            read is mapped                         mate is mapped                                  not paired correctly
            elsif( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $flag & $$BAMFLAGS{'read_paired'} ) && $mref ne '=' )
            {
               $candidates{ $name } = $flag & $$BAMFLAGS{'1st_in_pair'} ? 2 : 1; #mate is recorded
               
               my $pos = $sam[ 3 ];
               print $afh qq[$ref\t$pos\t].($pos+$readLen).qq[\t$name\n];
            }
        }
        if( $currentChr ne $ref ){print qq[Reading chromosome: $ref\n];$currentChr = $ref;}
    }
    close( $bfh );
    
    return \%candidates;
}

sub _checkDiscoveryOutput
{
    my $file = shift;
    
    open( my $tfh, $file ) or die $!;
    my $line = <$tfh>;chomp( $line );
    if( $line !~ /^$HEADER/ ){die qq[Malformed header of input file: $file\n];}
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
    
    my $header = _getVcfHeader($vcf_out);
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
                my $ci1 = $s[ 1 ] - $pos;
                my $ci2 = $s[ 2 ] - $pos;
                my $dir = defined( $s[ 5 ] ) && length( $s[ 5 ] ) > 0 ? $s[ 5 ] : 'NA';
                
                print $bfh qq[$s[0]\t$pos\t].($pos+1).qq[\t$type==$sample\t$s[4]\t$dir\n];
                
                my %out;
                $out{CHROM}  = $s[ 0 ];
                $out{POS}    = $pos;
                $out{ID}     = '.';
                $out{ALT}    = ['<INS:ME>'];
                $out{REF}    = $refbase;
                $out{QUAL}   = $s[ 4 ];
                $out{FILTER} = ['NOT_VALIDATED'];
                $out{INFO} = { SVTYPE=>'INS', NOT_VALIDATED=>undef, MEINFO=>qq[$type,$s[1],$s[2],$dir] };
                $out{FORMAT} = ['GT', 'GQ'];
                
                $out{gtypes}{$s[3]}{GT} = qq[<INS:ME>/<INS:ME>];
                $out{gtypes}{$s[3]}{GQ} = qq[$s[4]];
                
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
                my $ci1 = $s[ 1 ] - $pos;
                my $ci2 = $s[ 2 ] - $pos;
                my $dir = defined( $s[ 5 ] ) && length( $s[ 5 ] ) > 0 ? $s[ 5 ] : 'NA';
                
                print $bfh qq[$s[0]\t$pos\t].($pos+1).qq[\t$type==$sample\t$s[4]\t$dir\n];
                
                my %out;
                $out{CHROM}  = $s[ 0 ];
                $out{POS}    = $pos;
                $out{ID}     = '.';
                $out{ALT}    = [$refbase,'<INS:ME>'];
                $out{REF}    = $refbase;
                $out{QUAL}   = $s[ 4 ];
                $out{FILTER} = ['NOT_VALIDATED'];
                $out{INFO} = { SVTYPE=>'INS', NOT_VALIDATED=>undef, MEINFO=>qq[$type,$s[1],$s[2],$dir] };
                $out{FORMAT} = ['GT','GQ'];
                
                $out{gtypes}{$s[3]}{GT} = qq[$refbase/<INS:ME>];
                $out{gtypes}{$s[3]}{GQ} = qq[$s[4]];
                
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

#creates all the tags necessary for the VCF files
sub _getVcfHeader
{
    my $vcf_out = shift;
    
    $vcf_out->add_header_line({key=>'source',value=>'RetroSeq v'.$VERSION});
    
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    $vcf_out->add_header_line( {key=>'INFO',ID=>'SVTYPE',Number=>'1',Type=>'String', Description=>'Type of structural variant'} );
    
    ##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
    $vcf_out->add_header_line( {key=>'INFO',ID=>'MEINFO',Number=>'4',Type=>'String', Description=>'Mobile element info of the form NAME,START,END,POLARITY'} );
    
    ##ALT=<ID=INS:ME,Description="Insertion of a mobile element">
    $vcf_out->add_header_line( {key=>'ALT', ID=>'INS:ME', Type=>'String', Description=>"Insertion of a mobile element"} );
    
    ##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
    $vcf_out->add_header_line({key=>'FORMAT',ID=>'GT',Number=>'1',Type=>'String',Description=>"Genotype"});
    
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
    $vcf_out->add_header_line( {key=>'FORMAT', ID=>'GQ', Number=>'1', Type=>'Float', Description=>'Genotype quality'} );
    
    $vcf_out->add_header_line( {key=>'INFO', ID=>'NOT_VALIDATED', Number=>'0', Type=>'Flag', Description=>'Not validated experimentally'} );
    
    return $vcf_out->format_header();
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
    
    if( length( $base ) != 1 && $base !~ /acgtnACGTN/ ){die qq[Failed to get reference base at $chr:$pos-$pos: $base\n]}
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
    }
    print $ofh $FOOTER.qq[\n];
    close( $ofh );
}

