package Utilities;

use strict;
use warnings;
use Carp;

my $FILTER_WINDOW = 50;
my $BREAKPOINT_WINDOW = 250;

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

#status that can be returned from calling
our $UNKNOWN_FAIL = 0;
our $HOM_DEPTH_TOO_HIGH = 1;
our $NOT_ENOUGH_READS_CLUSTER = 2;
our $NOT_ENOUGH_READS_FLANKS = 3;
our $NOT_ENOUGH_READS_BLUE = 4;
our $NEITHER_SIDE_RATIO_PASSES = 5;
our $ONE_SIDED_RATIO_PASSES = 6;
our $DISTANCE_THRESHOLD = 7;
our $PASS = 8;
our $FILTER_NAMES = ['Fail', 'HomBreakDepth', 'MinReads', 'MinReadsFlanks', 'MinReadsBlue', 'BothRatioFail', 'OneRatioFail', 'DistanceThreshold' ];
our $FILTER_DESC = ['Unknown fail', 'HomBreakDepth', 'Not enough supporting reads', 'Not enough supporting reads on either flanking sides', 'Not enough supporting multi-mapped reads', 'Neither side has required ratio of fwd:rev reads', 'One side has required ratio of fwd:rev reads', 'Distance between 3\' and 5\' reads is greater than threshold' ];

our $VERSION = 0.2;

sub filterOutRegions
{
    my $inputBED = shift;
    my $filterBED = shift;
    my $outputBED = shift;

    my @calls;
    open( my $ifh, $inputBED ) or die $!;
    while( my $line = <$ifh> ){chomp($line);push(@calls,$line);}
    close( $ifh );

    open( my $ffh, $filterBED ) or die $!;
    my $filtered = 0;
    while( my $region = <$ffh> )
    {
        chomp( $region );
        my @s = split( /\t/, $region );
        foreach(my $i=0;$i<@calls;$i++)
        {
	    next unless $calls[$i];
            my @s1 = split( /\t/, $calls[$i] );
            croak qq[Badly formatted call: $calls[$i]] unless @s1==5;
            if( $s[ 0 ] eq $s1[ 0 ] && ( 
                ( $s1[ 1 ] > $s[ 1 ] && $s1[ 2 ] < $s[ 2 ] ) #TE call is totally enclosed in region
                ||
                ( abs( $s1[ 1 ] - $s[ 1 ] ) < $FILTER_WINDOW )
                ||
                ( abs( $s1[ 1 ] - $s[ 2 ] ) < $FILTER_WINDOW )
                ||
                ( abs( $s1[ 2 ] - $s[ 1 ] ) < $FILTER_WINDOW )
                ||
                ( abs( $s1[ 2 ] - $s[ 2 ] ) < $FILTER_WINDOW )
                ||
		( $s1[ 1 ] > $s[ 1 ] && $s1[ 1 ] < $s[ 2 ] )
		||
		( $s1[ 2 ] > $s[ 1 ] && $s1[ 2 ] < $s[ 2 ] )
		||
                ( $s[ 1 ] > $s1[ 1 ] && $s[ 2 ] < $s1[ 2 ] )  #TE in ref in enclosed within the called TE region
                )
              )
            {
                print qq[Excluding region: $calls[ $i ]\n];
                undef( $calls[ $i ] );
                $filtered ++;
            }
        }
    }
    close( $ffh );

    open( my $ofh, qq[>$outputBED] ) or die $!;
    foreach my $call ( @calls ){if(defined( $call ) ){ print $ofh qq[$call\n]; } }
    close( $ofh );
    
    return $filtered;
}

sub getCandidateBreakPointsDepth
{
    my $chr = shift;
    my $start = shift;
    my $end = shift;
	my @bams = @{ $_[ 0 ] };shift;
    my $ref = shift;
    my $dfh = shift; #handle to output stream for logging calls that fall out

	my $bamStr = '';
	foreach my $bam( @bams ){$bamStr.=qq[ $bam];}
    
    open( my $tfh, qq[samtools mpileup -r $chr:$start-$end -f $ref -u $bamStr | bcftools view - | ] ) or die $!;
	my %depths;
	while( <$tfh> )
	{
	    chomp;next if($_=~/^#/);my @s = split( /\t/, $_ );
	        
	    if( $s[ 4 ] =~ /^X$/ ) #not a snp call
	    {
	        #get the depth at the position
	        my @tags = split( /;/, $s[ 7 ] );
	        foreach( @tags ){if($_=~/^(DP=)(\d+)$/){$depths{ $s [ 1 ] } = $2;}}
	    }
	}
	close( $tfh );
	
	my @res = _local_min_max( %depths );
	if( !@res || !$res[ 0 ] ){print qq[WARNING: no max/min returned for $chr:$start-$end\n];print $dfh qq[$chr\t$start\t$end\tnodepth\n\n];return undef;}
	my %min = %{$res[ 0 ]};
	my @positions = keys( %min );
	
	#sort by min depths
	@positions = sort {$min{$a}<=>$min{$b}} @positions;
	
	return [\@positions,\%min];
}

sub _local_min_max
{
    my %depths = @_;
    return undef unless keys( %depths ) > 1;
    
    my %minima = ();
    my %maxima = ();
    my $prev_cmp = 0;
    
    my @positions = keys( %depths );
    for( my $i=0;$i<@positions-1;$i++)
    {
        my $cmp = $depths{$positions[$i]} <=> $depths{$positions[$i+1]};
        if ($cmp && $cmp != $prev_cmp) 
        {
            $minima{ $positions[ $i ] } = $depths{ $positions[ $i ] };
            $maxima{ $positions[ $i ] } = $depths{ $positions[ $i ] };
            $prev_cmp = $cmp;
        }
    }
    
    $minima{ $positions[ -1 ] } = $depths{ $positions[ -1 ] } if $prev_cmp >= 0;
    $maxima{ $positions[ -1 ] } = $depths{ $positions[ -1 ] } if $prev_cmp >= 0;
    
    return (\%minima, \%maxima);
}

sub testBreakPoint
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 9;
    
    my $chr = shift;
    my $refPos = shift;
    die qq[Non integer value passed as refPos] unless defined( $refPos ) && $refPos =~ /^\d+$/;
	my @bams = @{ $_[ 0 ] };shift;
    my $minMapQ = shift;
    my $originalCall = shift;
    my $dfh = shift; #file handle to print out info on failed calls
    my $ignoreRGs = shift; #file of lines with "RG:tag\t"
    my $minReads = shift;
    my $genotypeMode = shift; #0/1 saying whether to operate in genotyping mode (i.e. less stringent criteria)
    
    my @originalCallA = split( /\t/, $originalCall );
    
    #test to see if lots of rp's either side
    my $lhsFwdBlue = 0; my $lhsRevBlue = 0; my $rhsFwdBlue = 0; my $rhsRevBlue = 0;
    my $lhsFwdGreen = 0; my $lhsRevGreen = 0; my $rhsFwdGreen = 0; my $rhsRevGreen = 0;
	
    #store the last blue read before the b/point, and first blue read after the b/point
	my $lastBluePos = 0;my $firstBluePos = 100000000000;
	
	my $cmdpre;
    if( @bams > 1 )
    {
        if( _mergeRegionBAMs( \@bams, $chr, $refPos-$BREAKPOINT_WINDOW - 500, $refPos+$BREAKPOINT_WINDOW + 500, qq[/tmp/$$.region.bam] ) )
        {
            $cmdpre = qq[samtools view /tmp/$$.region.bam $chr:];
        }else{die qq[Failed to extract region $chr:$refPos BAM];}
    }
    else
    {
        $cmdpre = qq[samtools view $bams[0] $chr:];
    }
    
    #compute the distance from the breakpoint to the last fwd read and the first rev read (and then add this to the read window buffer)
    
    
	#also check the orientation of the supporting reads (i.e. its not just a random mixture of f/r reads overlapping)
	my $cmd = $cmdpre.($refPos-$BREAKPOINT_WINDOW).qq[-].($refPos+$BREAKPOINT_WINDOW).qq[ | ].(defined($ignoreRGs) ? qq[ grep -v -f $ignoreRGs |] : qq[]);
	open( my $tfh, $cmd ) or die $!;
	while( my $sam = <$tfh> )
	{
	    chomp( $sam );
	    my @s = split( /\t/, $sam );
        next unless $s[ 4 ] > $minMapQ;
        #        the mate is mapped                       or mate ref name is different chr
        if( !($s[ 1 ] & $$BAMFLAGS{'mate_unmapped'}) && ( $s[ 6 ] ne '=' ) )
        {
            if( ( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) )  #rev strand
            {
                if( $s[ 3 ] < $refPos ){$lhsRevBlue++;}else{$rhsRevBlue++;$firstBluePos = $s[ 3 ] if( $s[ 3 ] < $firstBluePos );}
            }
            else
            {
                if( $s[ 3 ] < $refPos ){$lhsFwdBlue++;$lastBluePos = $s[ 3 ] + length( $s[ 9 ] ) if( ( $s[ 3 ] + length( $s[ 9 ] ) ) > $lastBluePos );}else{$rhsFwdBlue++;}
            }
        }
        #        the mate is unmapped
        elsif( $s[ 1 ] & $$BAMFLAGS{'mate_unmapped'} )
        {
            if( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) #rev strand
            {
                if( $s[ 3 ] < $refPos ){$lhsRevGreen++;}else{$rhsRevGreen++;}
            }
            else
            {
                if( $s[ 3 ] < $refPos ){$lhsFwdGreen++;}else{$rhsFwdGreen++;}
            }
        }
    }
    
    unlink( qq[/tmp/$$.region.bam] );
    
    #check there are supporting read pairs either side of the depth minima
    my $lhsRev = $lhsRevGreen + $lhsRevBlue;my $rhsRev = $rhsRevGreen + $rhsRevBlue;my $lhsFwd = $lhsFwdGreen + $lhsFwdBlue;my $rhsFwd = $rhsFwdGreen + $rhsFwdBlue;
    my $dist = $firstBluePos - $lastBluePos;
    
    my $minBlue = int($minReads / 2);

    print $dfh qq[TEST: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
    
    if( ! $genotypeMode ) #calling mode
    {
        my $callString = qq[$chr\t$refPos\t].($refPos+1).qq[\t$originalCallA[ 3 ]\t$originalCallA[ 4 ]];
        
        #want to be more descriptive with filter values to return e.g. total reads too low, one-side ok only, sides OK but distance too large
        my $lhsRatioPass = ( $lhsRevBlue == 0 ) || ( $lhsRevBlue > 0 && $lhsFwdBlue / $lhsRevBlue > 2 );
        my $rhsRatioPass = ( $rhsFwdBlue == 0 ) || ( $rhsFwdBlue > 0 && $rhsRevBlue / $rhsFwdBlue > 2 );
        
        if( ( $lhsFwd + $rhsRev ) < $minReads )
        {
            return [$NOT_ENOUGH_READS_CLUSTER, $callString, 0];
        }
        elsif( $lhsFwd < $minReads || $rhsRev < $minReads )
        {
            return [$NOT_ENOUGH_READS_FLANKS, $callString, 0];
        }
        elsif( $lhsFwdBlue < $minBlue || $rhsRevBlue < $minBlue )
        {
            return [$NOT_ENOUGH_READS_BLUE, $callString, 0];
        }
        elsif( ! $lhsRatioPass && ! $rhsRatioPass )
        {
            return [$NEITHER_SIDE_RATIO_PASSES, $callString, 0];
        }
        elsif( ! ( $lhsRatioPass && $rhsRatioPass ) )
        {
            return [$ONE_SIDED_RATIO_PASSES, $callString, 0];
        }
        elsif( $dist > 120 )
        {
            return [$DISTANCE_THRESHOLD, $callString, 0];
        }
        else
        {
            my $ratio = ( $lhsRev + $rhsFwd ) / ( $lhsFwd + $rhsRev );
            return [$PASS, $callString, $ratio];
        }
        
        return [$UNKNOWN_FAIL, $callString, 10000 ];
=pod	            
                if( $lhsFwdBlue >= $minBlue && $rhsRevBlue >= $minBlue && $lhsFwd >= $minReads && $rhsRev >= $minReads && ( $lhsRevBlue == 0 || $lhsFwdBlue / $lhsRevBlue > 2 ) && ( $rhsFwdBlue == 0 || $rhsRevBlue / $rhsFwdBlue > 2 ) && $dist < 120 )
                {
                    my $ratio = ( $lhsRev + $rhsFwd ) / ( $lhsFwd + $rhsRev ); #objective function is to minimise this value (i.e. min depth, meets the criteria, and balances the 3' vs. 5' ratio best)
                    my $callString = qq[$chr\t$refPos\t].($refPos+1).qq[\t$originalCallA[ 3 ]\t$originalCallA[ 4 ]\n];
                    print $dfh qq[PASS: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
                    
                    return [$ratio,$callString];
                }
                else
                {
                    print $dfh qq[FILTER: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
                    return undef;
                }
=cut
            }
            else #genotyping mode - less stringent num of reads required
            {
                if( ( $lhsFwdBlue >= $minBlue && $lhsFwd >= $minReads && ( $lhsRevBlue == 0 || $lhsFwdBlue / $lhsRevBlue > 2 ) )#&& $dist < 120 )
                    ||
                    ( $rhsRevBlue >= $minBlue && $rhsRev >= $minReads && ( $rhsFwdBlue == 0 || $rhsRevBlue / $rhsFwdBlue > 2 ) )#&& $dist < 120 )
                  )
                {
                    my $ratio = ( $lhsRev + $rhsFwd ) / ( $lhsFwd + $rhsRev ); #objective function is to minimise this value (i.e. min depth, meets the criteria, and balances the 3' vs. 5' ratio best)
                    my $callString = qq[$chr\t$refPos\t].($refPos+1).qq[\t$originalCallA[ 3 ]\t$originalCallA[ 4 ]\n];
                    print $dfh qq[PASS: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
                    
                    return [$lhsFwd+$rhsRev]; #return number of reads supporting call
                }
                else
                {
                    print $dfh qq[FILTER: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
                    return undef;
                }
            }
}

sub genotypeRegion
{
    croak qq[Incorrect number of arguments: ].scalar(@_) unless @_ == 5;
    
    my $chr = shift;
    my $position = shift;
    my $bam = shift;
    my $minMapQ = shift;
    my $minReads = shift;
    
    die qq[cant find bam: $bam\n] unless -f $bam;
    
    #also check the orientation of the supporting reads (i.e. its not just a random mixture of f/r reads overlapping)
    my $cmd = qq[samtools view $bam $chr:].($position-$BREAKPOINT_WINDOW).qq[-].($position+$BREAKPOINT_WINDOW).qq[ | ];
	open( my $tfh, $cmd ) or die $!;
	my $lhsRev = 0; my $rhsRev = 0; my $lhsFwd = 0; my $rhsFwd = 0;
	while( <$tfh> )
	{
	    chomp;
	    my @s = split( /\t/, $_ );
	    next unless $s[ 4 ] > $minMapQ;
	    if( !( $s[ 1 ] & $$BAMFLAGS{'read_paired'} ) && ( $s[ 1 ] & $$BAMFLAGS{'mate_unmapped'} || $s[ 2 ] ne $s[ 6 ] ) )
	    {
	        if( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) #rev strand
	        {
	            if( $s[ 3 ] < $position ){$lhsRev++;}else{$rhsRev++;}
	        }
	        else
	        {
	            if( $s[ 3 ] < $position ){$lhsFwd++;}else{$rhsFwd++;}
	        }
	    }
	 }
	 
	 #N.B. Key difference is that only 1 side is required to have the correct ratio of fw:rev reads
	 if( $lhsFwd > $minReads && $rhsRev > $minReads && ( ( $lhsRev == 0 || $lhsFwd / $lhsRev > 2 ) || ( $rhsFwd == 0 || $rhsRev / $rhsFwd > 2 ) ) )
	 {
	     return ($lhsFwd+$rhsRev);
	 }
	 return undef;
}

sub getCandidateBreakPointsDirVote
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 5;
    
    my $chr = shift;
    my $start = shift;
    my $end = shift;
	my @bams = @{ $_[ 0 ] };shift;
    my $minQ = shift;
    
    if( $start !~ /^\d+$/ || $end !~ /^\d+$/ ){die qq[ERROR: Invalid parameters passed to getCandidateBreakPointsDir: $chr $start $end\n];}
    
    my %fwdCount;
    my %revCount;
    
    my $cmd;
    if( @bams > 1 )
    {
        if( _mergeRegionBAMs( \@bams, $chr, $start, $end, qq[/tmp/$$.region.bam] ) )
        {
            $cmd = qq[samtools view /tmp/$$.region.bam |];
        }else{die qq[Failed to extract region $chr:$start-$end BAM];}
    }
    else
    {
        $cmd = qq[samtools view $bams[0] $chr:$start-$end |];
    }
    
    open( my $tfh, $cmd ) or die $!;
    my $fwdCurrentPos = $start;
    my $revCurrentPos = $start;
    $fwdCount{ $fwdCurrentPos - 1 } = 0;
    $revCount{ $revCurrentPos - 1 } = 0;
    my %endCounts;
    while( my $sam = <$tfh> )
    {
        chomp( $sam );
        my @samL = split( /\t/, $sam );
        my $flag = $samL[ 1 ];
        
        next if $samL[ 4 ] < $minQ;
        next if ( $flag & $$BAMFLAGS{'duplicate'} );

        #update the counts to the current position
        my $gap = 0; #apply a gap-open penalty to the score to avoid odd spurious read throwing off counts
        for(my $i=$fwdCurrentPos;$i<$samL[3];$i++ )
        {
            $gap++;
            my $add = 0;if( $endCounts{ $i } ){$add=$endCounts{ $i };}
            if($gap%100==0&&$fwdCount{$i-1}>0)
            {
                $fwdCount{$i}=$fwdCount{$i-1} - 1 + $add;
            }
            else
            {
                $fwdCount{$i}=$fwdCount{$i-1} + $add;
            }
        }
        for(my $i=$revCurrentPos;$i<$samL[3];$i++ ){$gap++;if($gap%100==0&&$revCount{$i-1}<0){$revCount{$i}=$revCount{$i-1} + 1;}else{$revCount{$i}=$revCount{$i-1};}}
        
        #       paired technology                       not paired correctly                    is on fwd strand
        if( ( $flag & $$BAMFLAGS{'paired_tech'} ) && !( $flag & $$BAMFLAGS{'read_paired'} ) && !( $flag & $$BAMFLAGS{ 'reverse_strand' } ) ) #fwd supporting read
        {
            no warnings 'uninitialized'; #perl warns when the value is 0
            my $endPos = $samL[ 3 ] + length( $samL[ 9 ] );
            if( $endCounts{ $endPos } ){$endCounts{$endPos}++;}else{$endCounts{$endPos} = 1;};
        }
        #       paired technology                       not paired correctly                    is on rev strand
        elsif( ( $flag & $$BAMFLAGS{'paired_tech'} ) && !( $flag & $$BAMFLAGS{'read_paired'} ) && ( $flag & $$BAMFLAGS{ 'reverse_strand' } ) ) #rev supporting read
        {
            no warnings 'uninitialized'; #perl warns when the value is 0
            $revCount{ $samL[ 3 ] } = $revCount{ $samL[ 3 ] - 1 } - 1;
            $revCurrentPos = $samL[ 3 ] + 1;
        }
    }
    close( $tfh );
    
    #offset the rev array
    my @s = sort {$a<=>$b} ( values( %revCount ) );
    my $min = $s[ 0 ];
    my @keys = keys(%revCount);
    foreach my $key(@keys){$revCount{$key}+=abs($min);}
    
    #now find the position where the sum of the two counts maximises
    my $maxPos = $start;my @maxPoss;my $maxVal = $fwdCount{$start}+$revCount{$start};
    
    for(my $i=$start;$i<$revCurrentPos;$i++)
    {
        my $sum = $fwdCount{$i}+$revCount{$i};
        print qq[$i\t$sum\n];
        if( $fwdCount{$i}+$revCount{$i} >= $maxVal )
        {
            $maxPos = $i;
            @maxPoss = ();
            push( @maxPoss, $i );
            $maxVal = $fwdCount{$i}+$revCount{$i};
        }
        elsif( $fwdCount{$i}+$revCount{$i} == $maxVal )
        {
            push( @maxPoss, $i );
        }
    }
    
    return undef if( ! @maxPoss );
    
    $maxPos = $maxPoss[ int(scalar(@maxPoss)/2) ];
    
    #get the depth at the position
    my $depth = 100;
    $cmd = '';
    if( @bams > 1 )
    {
        $cmd = qq[samtools mpileup -r $chr:$maxPos-$maxPos /tmp/$$.region.bam | ];
    }
    else
    {
        $cmd = qq/samtools mpileup -r $chr:$maxPos-$maxPos $bams[0] | /;
    }
    
    open( my $ifh, $cmd ) or die qq[failed on samtools cmd: $cmd\n];
    my $mpileupOut = <$ifh>;
    @s = split( /\t/, $mpileupOut );
    
    return [$maxPos, length($s[4]) ];
}

#convert the individual read calls to calls for putative TE insertion calls
#output is a BED file and a VCF file of the calls
sub convertToRegionBED
{
	my $calls = shift;
	my $minReads = shift;
	my $id = shift;
	my $max_gap = shift;
	my $outputbed = shift;
	
	open( my $ofh, qq[>$outputbed] ) or die $!;
	open( my $cfh, $calls ) or die $!;
	my $lastEntry = undef;
	my $regionStart = 0;
	my $regionEnd = 0;
	my $regionChr = 0;
	my $reads_in_region = 0;
	my %startPos; #all of the reads must start from a different position (i.e. strict dup removal)
	while( <$cfh> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		if( ! defined $lastEntry )
		{
			$regionStart = $s[ 1 ];
			$regionEnd = $s[ 2 ];
			$regionChr = $s[ 0 ];
			$reads_in_region = 1;
			$startPos{ $s[ 1 ] } = 1;
		}
		elsif( $s[ 0 ] ne $regionChr )
		{
			#done - call the region
			my @s1 = split( /\t/, $lastEntry );
			$regionEnd = $s1[ 2 ];
			my $size = $regionEnd - $regionStart; $size = 1 unless $size > 0;
			print $ofh "$regionChr\t$regionStart\t$regionEnd\t$id\t$reads_in_region\n" if( $reads_in_region >= $minReads );
			
			$reads_in_region = 1;
			$regionStart = $s[ 1 ];
			$regionEnd = $s[ 2 ];
			$regionChr = $s[ 0 ];
			%startPos = ();
			$startPos{ $s[ 1 ] } = 1;
		}
		elsif( $s[ 1 ] - $regionEnd > $max_gap )
		{
			#call the region
			my @s1 = split( /\t/, $lastEntry );
			$regionEnd = $s1[ 2 ];
			my $size = $regionEnd - $regionStart; $size = 1 unless $size > 0;
			print $ofh "$regionChr\t$regionStart\t$regionEnd\t$id\t$reads_in_region\n" if( $reads_in_region >= $minReads );
			
			$reads_in_region = 1;
			$regionStart = $s[ 1 ];
			$regionChr = $s[ 0 ];
			
			%startPos = ();
			$startPos{ $s[ 1 ] } = 1;
		}
		else
		{
			#read is within the region - increment
			if( ! defined( $startPos{ $s[ 1 ] } ) )
			{
				$reads_in_region ++;
				$regionEnd = $s[ 2 ];
				$startPos{ $s[ 1 ] } = 1;
			}
		}
		$lastEntry = $_;
	}
	my $size = $regionEnd - $regionStart; $size = 1 unless $size > 0;
	
	print $ofh "$regionChr\t$regionStart\t$regionEnd\t$id\t$reads_in_region\n" if( $reads_in_region >= $minReads );
	close( $cfh );
	close( $ofh );
	
	return 1;
}

sub annotateCallsBED
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 2;
    
    my $sortedCallsBED = shift;
    my $sortedReadsBED = shift;
    print qq[$sortedCallsBED\t$sortedReadsBED\n];
    #read through the supporting reads bed and pick up the reads in the vicinity of the calls
    open( my $sfh, $sortedReadsBED ) or die qq[ERROR: Failed to open sorted reads BED: $sortedReadsBED];
    open( my $cfh, $sortedCallsBED ) or die qq[ERROR: Failed to open calls bed: $sortedCallsBED];
    open( my $ofh, qq[>$sortedCallsBED.dir] ) or die qq[ERROR: Failed to create new BED file];
    my $currentCall = <$cfh>;
    if( ! $currentCall ){print qq[Found zero calls to orientate in $sortedCallsBED - returning\n];close($ofh);return;}
    
    chomp($currentCall);my ($techr,$testart,$teend,$tename,$tescore) = split( /\t/, $currentCall );
    my $currentCallLHS = 0;my $currentCallRHS = 0;my $withinCall = 0;
    my $fwd = 1;
    while( <$sfh> )
    {
        chomp;
        my ($chr,$start,$end,$name,$orientation) = split( /\t/, $_ );
        
            if( $chr eq $techr )
            {
                if( $end < $testart && $start > $testart - $BREAKPOINT_WINDOW )
                {
                    if( $orientation eq '+' ){$currentCallLHS++;}else{$currentCallLHS--;}
                    $withinCall = 1;
                }
                elsif( $chr eq $techr && $start > $teend && $end < $teend + $BREAKPOINT_WINDOW )
                {
                    if( $orientation eq '+' ){$currentCallRHS++;}else{$currentCallRHS--;}
                    $withinCall = 1;
                }
            }
            
            if( $withinCall && ($end > $teend + $BREAKPOINT_WINDOW || $techr ne $chr ) )
            {
                $fwd = 1;
                if( $currentCallLHS > 0 && $currentCallRHS < 0 ){$fwd = 0;} #if it looks like a reverse insertions
                print $ofh qq[$techr\t$testart\t$teend\t$tename\t$tescore\t].($fwd > -1 ? qq[+\n] : qq[-\n] );
                $withinCall = 0;
                $currentCallLHS = 0;
                $currentCallRHS = 0;
                $currentCall = <$cfh>;
                last if( ! $currentCall );
                chomp($currentCall);($techr,$testart,$teend,$tename,$tescore) = split( /\t/, $currentCall );
            }
    }
    if( $withinCall ){if( $currentCallLHS > 0 && $currentCallRHS < 0 ){$fwd = 0;}print $ofh qq[$techr\t$testart\t$teend\t$tename\t$tescore\t].($fwd > -1 ? qq[+\n] : qq[-\n]);}
    
    close( $ofh );
    close( $sfh );
    
    #rename the output file to the input file
    rename(qq[$sortedCallsBED.dir],qq[$sortedCallsBED]) == 1 or die qq[Failed to rename orientated BED file: $sortedCallsBED.dir to $sortedCallsBED\n];
}

sub getBAMSampleName
{
    my @bams = @{ $_[ 0 ] };
    
    my %samples;
    foreach my $bam ( @bams )
    {
        open( my $bfh, qq[samtools view -H $bam |] );
        while( my $line = <$bfh> )
        {
            chomp( $line );
            next unless $line =~ /^\@RG/;
            my @s = split(/\t/, $line );
            foreach my $tag (@s)
            {
                if( $tag && $tag =~ /^(SM):(\S+)/ && ! $samples{ $2 } ){$samples{ $2 } = 1;}
            }
        }
    }
    
    my $sampleName = 'unknown';
    if( %samples && keys( %samples ) > 0 )
    {
        $sampleName = join( "_", keys( %samples ) ); #if multiple samples - join the names into a string
        print qq[Found sample: $sampleName\n];
    }
    else
    {
        print qq[WARNING: Cant determine sample name from BAM - setting to unknown\n];
        $sampleName = 'unknown';
    }
    
    return $sampleName; 
}

sub _mergeRegionBAMs
{
    my @bams = @{ $_[ 0 ] };shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $output = shift;
    
    if( @bams > 1 )
    {
        my $bamStr = '';
        foreach my $bam( @bams ){$bamStr.=qq[ $bam];}
        
        #buffer the regions from the bams and then merge into a single bam
        my $str = '';
        for(my $i=0;$i<@bams;$i++)
        {
            my $ind_bam = qq[/tmp/$$.$i.ind_region.bam];
            system( qq[samtools view -b -h $bams[$i] $chr:$start-$end > $ind_bam] ) == 0 or die qq[Failed to get region $chr:$start-$end for $bams[$i] $ind_bam\n];
            $str .= qq[ $ind_bam];
        }
        system( qq[samtools merge -f $output $str;samtools index $output] ) == 0 or die qq[Failed to merge region BAMs for region $chr:$start-$end : $str\n];
        for(my $i=0;$i<@bams;$i++)
        {
            unlink( qq[/tmp/$$.$i.ind_region.bam] );
        }
    }
    else
    {
        system( qq[samtools view $bams[0] $chr:$start-$end > $output] ) == 0 or die qq[Failed to extract region $chr:$start-$end BAM: $bams[0]\n];
    }
    
    return 1;
}

sub checkBinary
{
    my $binary = shift;
    my $version = undef;
    if( @_ == 1 )
    {
        $version = shift;
    }
    
    if( ! `which $binary` )
    {
        croak qq[Error: Cant find required binary $binary\n];
    }
    
    if( $version )
    {
        my @v = split( /\./, $version );
        my $hasOutput = 0;
        open( my $bfh, qq[ $binary 2>&1 | ] ) or die "failed to run $binary\n";
        while(<$bfh>)
        {
            chomp;
            if( lc($_) =~ /version/ && $_ =~ /(\d+)\.(\d+)\.(\d+)/ )
            {
                if( $1.'.'.$2 < $v[ 0 ].'.'.$v[ 1 ] )
                {
                    die qq[\nERROR: $binary version $version is required - your version is $1.$2.$3\n];
                }
                $hasOutput = 1;
            }
        }
        close( $bfh );
        die qq[ERROR: Cant determine version number of $binary\n] unless $hasOutput;
    }
}

#creates all the tags necessary for the VCF files
sub getVcfHeader
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

    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
    $vcf_out->add_header_line( {key=>'FORMAT', ID=>'FL', Number=>'1', Type=>'Integer', Description=>'Call Status - for reference calls a flag to say if the call failed a particular filter. Filters are ordered by priority in calling (higher number indicates closer to being called). 1 - depth too high in region, 2 - not enough reads in cluster, 3 - not enough total flanking reads, 4 - not enough inconsistently mapped reads, 5 - neither side passes ratio test, 6 - one side passes ratio test, 7 - distance too large at breakpoint, 8 - PASSED all filters'} );
    
    $vcf_out->add_header_line( {key=>'INFO', ID=>'NOT_VALIDATED', Number=>'0', Type=>'Flag', Description=>'Not validated experimentally'} );
    $vcf_out->add_header_line( {key=>'INFO', ID=>'1000G', Number=>'0', Type=>'Flag', Description=>'Overlaps with 1000G MEI call'} );
    $vcf_out->add_header_line( {key=>'INFO', ID=>'REPEATMASKER', Number=>'0', Type=>'Flag', Description=>'Overlaps with a reference ME element'} );

    return $vcf_out->format_header();
}

1;
