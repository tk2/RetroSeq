package RetroSeq::Utilities;
=pod
This file is part of RetroSeq.

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RetroSeq.  If not, see <http://www.gnu.org/licenses/>.
=cut

use strict;
use warnings;
use Carp;
use File::Basename;

my $FILTER_WINDOW = 200;
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
our $UNKNOWN_FAIL = -1;
our $INV_BREAKPOINT = 0;
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

our $VERSION = 1.34;

sub filterOutRegions
{
    my $inputBED = shift;
    my $filterBED = shift;
    my $outputBED = shift;
    
    #generate a new BED from the input BED using the mid points (approx. breakpoint) and use this to filter off
    open( my $tfh, qq[>$$.tofilter.bed] ) or die $!;
    open( my $ifh, qq[$inputBED] ) or die $!;
    my $start = 0;
    while( my $l = <$ifh> )
    {
        chomp( $l );
        my @s = split( /\s+/, $l );
        print $tfh qq[$s[0]\t].int((($s[1]+$s[2])/2)).qq[\t].int((($s[1]+$s[2])/2)+1).qq[\t$s[1]\n];
        $start ++;
    }
    close( $ifh );close( $tfh );
    
    #run windowBED to do this - much quicker than the code below!
    my %keep;
    open( my $kfh, qq[bedtools window -a $$.tofilter.bed -b $filterBED -w 30 -v | ] ) or die $!;
    while( my $l = <$kfh> )
    {
        chomp( $l );
        my @s = split( /\s+/, $l );
        $keep{ $s[ 3 ] } = 1;
    }
    close( $kfh );
    
    #now get these calls only and return them
    open( my $ofh, qq[>$outputBED] ) or die $!;
    open( $ifh, qq[$inputBED] ) or die $!;
    while( my $l = <$ifh> )
    {
        chomp( $l );
        my @s = split( /\s+/, $l );
        print $ofh qq[$l\n] if( $keep{ $s[ 1 ] } );
    }
    close( $ifh );
    
    my $retained = scalar( keys( %keep ) );
    print qq[Started with $start Retained $retained\n];
    return $retained;
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
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 10;
    
    my $chr = shift;
    my $refPos = shift;
    die qq[Non integer value passed as refPos] unless defined( $refPos ) && $refPos =~ /^\d+$/;
	my @bams = @{ $_[ 0 ] };shift;
    my $minMapQ = shift;
    my $originalCall = shift;
    my $dfh = shift; #file handle to print out info on failed calls
    my $ignoreRGs = shift; #file of lines with "RG:tag\t"
    my $minReads = shift;
    my $minSoftClip = shift;
    my $genotypeMode = shift; #0/1 saying whether to operate in genotyping mode (i.e. less stringent criteria)
    
    my @originalCallA = split( /\t/, $originalCall );
    my $lhsWindow = $refPos - $originalCallA[ 1 ]; $lhsWindow = 400 if( $lhsWindow < 400 );
    my $rhsWindow = $originalCallA[ 2 ] - $refPos; $rhsWindow = 400 if( $rhsWindow < 400 );
    
    #test to see if lots of rp's either side
    my $lhsFwdBlue = 0; my $lhsRevBlue = 0; my $rhsFwdBlue = 0; my $rhsRevBlue = 0;
    my $lhsFwdGreen = 0; my $lhsRevGreen = 0; my $rhsFwdGreen = 0; my $rhsRevGreen = 0;
	
    #store the last blue read before the b/point, and first blue read after the b/point
	my $lastBluePos = 0;my $firstBluePos = 100000000000;
	
	my $cmdpre;
    if( @bams > 1 )
    {
        if( _mergeRegionBAMs( \@bams, $chr, $refPos-$lhsWindow-400, $refPos+$rhsWindow+400, qq[/tmp/$$.region.bam] ) )
        {
            $cmdpre = qq[samtools view /tmp/$$.region.bam $chr:];
        }else{die qq[Failed to extract region $chr:$refPos BAM];}
    }
    else
    {
        $cmdpre = qq[samtools view $bams[0] $chr:];
    }
    
    #compute the distance from the breakpoint to the last fwd read and the first rev read (and then add this to the read window buffer)
    
	#also check the orientation of the supporting reads (i.e. its not just a random mixture of f/r reads overlapping) AND count the number of FF and RR read pairs to remove biases near inversion breakpoints
	my $numSameOrientation = 0;
	my $cmd = $cmdpre.($refPos-$lhsWindow).qq[-].($refPos+$rhsWindow).qq[ | ].(defined($ignoreRGs) ? qq[ grep -v -f $ignoreRGs |] : qq[]);
	open( my $tfh, $cmd ) or die $!;
	my $totalLFwd = 0; my $totalRRev = 0;my $totalSup = 0;
	while( my $sam = <$tfh> )
	{
	    chomp( $sam );
	    my @s = split( /\t/, $sam );
	    
	    my $supporting = isSupportingClusterRead( $s[ 1 ], $s[ 8 ], $s[ 4 ], $minMapQ, $minSoftClip, $s[ 5 ] );
	    if( $supporting  == 2 )
	    {$totalSup++;
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
        elsif( $supporting == 1 )
        {$totalSup++;
            if( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) #rev strand
            {
                if( $s[ 3 ] < $refPos ){$lhsRevGreen++;}else{$rhsRevGreen++;}
            }
            else
            {
                if( $s[ 3 ] < $refPos ){$lhsFwdGreen++;}else{$rhsFwdGreen++;}
            }
        }
        
        #       read is mapped                                  mate is mapped                           not paired correctly                  ins size < 50k               both mates are mapped to same strand
        if( ! ($s[ 1 ] & $$BAMFLAGS{'unmapped'} ) && !($s[ 1 ] & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $s[ 1 ] & $$BAMFLAGS{'read_paired'} ) && $s[ 8 ] < 50000 && ( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) == ( $s[ 1 ] & $$BAMFLAGS{'mate_reverse'} ) )
        {
            $numSameOrientation ++;
        }
    }
    unlink( qq[/tmp/$$.region.bam] );
    
    #if it appears to be an inversion breakpoint
    my $callString = qq[$chr\t$refPos\t].($refPos+1).qq[\t$originalCallA[ 3 ]\t$originalCallA[ 4 ]];
    print qq[No same orientation: $numSameOrientation\n];
    if( $numSameOrientation > $originalCallA[ 4 ] )
    {
        return [$INV_BREAKPOINT, $callString, 10000 ];
    }
    
    #check there are supporting read pairs either side of the depth minima
    my $lhsRev = $lhsRevGreen + $lhsRevBlue;my $rhsRev = $rhsRevGreen + $rhsRevBlue;my $lhsFwd = $lhsFwdGreen + $lhsFwdBlue;my $rhsFwd = $rhsFwdGreen + $rhsFwdBlue;
    my $dist = $firstBluePos - $lastBluePos;
    my $minBlue = int($minReads / 2);
    
    print $dfh qq[TEST: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
    
    if( ! $genotypeMode ) #calling mode
    {
        #want to be more descriptive with filter values to return e.g. total reads too low, one-side ok only, sides OK but distance too large
        #my $lhsRatioPass = ( $lhsRevBlue == 0 ) || ( $lhsRevBlue > 0 && $lhsFwdBlue / $lhsRevBlue > 2 );
        #my $rhsRatioPass = ( $rhsFwdBlue == 0 ) || ( $rhsFwdBlue > 0 && $rhsRevBlue / $rhsFwdBlue > 2 );
        my $lhsRatioPass = ( $lhsRev == 0 ) || ( $lhsRev > 0 && $lhsFwd / $lhsRev > 2 );
        my $rhsRatioPass = ( $rhsFwd == 0 ) || ( $rhsFwd > 0 && $rhsRev / $rhsFwd > 2 );
        
        if( ( $lhsFwd + $rhsRev ) < $minReads )
        {
            return [$NOT_ENOUGH_READS_CLUSTER, $callString, 0];
        }
        elsif( $lhsFwd < ($minReads/2) || $rhsRev < ($minReads/2) )
        {
            print qq[Code: $NOT_ENOUGH_READS_FLANKS\t$lhsFwd\t$rhsRev\n];
            return [$NOT_ENOUGH_READS_FLANKS, $callString, 0];
        }
=pod
        elsif( $lhsFwdBlue < $minBlue || $rhsRevBlue < $minBlue )
        {
            return [$NOT_ENOUGH_READS_BLUE, $callString, 0];
        }
=cut
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

sub checkBreakpointsSR
{
    croak qq[Incorrect number of arguments: ].scalar(@_) unless @_ == 6;
    my $calls = shift;
    my $bref = shift;my @bams = @{$bref};
    my $minReads = shift;
    my $ignoreRGs = shift;
    my $minMapQ = shift;
    my $outputF = shift;
    
    my $removed = 0;
    open( my $ifh, $calls ) or die qq[Failed to open calls file $calls: $!];
    open( my $ofh, qq[>$outputF] ) or die qq[Failed to create output calls file $outputF: $!];
    #go through each call and check there isnt a cluster of FF or RR reads around it (i.e. indicating its an inversion breakpoint)
    while( my $call = <$ifh> )
    {
        chomp($call);
        my @originalCallA = split( /\t/, $call );
        my $chr = $originalCallA[ 0 ];
        my $refPos = $originalCallA[ 2 ] == $originalCallA[ 1 ] ? $originalCallA[ 1 ] : $originalCallA[ 1 ] + int( ( $originalCallA[ 2 ] - $originalCallA[ 1 ] ) / 2 );
        my $lhsWindow = $refPos - $originalCallA[ 1 ]; $lhsWindow = 400 if( $lhsWindow < 400 );
        my $rhsWindow = $originalCallA[ 2 ] - $refPos; $rhsWindow = 400 if( $rhsWindow < 400 );
        
        my $cmdpre;
        if( @bams > 1 )
        {
            if( _mergeRegionBAMs( \@bams, $chr, $refPos-$lhsWindow-400, $refPos+$rhsWindow+400, qq[/tmp/$$.region.bam] ) )
            {
                $cmdpre = qq[samtools view /tmp/$$.region.bam $chr:];
            }else{die qq[Failed to extract region $chr:$refPos BAM];}
        }
        else
        {
            $cmdpre = qq[samtools view $bams[0] $chr:];
        }
        
        my %numSameOrientation;
        my $cmd = $cmdpre.($refPos-$lhsWindow).qq[-].($refPos+$rhsWindow).qq[ | ].(defined($ignoreRGs) ? qq[ grep -v -f $ignoreRGs |] : qq[]);
        open( my $tfh, $cmd ) or die $!;
        while( my $sam = <$tfh> )
        {
            chomp( $sam );
            my @s = split( /\t/, $sam );
            next unless $s[ 4 ] > $minMapQ;
            
            #       read is mapped                                  mate is mapped                           not paired correctly                  ins size < 50k               both mates are mapped to same strand
            if( ! ($s[ 1 ] & $$BAMFLAGS{'unmapped'} ) && !($s[ 1 ] & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $s[ 1 ] & $$BAMFLAGS{'read_paired'} ) && $s[ 8 ] < 50000 && ( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) == ( $s[ 1 ] & $$BAMFLAGS{'mate_reverse'} ) )
            {
                $numSameOrientation{ $s[ 0 ] } = 1;
                print qq[OR: $s[0]\t$s[1]\n];
            }
        }
        
        unlink( qq[/tmp/$$.region.bam] );
        
        if( keys( %numSameOrientation ) < ( $originalCallA[ 4 ] * 0.5 ) )
        {
            print $ofh qq[$call\n];
        }
        else{my $num = scalar( keys( %numSameOrientation ) );print qq[Removing call $call. Same orientation: $num\n];$removed++;}
    }
    
    return $removed;
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
        last if( ! defined( $fwdCount{$i} ) || ! defined( $revCount{$i} ) );
        my $sum = $fwdCount{$i}+$revCount{$i};
        #print qq[$i\t$sum\n];
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
    
    if( $mpileupOut )
    {
        @s = split( /\t/, $mpileupOut );
        return [$maxPos, length($s[4]) ]; #also need to return the size of the window to the two clusters either side of the breakpoint
    }
    else
    {
        return [$maxPos,5000];
    }
}

#cluster the reads into fwd clusters and rev clusters
sub convertToRegionBedPairs
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 8;
    
    my $bedIn = shift;
	my $minReads = shift;
	my $id = shift;
	my $max_intra_gap = shift;
	my $max_fwd_rev_gap = shift;
	my $doDupRemoval = shift; #1 to not count anchors starting at same pos, 0 to use all anchors
	my $allowOneSided = shift; #allow calls that only have support from one side
	my $outputbed = shift;
	
	my %fwdClusters; #{chr} {endpos} = bed
	my %revClusters; #{chr} {startpos} = bed
	
	#dup removal - uniq fwd start pos, uniq rev end pos
	my %startPosFwd;
	my %endPosRev;
	
	open( my $cfh, $bedIn ) or die $!;
	my $regionStartFwd = 0;
	my $regionEndFwd = 0;
	my $regionChrFwd = undef;
	my $reads_in_regionFwd = 0;
	my $reads_with_idFwd = 0; #number of reads that align to the TE element being called here (i.e. the input could be a mixture of unknown and reads that align to the TE type being called)
	
	my $regionStartRev = 0;
	my $regionEndRev = 0;
	my $regionChrRev = undef;
	my $reads_in_regionRev = 0;
	my $reads_with_idRev = 0;
	while( my $line = <$cfh> )
	{
	    chomp( $line );
	    my @s = split( /\t/, $line ); #chr, start, stop, name, orientation, mate orientation
	    if( $s[ 4 ] eq '+' )
	    {
	        if( ! defined $regionChrFwd )
            {
                $regionStartFwd = $s[ 1 ];
                $regionEndFwd = $s[ 2 ];
                $regionChrFwd = $s[ 0 ];
                $reads_in_regionFwd = 1;
                $startPosFwd{ $s[ 1 ] } = 1;
                $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 0 ] ne $regionChrFwd ) #new chr
            {
                #done - call the region
                my $size = $regionEndFwd - $regionStartFwd; $size = 1 unless $size > 0;
                $fwdClusters{ $regionChrFwd }{ $regionEndFwd } = [ $regionChrFwd, $regionStartFwd, $regionEndFwd, $reads_in_regionFwd ] if( $reads_in_regionFwd >= $minReads && $reads_with_idFwd >= ( $minReads * 0.5 ) );
                
                $reads_in_regionFwd = 1;
                $regionStartFwd = $s[ 1 ];
                $regionEndFwd = $s[ 2 ];
                $regionChrFwd = $s[ 0 ];
                
                %startPosFwd = ();
                $reads_with_idFwd = 0;
                $startPosFwd{ $s[ 1 ] } = 1;
                $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 1 ] - $regionEndFwd > $max_intra_gap )
            {
                #call the region
print "POS $regionChrFwd\t$regionStartFwd\t$regionEndFwd\t$id\t$reads_in_regionFwd\t$reads_with_idFwd\n" if( $reads_in_regionFwd >= $minReads );
                $fwdClusters{ $regionChrFwd }{ $regionEndFwd } = [ $regionChrFwd, $regionStartFwd, $regionEndFwd, $reads_in_regionFwd ] if( $reads_in_regionFwd >= $minReads && $reads_with_idFwd >= ( $minReads * 0.5 ) );
                
                $reads_in_regionFwd = 1;
                $regionStartFwd = $s[ 1 ];
                $regionChrFwd = $s[ 0 ];
                $regionEndFwd = $s[ 2 ];
                
                %startPosFwd = ();
                $reads_with_idFwd = 0;
                $startPosFwd{ $s[ 1 ] } = 1;
                $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
            }
            else
            {
                #read is within the region - increment
                if( $doDupRemoval == 0 || ( $doDupRemoval == 1 && ! defined( $startPosFwd{ $s[ 1 ] } ) ) )
                {
                    $reads_in_regionFwd ++;
                    $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
                    $regionEndFwd = $s[ 2 ];
                    $startPosFwd{ $s[ 1 ] } = 1;
                }
            }
	    }else
	    {
	        if( ! defined $regionChrRev )
            {
                $regionStartRev = $s[ 1 ];
                $regionEndRev = $s[ 2 ];
                $regionChrRev = $s[ 0 ];
                $reads_in_regionRev = 1;
                $endPosRev{ $s[ 2 ] } = 1;
                $reads_with_idRev ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 0 ] ne $regionChrRev ) #on next chr
            {
                #done - call the region
                $revClusters{ $regionChrRev }{ $regionStartRev } = [ $regionChrRev, $regionStartRev, $regionEndRev, $reads_in_regionRev ] if( $reads_in_regionRev >= $minReads && $reads_with_idRev >= ( $minReads * 0.5 ) );
                
                $reads_in_regionRev = 1;
                $regionStartRev = $s[ 1 ];
                $regionEndRev = $s[ 2 ];
                $regionChrRev = $s[ 0 ];
                
                %endPosRev = ();
                $reads_with_idRev = 0;
                $endPosRev{ $s[ 2 ] } = 1;
                $reads_with_idRev ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 1 ] - $regionEndRev > $max_intra_gap )
            {
                #call the region
print "NEG $regionChrRev\t$regionStartRev\t$regionEndRev\t$id\t$reads_in_regionRev\t$reads_with_idRev\n" if( $reads_in_regionRev >= $minReads );
                $revClusters{ $regionChrRev }{ $regionStartRev } = [ $regionChrRev, $regionStartRev, $regionEndRev, $reads_in_regionRev ] if( $reads_in_regionRev >= $minReads && $reads_with_idRev >= ( $minReads * 0.5 ) );
                
                $reads_in_regionRev = 1;
                $regionStartRev = $s[ 1 ];
                $regionEndRev = $s[ 2 ];
                $regionChrRev = $s[ 0 ];
                
                %endPosRev = ();
                $reads_with_idRev = 0;
                $endPosRev{ $s[ 1 ] } = 1;
            }
            else
            {
                #read is within the region - increment
                if( $doDupRemoval == 0 || ( $doDupRemoval == 1 && ! defined( $endPosRev{ $s[ 2 ] } ) ) )
                {
                    $reads_in_regionRev ++;
                    $reads_with_idRev ++ if( $s[ 3 ] eq $id );
                    $regionEndRev = $s[ 2 ];
                    $endPosRev{ $s[ 2 ] } = 1;
                }
            }
	    }
	}
	print "POS $regionChrFwd\t$regionStartFwd\t$regionEndFwd\t$id\t$reads_in_regionFwd\n" if( $reads_in_regionFwd >= $minReads );
    $fwdClusters{ $regionChrFwd }{ $regionEndFwd } = [ $regionChrFwd, $regionStartFwd, $regionEndFwd, $reads_in_regionFwd ] if( $reads_in_regionFwd >= $minReads && $reads_with_idFwd >= ( $minReads * 0.5 ) );
	print "NEG $regionChrRev\t$regionStartRev\t$regionEndRev\t$id\t$reads_in_regionRev\n" if( $reads_in_regionRev >= $minReads );
    $revClusters{ $regionChrRev }{ $regionStartRev } = [ $regionChrRev, $regionStartRev, $regionEndRev, $reads_in_regionRev ] if( $reads_in_regionRev >= $minReads && $reads_with_idRev >= ( $minReads * 0.5 ) );
    
	close( $cfh );
	
	if( ! %fwdClusters || ! %revClusters ){return 0;}
	
	my $regionsCalled = 0;
	#now pair up the clusters by closest end/start positions
	open( my $ofh, qq[>$outputbed] ) or die $!;
	foreach my $fwdChr( keys( %fwdClusters ) )
	{
	    my %fwdChrClusters = %{$fwdClusters{$fwdChr}};
	    
	    next if( ! $revClusters{$fwdChr} );
	    my %revChrClusters = %{$revClusters{$fwdChr}};
	    
	    my @revStartPositions = sort( {$a<=>$b} keys( %{$revClusters{$fwdChr}} ) );
	    my $revStartPositionsIndex = 0;
	    my @fwdEndPositions = sort( {$a<=>$b} keys( %fwdChrClusters ) );
	    my $fwdEndPositionsIndex = 0;
	    while( $revStartPositionsIndex < @revStartPositions && $fwdEndPositionsIndex < @fwdEndPositions )
	    {
print qq[TEST: $fwdEndPositions[ $fwdEndPositionsIndex ] $fwdEndPositionsIndex $revStartPositions[ $revStartPositionsIndex ] $revStartPositionsIndex\n];
print abs( $revStartPositions[ $revStartPositionsIndex ] - $fwdEndPositions[ $fwdEndPositionsIndex ] ).qq[\n];
            #       if end of fwd is close enough to start of rev cluster
	        if( abs( $revStartPositions[ $revStartPositionsIndex ] - $fwdEndPositions[ $fwdEndPositionsIndex ] ) < $max_fwd_rev_gap
	            ||
	            #       forward cluster contained within reverse cluster
	            $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ] >= $revStartPositions[ $revStartPositionsIndex ] && $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ] <= $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ]
	            ||
	            #reverse cluster contained within forward cluster
	            $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ] <= $revStartPositions[ $revStartPositionsIndex ] && $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ] >= $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ]
            )
	        {
	            print qq[Cluster pair: ].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ].qq[\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ].qq[\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ].qq[\n];
	            
                my $reads = $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 3 ] + $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 3 ];
                if( $reads >= $minReads * 2 )
                {
                    #work out 5' and 3' boundaries of the cluster
                    my $start = $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ] < $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 1 ] ? $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ] : $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 1 ];
                    my $end = $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ] > $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ] ? $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ] : $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ];
                    
                    $regionsCalled ++;
                    print $ofh qq[$fwdChr\t$start\t$end\t$id\t$reads\t].$PASS.qq[\n];
                    print qq[CLUSTER: $fwdChr\t$start\t$end\t$id\t$reads\t].$PASS.qq[\n];
                    
                    #print $ofh qq[$fwdChr\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    #print qq[CLUSTER $fwdChr\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    
                    undef( $fwdEndPositions[ $fwdEndPositionsIndex ] );
                    undef( $revStartPositions[ $revStartPositionsIndex ] );
                }
                $revStartPositionsIndex ++;
                $fwdEndPositionsIndex ++;
                next;
	        }
=pod	        
	        elsif( $allowOneSided == 1 )
	        {
	            if($fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 3 ] > $minReads * 2 )
	            {
	                $regionsCalled ++;
	                my $reads = $fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 3 ];
	                print $ofh qq[$fwdChr\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ].qq[\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    print qq[CLUSTER_F $fwdChr\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 1 ].qq[\t].$fwdClusters{ $fwdChr }{ $fwdEndPositions[ $fwdEndPositionsIndex ] }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    $fwdEndPositionsIndex ++;
	            }
	            elsif( $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 3 ] > $minReads * 2 )
	            {
	                $regionsCalled ++;
	                my $reads = $revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 3 ];
	                print $ofh qq[$fwdChr\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    print qq[CLUSTER_R $fwdChr\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revStartPositions[ $revStartPositionsIndex ] }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
	                $revStartPositionsIndex ++;
	            }
	        }
=cut	        
	        if( $revStartPositions[ $revStartPositionsIndex ] - $fwdEndPositions[ $fwdEndPositionsIndex ] > $max_fwd_rev_gap )
	        {
	            $fwdEndPositionsIndex ++;
	        }
	        else{$revStartPositionsIndex ++;}
	    }
	    
	    if( $allowOneSided == 1 )
	    {
#=pod	    
            foreach my $fwdClusterEndPos( @fwdEndPositions )
            {
                if( $fwdClusterEndPos && $fwdClusters{ $fwdChr }{ $fwdClusterEndPos }[ 3 ] > $minReads )
                {
                    my $reads = $fwdClusters{ $fwdChr }{ $fwdClusterEndPos }[ 3 ];
                    print $ofh qq[$fwdChr\t].$fwdClusters{ $fwdChr }{ $fwdClusterEndPos }[ 1 ].qq[\t].$fwdClusters{ $fwdChr }{ $fwdClusterEndPos }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    print qq[CLUSTER_F $fwdChr\t].$fwdClusters{ $fwdChr }{ $fwdClusterEndPos }[ 1 ].qq[\t].$fwdClusters{ $fwdChr }{ $fwdClusterEndPos }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    $regionsCalled ++;
                }
            }
            
            foreach my $revClusterEndPos( @revStartPositions )
            {
                if( $revClusterEndPos && $revClusters{ $fwdChr }{ $revClusterEndPos }[ 3 ] > $minReads )
                {
                    my $reads = $revClusters{ $fwdChr }{ $revClusterEndPos }[ 3 ];
                    print $ofh qq[$fwdChr\t].$revClusters{ $fwdChr }{ $revClusterEndPos }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revClusterEndPos }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    print qq[CLUSTER_R $fwdChr\t].$revClusters{ $fwdChr }{ $revClusterEndPos }[ 1 ].qq[\t].$revClusters{ $fwdChr }{ $revClusterEndPos }[ 2 ].qq[\t$id\t$reads\t].$PASS.qq[\n];
                    $regionsCalled ++;
                }
            }
#=cut
        }
	}
	
	close( $ofh );
	
	return $regionsCalled;
}

#cluster the reads into fwd clusters and rev clusters
sub convertToRegionBedPairsWindowBED
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 8;
    
    my $bedIn = shift;
	my $minReads = shift;
	my $id = shift;
	my $max_intra_gap = shift;
	my $max_fwd_rev_gap = shift;
	my $doDupRemoval = shift; #1 to not count anchors starting at same pos, 0 to use all anchors
	my $allowOneSided = shift; #allow calls that only have support from one side
	my $outputbed = shift;
	
	#dup removal - uniq fwd start pos, uniq rev end pos
	my %startPosFwd;
	my %endPosRev;
	
	open( my $cfh, $bedIn ) or die $!;
	my $regionStartFwd = 0;
	my $regionEndFwd = 0;
	my $regionChrFwd = undef;
	my $reads_in_regionFwd = 0;
	my $reads_with_idFwd = 0; #number of reads that align to the TE element being called here (i.e. the input could be a mixture of unknown and reads that align to the TE type being called)
	
	my $regionStartRev = 0;
	my $regionEndRev = 0;
	my $regionChrRev = undef;
	my $reads_in_regionRev = 0;
	my $reads_with_idRev = 0;
	open( my $posfh, qq[>$$.$id.pos.bed] ) or die qq[failed to create temp POS file\n];
	open( my $negfh, qq[>$$.$id.neg.bed] ) or die qq[failed to create temp NEG file\n];
	while( my $line = <$cfh> )
	{
	    chomp( $line );
	    my @s = split( /\t/, $line ); #chr, start, stop, name, orientation, mate orientation
	    if( $s[ 4 ] eq '+' )
	    {
	        if( ! defined $regionChrFwd )
            {
                $regionStartFwd = $s[ 1 ];
                $regionEndFwd = $s[ 2 ];
                $regionChrFwd = $s[ 0 ];
                $reads_in_regionFwd = 1;
                $startPosFwd{ $s[ 1 ] } = 1;
                $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 0 ] ne $regionChrFwd ) #new chr
            {
                #done - call the region
                my $size = $regionEndFwd - $regionStartFwd; $size = 1 unless $size > 0;
                if( $reads_in_regionFwd >= $minReads )
                {
                    if( $reads_with_idFwd >= ( $minReads * 0.5 ) )
                    {
                        print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\t$id\t+\t$reads_in_regionFwd\n];
                    }
                    else
                    {
                        print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\tunknown\t+\t$reads_in_regionFwd\n];
                    }
                }
                
                $reads_in_regionFwd = 1;
                $regionStartFwd = $s[ 1 ];
                $regionEndFwd = $s[ 2 ];
                $regionChrFwd = $s[ 0 ];
                
                %startPosFwd = ();
                $reads_with_idFwd = 0;
                $startPosFwd{ $s[ 1 ] } = 1;
                $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 1 ] - $regionEndFwd > $max_intra_gap )
            {
                #call the region
                if( $reads_in_regionFwd >= $minReads )
                {
                    if( $reads_with_idFwd >= ( $minReads * 0.5 ) )
                    {
                        print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\t$id\t+\t$reads_in_regionFwd\n];
                    }
                    else
                    {
                        print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\tunknown\t+\t$reads_in_regionFwd\n]
                    }
                }
                
                $reads_in_regionFwd = 1;
                $regionStartFwd = $s[ 1 ];
                $regionChrFwd = $s[ 0 ];
                $regionEndFwd = $s[ 2 ];
                
                %startPosFwd = ();
                $reads_with_idFwd = 0;
                $startPosFwd{ $s[ 1 ] } = 1;
                $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
            }
            else
            {
                #read is within the region - increment
                if( $doDupRemoval == 0 || ( $doDupRemoval == 1 && ! defined( $startPosFwd{ $s[ 1 ] } ) ) )
                {
                    $reads_in_regionFwd ++;
                    $reads_with_idFwd ++ if( $s[ 3 ] eq $id );
                    $regionEndFwd = $s[ 2 ];
                    $startPosFwd{ $s[ 1 ] } = 1;
                }
            }
	    }else
	    {
	        if( ! defined $regionChrRev )
            {
                $regionStartRev = $s[ 1 ];
                $regionEndRev = $s[ 2 ];
                $regionChrRev = $s[ 0 ];
                $reads_in_regionRev = 1;
                $endPosRev{ $s[ 2 ] } = 1;
                $reads_with_idRev ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 0 ] ne $regionChrRev ) #on next chr
            {
                #done - call the region
                if( $reads_in_regionRev >= $minReads )
                {
                    if( $reads_with_idRev >= ( $minReads * 0.5 ) )
                    {
                        print $negfh qq[$regionChrRev\t$regionStartRev\t$regionEndRev\t$id\t-\t$reads_in_regionRev\n]
                    }
                    else
                    {
                        print $negfh qq[$regionChrRev\t$regionStartRev\t$regionEndRev\tunknown\t-\t$reads_in_regionRev\n]
                    }
                }
                $reads_in_regionRev = 1;
                $regionStartRev = $s[ 1 ];
                $regionEndRev = $s[ 2 ];
                $regionChrRev = $s[ 0 ];
                
                %endPosRev = ();
                $reads_with_idRev = 0;
                $endPosRev{ $s[ 2 ] } = 1;
                $reads_with_idRev ++ if( $s[ 3 ] eq $id );
            }
            elsif( $s[ 1 ] - $regionEndRev > $max_intra_gap )
            {
                #call the region
                if( $reads_in_regionRev >= $minReads )
                {
                    if( $reads_with_idRev >= ( $minReads * 0.5 ) )
                    {
                        print $negfh qq[$regionChrRev\t$regionStartRev\t$regionEndRev\t$id\t-\t$reads_in_regionRev\n];
                    }
                    else
                    {
                        print $negfh qq[$regionChrRev\t$regionStartRev\t$regionEndRev\tunknown\t-\t$reads_in_regionRev\n];
                    }
                }
                $reads_in_regionRev = 1;
                $regionStartRev = $s[ 1 ];
                $regionEndRev = $s[ 2 ];
                $regionChrRev = $s[ 0 ];
                
                %endPosRev = ();
                $reads_with_idRev = 0;
                $endPosRev{ $s[ 1 ] } = 1;
            }
            else
            {
                #read is within the region - increment
                if( $doDupRemoval == 0 || ( $doDupRemoval == 1 && ! defined( $endPosRev{ $s[ 2 ] } ) ) )
                {
                    $reads_in_regionRev ++;
                    $reads_with_idRev ++ if( $s[ 3 ] eq $id );
                    $regionEndRev = $s[ 2 ];
                    $endPosRev{ $s[ 2 ] } = 1;
                }
            }
	    }
	}
	
	if( $reads_in_regionFwd >= $minReads )
	{
	    if( $reads_with_idFwd >= ( $minReads * 0.5 ) )
	    {
	        print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\t$id\t+\t$reads_in_regionFwd\n]
	    }
	    else{print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\tunknown\t+\t$reads_in_regionFwd\n]}
	}
	if( $reads_in_regionRev >= $minReads )
	{
	    if( $reads_with_idRev >= ( $minReads * 0.5 ) )
	    {
	        print $negfh qq[$regionChrRev\t$regionStartRev\t$regionEndRev\t$id\t-\t$reads_in_regionRev\n];
	    }
	    else{print $negfh qq[$regionChrRev\t$regionStartRev\t$regionEndRev\tunknown\t-\t$reads_in_regionRev\n];}
	}
    
	close( $posfh );close( $negfh );
	
	my $updownWindow = $max_fwd_rev_gap / 2;
	system( qq[bedtools window -a $$.$id.pos.bed -b $$.$id.neg.bed -l $updownWindow -r $updownWindow > $$.$id.wb.out] ) == 0 or die qq[bedtools window failed\n];
	
	#now pair up the clusters by closest end/start positions - use windowBED to cluster
	open( my $ofh, qq[>$outputbed] ) or die $!;
	open( my $ifh, qq[$$.$id.wb.out] ) or die qq[failed to read bedtools window output: $!];
	my $regionsCalled = 0;
	while( my $line = <$ifh> )
	{
	    my @s = split( /\t/, $line );
	    
	    #check if at least one of the clusters is from the TE type
	    if( $s[ 3 ] ne $id && $s[ 9 ] ne $id )
	    {
	        print qq[Skipping cluster: unknown type\n];
	        next;
	    }
	    
	    #do some sanity checks??
	    my $totalReads = $s[ 5 ] + $s[ 11 ];
	    my $end = $s[2]>$s[8]?$s[2]:$s[8];
	    my $start = $s[1]<$s[7]?$s[1]:$s[7];
	    
	    if( $s[ 3 ] ne $id || $s[ 9 ] ne $id )
	    {
	        print qq[HYBRID: $s[0]\t$start\t$end\t$id\t$totalReads\t$PASS\n];
	    }
	    print $ofh qq[$s[0]\t$start\t$end\t$id\t$totalReads\t$PASS\n];
	    $regionsCalled ++;
	}
	close( $ofh );
	
	return $regionsCalled;
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
        print qq[WARNING: Cant determine sample name from BAM - setting to file name\n];
        $sampleName = basename( $bams[ 0 ] );
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
                return;
            }
        }
        close( $bfh );
        die qq[ERROR: Cant determine version number of $binary\n];
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

sub calculateCigarBreakpoint
{
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 2;
    
    my $sam_pos = shift;
    my $cigar = shift;
    
    #clip is 5'
    if( $cigar =~ /^\d+S/ )
    {
        return $sam_pos;
    }
    #clip is 3'
    else
    {
        if( $cigar =~ /(\d\d)M/ ){return ($sam_pos + $1);}else{return $sam_pos;}
    }
}

=pod
return 1 - unmapped mate
return 2 - mate mapped but incorrectly
=cut
sub isSupportingClusterRead
{
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 6;
    
    my $flag = shift;
    my $insert = shift;
    my $mapQ = shift;
    my $minQual = shift;
    my $minSoftClip = shift;
    my $cigar = shift;
    
    #            read is not a duplicate        map quality is >= minimum
    if( ! ( $flag & $$BAMFLAGS{'duplicate'} ) && $mapQ >= $minQual )
    {
        #           read is mapped                       mate is unmapped
        if( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ( $flag & $$BAMFLAGS{'mate_unmapped'} ) )
        {
            return 1;
        }
        #                               read mapped                             mate mapped                            ins size sensible          has single soft clip in cigar bigger than minimum clip size
        elsif( ( defined( $minSoftClip ) && $minSoftClip > 0 && ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ( $flag & $$BAMFLAGS{'read_paired'} ) && abs( $insert ) < 3000 && ($cigar=~tr/S/S/) == 1 && $cigar =~ /(\d+)(S)/ && $1 > $minSoftClip )
            ||
            #            read is mapped                         mate is mapped                                  not paired correctly                has large enough deviation from expected ins size (short libs assumed)
            ( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $flag & $$BAMFLAGS{'read_paired'} ) && ( abs( $insert ) > 30000 || $insert == 0 ) )
        )
        {
            return 2;
        }
    }
    return 0;
}

sub _removeDups
{
    my $input = shift;
    my $output = shift;
    
    open( my $ifh, $input ) or die $!;
    open( my $ofh, qq[>$output] ) or die $!;
    my $currentStart = -1;
    my $currentScore = -1;
    my $currentBest;
    while( my $l = <$ifh> )
    {
        chomp( $l );
        my @s = split( /\t/, $l );
        
        if( $currentStart == -1 || abs( $s[ 1 ] - $currentStart ) > 200 )
        {
            if( $currentStart != -1 )
            {
                print $ofh qq[$currentBest\n];
            }
            
            $currentStart = $s[ 2 ];
            $currentBest = $l;
            $currentScore = $s[ 4 ];
        }
        elsif( $s[ 4 ] > $currentScore )
        {
            $currentBest = $l;
            $currentScore = $s[ 4 ];
        }
    }
    print $ofh qq[$currentBest\n] if( defined( $currentBest ) );
    close( $ifh );
    close( $ofh );
    
    return 1;
}

1;
