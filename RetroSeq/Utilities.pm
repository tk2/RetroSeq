package RetroSeq::Utilities;
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
    'secondary'    => 0x0100,
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

our $VERSION = 1.5;

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
    die qq[Non integer value passed as refPos: $refPos] unless defined( $refPos ) && $refPos =~ /^\d+$/;
	my @bams = @{ $_[ 0 ] };shift;
    my $minMapQ = shift;
    my $originalCall = shift;
    my $dfh = shift; #file handle to print out info on failed calls
    my $ignoreRGs = shift; #file of lines with "RG:tag\t"
    my $minReads = shift;
    my $genotypeMode = shift; #0/1 saying whether to operate in genotyping mode (i.e. less stringent criteria)
    
    my @originalCallA = split( /\t/, $originalCall );
    my $lhsWindow = $refPos - $originalCallA[ 1 ]; $lhsWindow = 400 if( $lhsWindow < 400 );
    my $rhsWindow = $originalCallA[ 2 ] - $refPos; $rhsWindow = 400 if( $rhsWindow < 400 );
    
    #test to see if lots of rp's either side
    my $lhsFwdBlue = 0; my $lhsRevBlue = 0; my $rhsFwdBlue = 0; my $rhsRevBlue = 0;
    my $lhsFwdGreen = 0; my $lhsRevGreen = 0; my $rhsFwdGreen = 0; my $rhsRevGreen = 0;
    my $softClipSupporting = 0;
	
    #store the last blue read before the b/point, and first blue read after the b/point
	my $lastSupportingPos = 0;my $firstSupportingPos = 100000000000;
	
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
	my @lhsFwd;my @lhsRev;my @rhsFwd;my @rhsRev;my @soft;
	my $totalSpanningRPs = 0;my %spanningFrags;
	while( my $sam = <$tfh> )
	{
	    chomp( $sam );
	    my @s = split( /\t/, $sam );
	    my $supporting = isSupportingClusterRead( $s[ 1 ], $s[ 8 ], $s[ 4 ], $minMapQ, $s[ 5 ], $refPos, $s[ 3 ], length( $s[9] ) );
		
	    if( $supporting  == 2 )
	    {
			if( ( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) )  #rev strand
            {
				if( $s[ 5 ] =~ /^[0-9]{2}+S/ ){$rhsRevBlue++;push(@rhsRev, $s[0]);}
                elsif( $s[ 3 ] < $refPos ){$lhsRevBlue++;push(@lhsRev, $s[0]);}else{$rhsRevBlue++;$firstSupportingPos = $s[ 3 ] if( $s[ 3 ] < $firstSupportingPos );push(@rhsRev, $s[0]);}
            }
            else
            {
				if( $s[ 5 ] =~ /[0-9]{2}+S$/ ){$lhsFwdBlue++;push(@lhsFwd, $s[0]);}
                elsif( $s[ 3 ] < $refPos ){$lhsFwdBlue++;$lastSupportingPos = $s[ 3 ] + length( $s[ 9 ] ) if( ( $s[ 3 ] + length( $s[ 9 ] ) ) > $lastSupportingPos );push(@lhsFwd, $s[0]);}else{$rhsFwdBlue++;push(@rhsFwd, $s[0]);}
            }
        }
        #        the mate is unmapped
        elsif( $supporting == 1 )
        {
            if( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) #rev strand
            {
				if( $s[ 5 ] =~ /^[0-9]{2}+S/ ){$rhsRevGreen++;push(@rhsRev, $s[0]);}
                elsif( $s[ 3 ] < $refPos ){$lhsRevGreen++;push(@lhsRev, $s[0]);}else{$rhsRevGreen++;$firstSupportingPos = $s[ 3 ] if( $s[ 3 ] < $firstSupportingPos );push(@rhsRev, $s[0]);}
            }
            else
            {
				if( $s[ 5 ] =~ /[0-9]{2}+S$/ ){$lhsFwdGreen++;push(@lhsFwd, $s[0]);}
                elsif( $s[ 3 ] < $refPos ){$lhsFwdGreen++;$lastSupportingPos = $s[ 3 ] + length( $s[ 9 ] ) if( ( $s[ 3 ] + length( $s[ 9 ] ) ) > $lastSupportingPos );push(@lhsFwd, $s[0]);}else{$rhsFwdGreen++;push(@rhsFwd, $s[0]);}
            }
        }
        
        #       read is mapped                                  mate is mapped                           not paired correctly                  ins size < 50k               both mates are mapped to same strand
        if( ! ($s[ 1 ] & $$BAMFLAGS{'unmapped'} ) && !($s[ 1 ] & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $s[ 1 ] & $$BAMFLAGS{'read_paired'} ) && abs($s[ 8 ]) < 50000 && ( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) == ( $s[ 1 ] & $$BAMFLAGS{'mate_reverse'} ) && ( $s[ 2 ] eq $s[ 6 ] || $s[ 6 ] eq '=') ) #and both mates mapped to the same chr
        {
            $numSameOrientation ++;
        }
        
        #determine outer co-ordinate of read (ex clipping)
        if( $s[ 1 ] & $$BAMFLAGS{'paired_tech'} && !($s[ 1 ] & $$BAMFLAGS{'unmapped'} ) && !($s[1] & $$BAMFLAGS{'secondary'}) )
        {
            if( ( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) )
            {
                my $outer = getSAMendpos($s[3],$s[5]);
                $spanningFrags{$s[0]}{threeprime} = $outer if($outer > $refPos+50);
            }
            else
            {
                $spanningFrags{$s[0]}{fiveprime} = $s[3] if($s[3]<$refPos-50);
            }
        }
    }
    unlink( qq[/tmp/$$.region.bam] );
    foreach my $pair(keys(%spanningFrags)){if($spanningFrags{$pair}{threeprime}&&$spanningFrags{$pair}{fiveprime}){$totalSpanningRPs++; print qq[SPAN: $pair\n];};}print qq[\n];

    #if it appears to be an inversion breakpoint
    my $lhsRev = $lhsRevGreen + $lhsRevBlue;my $rhsRev = $rhsRevGreen + $rhsRevBlue;my $lhsFwd = $lhsFwdGreen + $lhsFwdBlue;my $rhsFwd = $rhsFwdGreen + $rhsFwdBlue;
    my $totalSupporting = ( $lhsFwd + $rhsRev + $softClipSupporting );
    my $callString = qq[$chr\t$refPos\t].($refPos+1).qq[\t$originalCallA[ 3 ]\t$totalSupporting];
    print qq[No same orientation: $numSameOrientation\n];
    if( $numSameOrientation > $totalSupporting )
    {
        return [$INV_BREAKPOINT, $callString, 10000 ];
    }
    
    #check there are supporting read pairs either side of the depth minima
    my $dist = $firstSupportingPos - $lastSupportingPos;
    my $minBlue = int($minReads / 2);
    
    print $dfh qq[TEST: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastSupportingPos\t$firstSupportingPos\t$dist\n];
    
    if( ! $genotypeMode ) #calling mode
    {
        #want to be more descriptive with filter values to return e.g. total reads too low, one-side ok only, sides OK but distance too large
        my $lhsRatioPass = ( $lhsRev == 0 ) || ( $lhsRev > 0 && $lhsFwd / $lhsRev >= 2 );
        my $rhsRatioPass = ( $rhsFwd == 0 ) || ( $rhsFwd > 0 && $rhsRev / $rhsFwd >= 2 ); 
        
        if( $totalSupporting < $minReads )
        {
            return [$NOT_ENOUGH_READS_CLUSTER, $callString, 0, $totalSupporting, $totalSpanningRPs];
        }
        elsif( $lhsFwd < ($minReads/2) || $rhsRev < ($minReads/2) )
        {
            print qq[Code: $NOT_ENOUGH_READS_FLANKS\t$lhsFwd\t$rhsRev\n];
            return [$NOT_ENOUGH_READS_FLANKS, $callString, 0, $totalSupporting, $totalSpanningRPs];
        }
        elsif( ! $lhsRatioPass && ! $rhsRatioPass )
        {
            return [$NEITHER_SIDE_RATIO_PASSES, $callString, 0, $totalSupporting, $totalSpanningRPs];
        }
        elsif( ! ( $lhsRatioPass && $rhsRatioPass ) )
        {
            return [$ONE_SIDED_RATIO_PASSES, $callString, 0, $totalSupporting, $totalSpanningRPs];
        }
        elsif( $dist > 220 )
        {
            print qq[Distance between 5' and 3' clusters: $dist\n];
            return [$DISTANCE_THRESHOLD, $callString, 0, $totalSupporting, $totalSpanningRPs ];
        }
        else
        {
            my $ratio = ( $lhsRev + $rhsFwd ) / ( $lhsFwd + $rhsRev );
            return [$PASS, $callString, $ratio, $totalSupporting, $totalSpanningRPs];
        }
        
        return [$UNKNOWN_FAIL, $callString, 10000, $totalSpanningRPs];
    }
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
    
    my %fwdCount; #pos->read count
    my %revCount; #pos->read count
    my %spanPairs; #readname->f/r->outerpos
    
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
    while( my $sam = <$tfh> )
    {
        chomp( $sam );
        my @samL = split( /\t/, $sam );
        my $flag = $samL[ 1 ];
        
        next if $samL[ 4 ] < $minQ;
        next if ( $flag & $$BAMFLAGS{'duplicate'} );
        
        #       paired technology                       is on fwd strand
        if( ( $flag & $$BAMFLAGS{'paired_tech'} ) && !( $flag & $$BAMFLAGS{ 'reverse_strand' } ) ) #fwd supporting read
        {
            my $endPos = getSAMendpos($samL[3],$samL[5]);
            if( $samL[5]=~/[0-9]+S$/)
            {
                if( defined( $fwdCount{$endPos} ) ){$fwdCount{$endPos}+=2;}else{$fwdCount{$endPos}=2;}
            }
            elsif( !( $flag & $$BAMFLAGS{'read_paired'} ) ) # not paired correctly
            {
                if( defined($fwdCount{$endPos} )){$fwdCount{$endPos}=$fwdCount{$endPos}+1;}else{$fwdCount{$endPos}=1;}
            }
            else{$spanPairs{$samL[0]}{f}=$samL[3];}
        }
        elsif( ( $flag & $$BAMFLAGS{'paired_tech'} ) && ( $flag & $$BAMFLAGS{ 'reverse_strand' } ) ) #rev supporting read
        {
            my $endPos = getSAMendpos($samL[3],$samL[5]);
            if( $samL[5]=~/^[0-9]+S/ )
            {
                if( defined( $revCount{$samL[ 3 ]} ) ){$revCount{$samL[ 3 ]}+=2;}else{$revCount{$samL[ 3 ]}=2;}
			}
            elsif( !( $flag & $$BAMFLAGS{'read_paired'} ) ) # not paired correctly
            {
                if( defined($revCount{$samL[ 3 ]} ) ){$revCount{$samL[3]}++;}else{$revCount{ $samL[ 3 ] }=1;}
            }
            else{$spanPairs{$samL[0]}{r}=$endPos;}
        }
    }
    close( $tfh );

    #now transform the fwdCount and revCount into a decaying profile of cumulative fwd and rev reads
    my $lastUpdatePos = -1;
    $fwdCount{$start}=0;
    for(my $pos=$start+1;$pos<=$end;$pos++)
    {
        if(defined($fwdCount{$pos})){$lastUpdatePos=$pos;$fwdCount{$pos}=$fwdCount{$pos-1}+$fwdCount{$pos};}
        else
        {
            $fwdCount{$pos}=$fwdCount{$pos-1};
            if($lastUpdatePos!=-1&&$fwdCount{$pos}>0&&($pos-$lastUpdatePos)%20==0){$fwdCount{$pos}--;} #decay the cumulative total if 20bp without an increase have passed
        }
    }
    
    #now the revCount decaying profile
    $lastUpdatePos = -1;
    $revCount{$end}=0;
    for(my $pos=$end-1;$pos>=$start;$pos--)
    {
        if(defined($revCount{$pos})){$lastUpdatePos=$pos;$revCount{$pos}=$revCount{$pos+1}+$revCount{$pos};}
        else
        {
            $revCount{$pos}=$revCount{$pos+1};
            if($lastUpdatePos!=-1&&$revCount{$pos}>0&&($lastUpdatePos-$pos)%20==0){$revCount{$pos}--;} #decay the cumulative total if 20bp without an increase have passed
        }
    }
    
    #now find the position where the sum of the two counts maximises
    my $maxPos = $start;my @maxPoss;my $maxVal = $fwdCount{$start}+$revCount{$start};
    
    for(my $i=$start;$i<$end;$i++)
    {
        my $sum = $fwdCount{$i}+$revCount{$i};
        if( $fwdCount{$i}+$revCount{$i} > $maxVal )
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
    if( scalar(@maxPoss)==1){return [$maxPoss[0]];}
    else
    {
        #breakpoint is not precise, therefore use spanning read pairs minimum to select the point
        print qq[Refining imprecise breakpoint: $maxPoss[0] ].$maxPoss[scalar(@maxPoss)-1].qq[\n];
        my %spanningPairs;
        foreach my $pair(keys(%spanPairs))
        {
            next if(! ($spanPairs{$pair}{f} && $spanPairs{$pair}{r}) );
            for(my $i=$spanPairs{$pair}{f};$i<$spanPairs{$pair}{r};$i++){$spanningPairs{$i}++;}
        }
        
        if( $spanningPairs{$maxPoss[0]} )
        {
            my $minSpanning=$spanningPairs{$maxPoss[0]};my $minPosition=$maxPoss[0];
            for(my $i=$maxPoss[0];$i<$maxPoss[scalar(@maxPoss)-1];$i++)
            {  
                if( defined($spanningPairs{$i}) && $spanningPairs{$i}<$minSpanning ){$minSpanning=$spanningPairs{$i};$minPosition=$i;}
            }
            
            return [$minPosition, $maxPoss[0], $maxPoss[scalar(@maxPoss)-1]];
        }
        else{return [int(($maxPoss[0]+$maxPoss[scalar(@maxPoss)-1])/2), $maxPoss[0], $maxPoss[scalar(@maxPoss)-1]];}
    }
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
	        print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\t$id\t+\t$reads_in_regionFwd\n];
	    }
	    else{print $posfh qq[$regionChrFwd\t$regionStartFwd\t$regionEndFwd\tunknown\t+\t$reads_in_regionFwd\n];}
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
	
	if( -s qq[$$.$id.pos.bed] == 0 || -s qq[$$.$id.neg.bed] == 0 )
	{
	    if( -s qq[$$.$id.pos.bed] == 0 && -s qq[$$.$id.neg.bed] > 0 )
	    {
	        system( qq[cat $$.$id.neg.bed >> $$.hybrid.neg.out] ) == 0 or die qq[cat neg failed\n];
	    }
	    elsif( -s qq[$$.$id.pos.bed] > 0 && -s qq[$$.$id.neg.bed] == 0 )
	    {
	        system( qq[cat $$.$id.pos.bed >> $$.hybrid.pos.out] ) == 0 or die qq[cat pos failed\n];
	    }
	    print qq[Zero potential TE call regions identified for $id\n];return 0;
	}
	
	my $updownWindow = $max_fwd_rev_gap / 2;
	system( qq[bedtools window -a $$.$id.pos.bed -b $$.$id.neg.bed -l $updownWindow -r $updownWindow > $$.$id.wb.out] ) == 0 or die qq[bedtools window failed\n];
	
	#collect the unused pos/neg clusters into separate file
	system( qq[bedtools window -a $$.$id.pos.bed -b $$.$id.neg.bed -v -l $updownWindow -r $updownWindow | awk -F"\t" '\$6>$minReads' >> $$.hybrid.pos.out] ) == 0 or die qq[bedtools window failed\n];
	system( qq[bedtools window -a $$.$id.neg.bed -b $$.$id.pos.bed -v -l $updownWindow -r $updownWindow | awk -F"\t" '\$6>$minReads' >> $$.hybrid.neg.out] ) == 0 or die qq[bedtools window failed\n];
	
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
    my $min_version = undef;
    my $max_version = undef;
    if( @_ > 0 ){$min_version = shift;}
    if( @_ > 0 ){$max_version = shift;}
    
    if( ! `which $binary` )
    {
        croak qq[Error: Cant find required binary $binary\n];
    }
    
    if( $min_version )
    {
        my @minv = split( /\./, $min_version );
        open( my $bfh, qq[ $binary 2>&1 | ] ) or die "failed to run $binary\n";
        while(<$bfh>)
        {
            chomp;
            if( lc($_) =~ /version/ && ( $_ =~ /(\d+)\.(\d+)\.(\d+)/ || $_ =~ /(\d+)\.(\d+)/ ) )
            {
                my $ver = $1.'.'.$2;
                if( $ver < $minv[ 0 ].'.'.$minv[ 1 ] )
                {
                    die qq[\nERROR: $binary version at least $min_version is required - your version is $ver\n];
                }
                if( $max_version )
                {
                    my @maxv = split( /\./, $min_version );
                    if( $ver > qq[$maxv[0].$maxv[1]] )
                    {
                        die qq[\nERROR: $binary version incompatible, please use $min_version to $max_version - your version is $ver\n];
                    }
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

    ##FORMAT=<ID=SP,Number=1,Type=Float,Description="Number of read pairs spanning breakpoint, useful for estimation of size of insertion">
    $vcf_out->add_header_line( {key=>'FORMAT', ID=>'SP', Number=>'1', Type=>'Float', Description=>'Number of correctly mapped read pairs spanning breakpoint, useful for estimation of size of insertion'} );
    
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
return 3 - soft clipped read supporting the breakpoint
=cut
sub isSupportingClusterRead
{
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 8;
    
    my $flag = shift;
    my $insert = shift;
    my $mapQ = shift;
    my $minQual = shift;
    my $cigar = shift;
    my $refpos = shift; #if the breakpiont position + read start is provided, then check if the reads soft clip OVER the breakpoint
    my $readpos = shift;
    my $readlength = shift;
    
    return 0 if( !( $flag & $$BAMFLAGS{'paired_tech'} ) );
    
    #            read is not a duplicate        map quality is >= minimum
    if( ! ( $flag & $$BAMFLAGS{'duplicate'} ) && $mapQ >= $minQual )
    {
        #           read is mapped                       mate is unmapped
        if( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ( $flag & $$BAMFLAGS{'mate_unmapped'} ) )
        {
            return 1;
        }
        #            read is mapped                         mate is mapped                                  not paired correctly                has large enough deviation from expected ins size (short libs assumed)
        elsif( ! ( $flag & $$BAMFLAGS{'unmapped'} ) && ! ( $flag & $$BAMFLAGS{'mate_unmapped'} ) && ! ( $flag & $$BAMFLAGS{'read_paired'} ) && ( abs( $insert ) > 30000 || $insert == 0 ) )
        {
            return 2;
        }
    }
    return 0;
}

=pod
return end position of SAM alignment
Input: SAM align start, cigar string
=cut
sub getSAMendpos
{
    die qq[Incorrect number of parameters: ].scalar(@_) unless @_ == 2;
    
    my $startPos = shift;
    my $cigar = shift;
    
    my $endPos=$startPos;
    while($cigar=~/([0-9]+[MIDSHP])/g)
    {
        my $entry=$1;
        if( $entry=~/M$/){$endPos+=substr($entry,0,length($entry)-1);}
    }
    return $endPos;
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
