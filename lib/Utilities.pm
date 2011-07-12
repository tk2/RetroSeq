package Utilities;

use strict;
use warnings;
use Carp;

my $FILTER_WINDOW = 50;
my $BREAKPOINT_WINDOW = 220;

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
    my $bam = shift;
    my $ref = shift;
    my $dfh = shift; #handle to output stream for logging calls that fall out
    
    open( my $tfh, qq[samtools mpileup -r $chr:$start-$end -f $ref -u $bam | bcftools view - | ] ) or die $!;
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
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 6;
    
    my $chr = shift;
    my $refPos = shift;
    die qq[Non integer value passed as refPos] unless defined( $refPos ) && $refPos =~ /^\d+$/;
    my $bam = shift;
    my $minMapQ = shift;
    my $originalCall = shift;
    my $dfh = shift; #file handle to print out info on failed calls
    my @originalCallA = split( /\t/, $originalCall );
    
    #test to see if lots of rp's either side
    my $lhsFwdBlue = 0; my $lhsRevBlue = 0; my $rhsFwdBlue = 0; my $rhsRevBlue = 0;
    my $lhsFwdGreen = 0; my $lhsRevGreen = 0; my $rhsFwdGreen = 0; my $rhsRevGreen = 0;
	
    #store the last blue read before the b/point, and first blue read after the b/point
	my $lastBluePos = 0;my $firstBluePos = 100000000000;
	
	#also check the orientation of the supporting reads (i.e. its not just a random mixture of f/r reads overlapping)
	my $cmd = qq[samtools view $bam $chr:].($refPos-$BREAKPOINT_WINDOW).qq[-].($refPos+$BREAKPOINT_WINDOW).qq[ | ];
	open( my $tfh, $cmd ) or die $!;
	while( my $sam = <$tfh> )
	{
	    chomp( $sam );
	    my @s = split( /\t/, $sam );
	            next unless $s[ 4 ] > $minMapQ;
	            #        the mate is mapped                        not paired correctly                        or mate ref name is different chr
	            if( !($s[ 1 ] & $$BAMFLAGS{'mate_unmapped'}) && ( ( !( $s[ 1 ] & $$BAMFLAGS{'read_paired'} ) ) || ( $s[ 6 ] ne '=' ) ) )
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
	        
	        #check there are supporting read pairs either side of the depth minima
	        my $lhsRev = $lhsRevGreen + $lhsRevBlue;my $rhsRev = $rhsRevGreen + $rhsRevBlue;my $lhsFwd = $lhsFwdGreen + $lhsFwdBlue;my $rhsFwd = $rhsFwdGreen + $rhsFwdBlue;
	        my $dist = $firstBluePos - $lastBluePos;
	        
	        if( $lhsFwdBlue >= 5 && $rhsRevBlue >= 5 && $lhsFwd > 10 && $rhsRev > 10 && ( $lhsRevBlue == 0 || $lhsFwdBlue / $lhsRevBlue > 2 ) && ( $rhsFwdBlue == 0 || $rhsRevBlue / $rhsFwdBlue > 2 ) && $dist < 120 )
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
}

sub genotypeRegion
{
    croak qq[Incorrect number of arguments: ].scalar(@_) unless @_ == 7;
    
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $bam = shift;
    my $minDepth = shift;
    my $minMapQ = shift;
    my $ref = shift;
    
    die qq[cant find bam: $bam\n] unless -f $bam;
    
    #check the min depth in the region is < minDepth
    my $regStart = $start - 200; my $regEnd = $end + 200;
    open( my $tfh, qq[samtools mpileup -f $ref -r $chr:$regStart-$regEnd $bam | tail -300 | head -200 | ] ) or die $!;
    my $minDepth_ = 100; my $minDepthPos = -1;
    while( <$tfh> )
    {
        chomp;my @s = split( /\t/, $_ );
	    my $d = ($s[8]=~tr/,\./x/); #count the depth of the bases that match the reference (i.e. sometimes at the breakpoint there are snps causing false depth)
	    if( $s[ 2 ] eq $s[ 3 ] && $d < $minDepth_ )
	    {
	        $minDepth_ = $d;$minDepthPos = $s[ 1 ];
	    }
	    elsif($s[ 7 ] < $minDepth_ ){$minDepth_ = $s[ 7 ];$minDepthPos = $s[ 1 ];}
    }
    close( $tfh );
    
    if( $minDepth_ > $minDepth ){return undef;} #no call - depth to high
    
    #also check the orientation of the supporting reads (i.e. its not just a random mixture of f/r reads overlapping)
    my $cmd = qq[samtools view $bam $chr:].($minDepthPos-300).qq[-].($minDepthPos+300).qq[ | ];
	open( $tfh, $cmd ) or die $!;
	my $lhsRev = 0; my $rhsRev = 0; my $lhsFwd = 0; my $rhsFwd = 0;
	while( <$tfh> )
	{
	    chomp;
	    my @s = split( /\t/, $_ );
	    next unless $s[ 4 ] > $minMapQ;
	    if( !( $s[ 1 ] & $$BAMFLAGS{'read_paired'} ) && ( $s[ 1 ] & $$BAMFLAGS{'mate_unmapped'} || $s[ 2 ] ne $s[ 6 ] ) )
	    {
	        #print qq[candidate: $s[1]\n];
	        if( $s[ 1 ] & $$BAMFLAGS{'reverse_strand'} ) #rev strand
	        {
	            if( $s[ 3 ] < $minDepthPos ){$lhsRev++;}else{$rhsRev++;}
	        }
	        else
	        {
	            if( $s[ 3 ] < $minDepthPos ){$lhsFwd++;}else{$rhsFwd++;}
	        }
	    }
	 }
	 
	 #N.B. Key difference is that only 1 side is required to have the correct ratio of fw:rev reads
	 if( $lhsFwd > 5 && $rhsRev > 5 && ( ( $lhsRev == 0 || $lhsFwd / $lhsRev > 2 ) || ( $rhsFwd == 0 || $rhsRev / $rhsFwd > 2 ) ) )
	 {
	     return ($lhsFwd+$rhsRev);
	 }
	 return undef;
}

sub getCandidateBreakPointsDir
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 5;
 
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $bam = shift;
    my $minQ = shift;
    
    if( $chr !~ /^\d+$/ || $start !~ /^\d+$/ || $end !~ /^\d+$/ ){die qq[ERROR: Invalid parameters passed to getCandidateBreakPointsDir: $chr $start $end\n];}
    if( ! -f $bam ){die qq[ERROR: Cant find bam file: $bam];}
    
    #get the reads over the region into a bam on local /tmp first
    my $cmd = qq[samtools view $bam -h -b $chr:].($start-1000).qq[-].($end+1000).qq[ > /tmp/$$.region.bam; samtools index /tmp/$$.region.bam];
    if( system( $cmd ) != 0 )
    {
        die qq[WARNING: Failed to extract region from BAM for call $chr:$start-$end];
        return undef;
    }
    
    #get the coverage of the forward orientated reads
    #                   get sam                                 filter out low map Q reads                               mate is unmapped    mate not mapped in proper pair             mapped of fwd strand
    $cmd = qq[samtools view -h /tmp/$$.region.bam $chr:$start-$end | gawk -F"\\t" '(\$1~/^\@/||\$4>=$minQ)'].q[ | gawk -F"\t" '($1~/^@/||and($2,0x0008)||and($2,0x0002)==0)' | gawk -F"\t" '$1~/^@/||and($2,0x0010)==0' | samtools view -bS - > /tmp/].qq[$$.fwd.bam];
    my $fwdpos = undef;
    if( ! system( $cmd ) )
    {
        #call pileup to get the depth profile                                           get the max depth position
        $cmd = qq[samtools mpileup -A /tmp/$$.fwd.bam | sort -k4,4n | tail -1 | ].q[awk -F"\t" '{print $2}' | ];
        open( my $tfh, $cmd ) or warn $!;
        $fwdpos = <$tfh>;
        if( ! defined( $fwdpos ) ){warn qq[WARNING: Failed to get the fwdpos for het call: $chr:$start-$end];return undef;}
        chomp( $fwdpos );
    }
    else{warn qq[Failed to run samtools over region to generate fwd reads bam: $cmd\n];return undef;}
    
    unlink( qq[/tmp/$$.fwd.bam] ) > 0 or warn qq[WARNING: Failed to delete /tmp/$$.fwd.bam files];
    
    if( !$fwdpos || $fwdpos !~ /^\d+$/ ){print qq[WARNING: Failed to estimate fwd heterozygous breakpoints for $chr:$start-$end\n];return undef;}
    
    #get the coverage of the reverse orientated reads
    #                   get sam                                 filter out low map Q reads                               mate is unmapped    mate not mapped in proper pair             mapped of rev strand
    $cmd = qq[samtools view /tmp/$$.region.bam -h $chr:$start-$end | gawk -F"\\t" '(\$1~/^\@/||\$4>=$minQ)'].q[ | gawk -F"\t" '($1~/^@/||and($2,0x0008)||and($2,0x0002)==0)' | gawk -F"\t" '$1~/^@/||and($2,0x0010)' | samtools view -bS - > /tmp/].qq[$$.rev.bam];
    my $revpos = undef;
    if( ! system( $cmd ) )
    {
        $cmd = qq[samtools mpileup -A /tmp/$$.rev.bam | sort -k4,4n | tail -1 | ].q[awk -F"\t" '{print $2}' | ];
        open( my $tfh, $cmd ) or warn $!;
        $revpos = <$tfh>;
        if( ! defined( $revpos ) ){warn qq[WARNING: Failed to get the revpos for het call: $chr:$start-$end];return undef;}
        chomp( $revpos );
    }
    else{warn qq[Failed to run samtools over region to generate rev reads bam: $cmd\n];}
    
    unlink( qq[/tmp/$$.rev.bam] ) > 0 or warn qq[WARNING: Failed to delete /tmp/$$.rev.bam files];
    
    if( !$revpos || $revpos !~ /^\d+$/ ){print qq[WARNING: Failed to estimate rev heterozygous breakpoints for $chr:$start-$end\n];return undef;}
    
    #find the point with the lowest depth between the candidate points
    $cmd = qq[samtools view -h -b /tmp/$$.region.bam -h $chr:].($fwdpos-200).qq[-].($revpos+200).qq[ > /tmp/$$.localised.bam];
    my $lowestDPos = undef;
    if( ! system( $cmd ) )
    {
        $cmd = qq[ samtools mpileup -A /tmp/$$.localised.bam | ].q[ awk -F"\t" '$2>$fwdpos&&$2<$revpos' | sort -k4,4n | head -1 | ].q[awk -F"\t" '{print $2}' | ];
        open( my $tfh, $cmd ) or warn $!;
        $lowestDPos = <$tfh>;
        if( ! defined( $lowestDPos ) ){warn qq[WARNING: Failed to get the lowest depth pos for het call: $chr:$start-$end];}
        else{chomp( $lowestDPos );}
    }else{warn qq[Failed to run samtools over region to generate localised reads bam: $cmd\n];}
    
    unlink( qq[/tmp/$$.region.bam], qq[/tmp/$$.region.bam.bai], qq[/tmp/$$.localised.bam] ) > 0 or warn qq[WARNING: Failed to delete /tmp/$$.region.bam files];
    
    my @positions = (int(( $fwdpos + $revpos ) / 2), int($fwdpos+(($revpos-$fwdpos)*0.25)), int($fwdpos+(($revpos-$fwdpos)*0.75)));
    push( @positions, $lowestDPos ) if( $lowestDPos );
    return \@positions;
}

sub getCandidateBreakPointsDirVote
{
    die qq[ERROR: Incorrect number of arguments supplied: ].scalar(@_) unless @_ == 5;
 
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $bam = shift;
    my $minQ = shift;

    if( $chr !~ /^\d+$/ || $start !~ /^\d+$/ || $end !~ /^\d+$/ ){die qq[ERROR: Invalid parameters passed to getCandidateBreakPointsDir: $chr $start $end\n];}
    if( ! -f $bam ){die qq[ERROR: Cant find bam file: $bam];}
    
    my %fwdCount;
    my %revCount;
    
    #get the coverage of the forward orientated reads
    #                   get sam                                 filter out low map Q reads                               mate is unmapped    mate not mapped in proper pair             mapped of fwd strand
    my $cmd = qq[samtools view $bam $chr:$start-$end |];
    
    open( my $tfh, $cmd ) or die $!;
    my $fwdCurrentPos = $start;
    my $revCurrentPos = $start;
    $fwdCount{ $fwdCurrentPos - 1 } = 0;
    $revCount{ $revCurrentPos - 1 } = 0;
    while( my $sam = <$tfh> )
    {
        chomp( $sam );
        my @samL = split( /\t/, $sam );
        my $flag = $samL[ 1 ];
        
        #update the counts to the current position
        for(my $i=$fwdCurrentPos;$i<$samL[3];$i++ ){$fwdCount{$i}=$fwdCount{$i-1};}
        for(my $i=$revCurrentPos;$i<$samL[3];$i++ ){$revCount{$i}=$revCount{$i-1};}
        
        next if $samL[ 4 ] < $minQ;
        next if ( $flag & $$BAMFLAGS{'duplicate'} );
        #       paired technology                       not paired correctly                    is on fwd strand
        if( ( $flag & $$BAMFLAGS{'paired_tech'} ) && !( $flag & $$BAMFLAGS{'read_paired'} ) && !( $flag & $$BAMFLAGS{ 'reverse_strand' } ) ) #fwd supporting read
        {
            no warnings 'uninitialized'; #perl warns when the value is 0
            $fwdCount{ $samL[ 3 ] } = $fwdCount{ $samL[ 3 ] - 1 } + 1;
            $fwdCurrentPos = $samL[ 3 ] + 1;
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
    
    for(my $i=$start;$i<$revCurrentPos&&$i<$fwdCurrentPos;$i++)
    {
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
    
    $maxPos = $maxPoss[ int(scalar(@maxPoss)/2) ];
    
    return [$maxPos];
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
    
    #read through the supporting reads bed and pick up the reads in the vicinity of the calls
    open( my $sfh, $sortedReadsBED ) or die qq[ERROR: Failed to open sorted reads BED: $sortedReadsBED];
    open( my $cfh, $sortedCallsBED ) or die qq[ERROR: Failed to open calls bed];
    open( my $ofh, qq[>$sortedCallsBED.dir] ) or die qq[ERROR: Failed to create new BED file];
    my $currentCall = <$cfh>;chomp($currentCall);my ($techr,$testart,$teend,$tename,$tescore) = split( /\t/, $currentCall );
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
            elsif( $start > $teend && $end < $teend + $BREAKPOINT_WINDOW )
            {
                if( $orientation eq '+' ){$currentCallRHS++;}else{$currentCallRHS--;}
                $withinCall = 1;
            }
            elsif( $withinCall && $end > $teend + $BREAKPOINT_WINDOW )
            {
                $fwd = 1;
                if( $currentCallLHS > 0 && $currentCallRHS < 0 ){$fwd = 0;} #if it looks like a reverse insertions
                print $ofh qq[$techr\t$testart\$teend\t$tename\t$tescore\t].$fwd ? qq[+\n] : qq[-\n];
                $withinCall = 0;
                $currentCall = <$cfh>;
                last if( ! $currentCall );
                chomp($currentCall);($techr,$testart,$teend,$tename,$tescore) = split( /\t/, $currentCall );
            }
        }
    }
    if( $withinCall ){print $ofh qq[$techr\t$testart\$teend\t$tename\t$tescore\t].$fwd ? qq[+\n] : qq[-\n];}
    
    exit;
    close( $sfh );
}

1;
