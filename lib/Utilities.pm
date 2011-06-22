package Utilities;

use strict;
use warnings;
use Carp;

my $FILTER_WINDOW = 50;
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
	if( !@res || !$res[ 0 ] ){print qq[WARNING: no max/min returned for $chr:$start-$end\n];print $dfh qq[$chr\t$start\t$end\tnodepth\n\n];next;}
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
	my $cmd = qq[samtools view $bam $chr:].($refPos-450).qq[-].($refPos+450).qq[ | ];
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
    
    #get the coverage of the forward orientated reads
    #                   get sam                                 filter out low map Q reads                               mate is unmapped    mate not mapped in proper pair             mapped of fwd strand
    my $cmd = qq[samtools view -h $bam $chr:$start-$end | gawk -F"\\t" '(\$1~/^\@/||\$4>=$minQ)'].q[ | gawk -F"\t" '($1~/^@/||and($2,0x0008)||and($2,0x0002)==0)' | gawk -F"\t" '$1~/^@/||and($2,0x0010)==0' | samtools view -bS - > /tmp/].qq[$$.fwd.bam];
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
    
    unlink( qq[/tmp/$$.fwd.bam] ) > 0 or print qq[WARNING: Failed to delete /tmp/$$.fwd.bam files];
    
    if( $fwdpos !~ /^\d+$/ ){print qq[WARNING: Failed to estimate fwd heterozygous breakpoints for $chr:$start-$end\n];return undef;}
    
    #get the coverage of the reverse orientated reads
    #                   get sam                                 filter out low map Q reads                               mate is unmapped    mate not mapped in proper pair             mapped of rev strand
    $cmd = qq[samtools view $bam -h $chr:$start-$end | gawk -F"\\t" '(\$1~/^\@/||\$4>=$minQ)'].q[ | gawk -F"\t" '($1~/^@/||and($2,0x0008)||and($2,0x0002)==0)' | gawk -F"\t" '$1~/^@/||and($2,0x0010)' | samtools view -bS - > /tmp/].qq[$$.rev.bam];
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
    
    unlink( qq[/tmp/$$.rev.bam] ) > 0 or print qq[WARNING: Failed to delete /tmp/$$.rev.bam files];
    
    if( $revpos !~ /^\d+$/ ){print qq[WARNING: Failed to estimate rev heterozygous breakpoints for $chr:$start-$end\n];return undef;}
    
    my $positions = [ int(( $fwdpos + $revpos ) / 2) ];
    return $positions;
}

1;