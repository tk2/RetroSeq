package Utilities;

use strict;
use warnings;
use Carp;

my $FILTER_WINDOW = 50;

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
=pod
sub testRegion
{
    my $chr = shift;
    my $refPos = shift;
    my $bam = shift;
    my $minMapQ
    
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
	        #print qq[$refPos\t$depth\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
	        if( $lhsFwdBlue >= 5 && $rhsRevBlue >= 5 && $lhsFwd > 10 && $rhsRev > 10 && ( $lhsRevBlue == 0 || $lhsFwdBlue / $lhsRevBlue > 2 ) && ( $rhsFwdBlue == 0 || $rhsRevBlue / $rhsFwdBlue > 2 ) && $dist < 120 )
	        {
	            my $ratio = ( $lhsRev + $rhsFwd ) / ( $lhsFwd + $rhsRev ); #objective function is to minimise this value (i.e. min depth, meets the criteria, and balances the 3' vs. 5' ratio best)
	            if( $ratio < $minRatio ){$minRatioCall = qq[$chr\t$refPos\t].($refPos+1).qq[\t$originalCallA[ 3 ]\t$originalCallA[ 4 ]\n];$minRatio = $ratio;}
	            $found = 1;
	            print $dfh qq[PASS: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$depth\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];
	        }
	        else{print $dfh qq[FILTER: $originalCallA[ 0 ]\t$refPos\t$refPos\t$originalCallA[3]_filter\t$depth\t$lhsFwdBlue\t$lhsRevBlue\t$lhsFwdGreen\t$lhsRevGreen\t$rhsFwdBlue\t$rhsRevBlue\t$rhsFwdGreen\t$rhsRevGreen\t$lastBluePos\t$firstBluePos\t$dist\n];}	        
	        
	        if( $found )
	        {
	            return $minRatioCall;
	        }
	        else{return undef;}
}
=cut

1;
