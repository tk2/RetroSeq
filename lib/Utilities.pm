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
        my $i = 0;
        foreach my $call (  @calls )
        {
            next unless $call;
            my @s1 = split( /\t/, $call );
            croak qq[Badly formatted call: $call] unless @s1==5;
            if( $s[ 0 ] eq $s1[ 0 ] && ( 
                ( $s1[ 1 ] > $s[ 1 ] && $s1[ 2 ] < $s[ 2 ] ) #totally enclosed in region
                ||
                ( abs( $s1[ 1 ] - $s[ 1 ] ) < $FILTER_WINDOW )
                ||
                ( abs( $s1[ 1 ] - $s[ 2 ] ) < $FILTER_WINDOW )
                ||
                ( abs( $s1[ 2 ] - $s[ 1 ] ) < $FILTER_WINDOW )
                ||
                ( abs( $s1[ 2 ] - $s[ 2 ] ) < $FILTER_WINDOW )
                )
              )
            {
                print qq[Excluding region: $calls[ $i ]\n];
                undef( $calls[ $i ] );
                $filtered ++;
            }
            $i ++;
        }
    }
    close( $ffh );

    open( my $ofh, qq[>$outputBED] ) or die $!;
    foreach my $call ( @calls ){if(defined( $call ) ){ print $ofh qq[$call\n]; } }
    close( $ofh );
    
    return $filtered;
}

1;
