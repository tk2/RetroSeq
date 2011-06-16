#!/usr/bin/perl -w

use strict;
use File::Basename;
use Cwd;
use Test::More qw(no_plan);

use lib dirname(__FILE__).'/../lib/';
use Utilities;

if( -f qq[data/filtered.bed]){unlink(qq[data/filtered.bed]);}

is(Utilities::filterOutRegions(qq[data/input.bed], qq[data/reference_set.bed], qq[data/filtered.bed]), 9, "region filter checks");

my $t = Utilities::getCandidateBreakPointsDir( 2, 74124481, 74126081, '/lustre/scratch102/user/tk2/RetroSeq/Human/NA18506.raw.bam', 30, 10 );
isa_ok( $t, 'ARRAY' );
my @pos = @{$t};
foreach(@pos){print qq[$_\n];}
