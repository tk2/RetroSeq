#!/usr/bin/perl -w

use strict;
use File::Basename;
use Cwd;
use Test::More qw(no_plan);

use lib dirname(__FILE__).'/../lib/';
use Utilities;

Utilities::checkBinary( q[samtools], qq[0.1.16] );
Utilities::checkBinary( q[bcftools] );
Utilities::checkBinary( q[exonerate], qq[2.2.0] );

if( -f qq[data/filtered.bed]){unlink(qq[data/filtered.bed]);}

is(Utilities::filterOutRegions(qq[data/input.bed], qq[data/reference_set.bed], qq[data/filtered.bed]), 9, "region filter checks");

my @bams = (qq[data/bam_header.bam]);
is(Utilities::getBAMSampleName(\@bams),'QTL210645', 'Sample identifier test' );

my $t = Utilities::getCandidateBreakPointsDir( 2, 74124481, 74126081, '/lustre/scratch102/user/tk2/RetroSeq/Human/striped/NA18506.raw.bam', 30 );
isa_ok( $t, 'ARRAY' );
my @pos = @{$t};
foreach(@pos){print qq[$_\n];}

$t = Utilities::getCandidateBreakPointsDir( 1, 4509962, 4510033, '/lustre/scratch102/user/tk2/RetroSeq/Human/striped/NA18506.raw.bam', 30 );
is( $t, undef, "Undefined breakpoint returned as expected" );
#my @pos = @{$t};
#foreach(@pos){print qq[$_\n];}

$t = Utilities::getCandidateBreakPointsDir( 1, 246219769, 246220104, '/lustre/scratch102/user/tk2/RetroSeq/Human/striped/NA18506.raw.bam', 30 );
isa_ok( $t, 'ARRAY' );
@pos = @{$t};
foreach(@pos){print qq[$_\n];}

@bams = ('/lustre/scratch102/user/tk2/RetroSeq/Human/striped/NA18506.raw.bam','/lustre/scratch102/user/tk2/RetroSeq/Human/striped/NA18507.raw.bam');

open( my $dfh, qq[>STDOUT] ) or die $!;
$t = Utilities::testBreakPoint(1,246219988,\@bams,20,qq[1	246219744	246220162	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(10,56902335,\@bams,20,qq[10	56902376	56902359	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(15,94819895,\@bams,20,qq[15	94818895	94820895	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(16,82912110,\@bams,20,qq[16	82911110	82913110	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(2,74125457,\@bams,20,qq[2	74124457	74126457	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(11,428063,\@bams,20,qq[11	428029	428070	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(20,2252212,\@bams,20,qq[20	2251212	2253212	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(20,2252177,\@bams,20,qq[20	2251212	2253212	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(1,7486242,\@bams,20,qq[11	428029	428070	Alu	10],$dfh);
isa_ok( $t, 'ARRAY' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::getCandidateBreakPointsDirVote(20,53560000,53561000,\@bams,20);
isa_ok( $t, 'ARRAY', "break point vote1" );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::getCandidateBreakPointsDirVote(11,25671300,25672000,\@bams,20);
isa_ok( $t, 'ARRAY', "break point vote2" );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::getCandidateBreakPointsDirVote(11,427600,428500,\@bams,20);
isa_ok( $t, 'ARRAY', "break point vote3" );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::getCandidateBreakPointsDirVote(20,34458000,34459000,\@bams,20);
isa_ok( $t, 'ARRAY', "break point vote3" );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(20,34458559,\@bams,20,qq[20	428029	428070	Alu	10],$dfh);
isa_ok( $t, 'ARRAY', 'break point test' );
foreach( @{$t} ){print qq[$_\n];}

$t = Utilities::testBreakPoint(6,29892876,\@bams,20,qq[20	428029	428070	Alu	10],$dfh);
isa_ok( $t, 'ARRAY', 'break point test1' );
foreach( @{$t} ){print qq[$_\n];}

close( $dfh );
