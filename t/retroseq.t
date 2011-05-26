#!/usr/bin/perl -w

use strict;
use File::Basename;
use Cwd;
use Test::More qw(no_plan);

use lib dirname(__FILE__).'/../lib/';
use Utilities;

if( -f qq[data/filtered.bed]){unlink(qq[data/filtered.bed]);}

is(Utilities::filterOutRegions(qq[data/input.bed], qq[data/reference_set.bed], qq[data/filtered.bed]), 3, "region filter checks");


