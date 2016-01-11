#!/usr/bin/perl

### Written by Andrey Rozenberg (jaera at yandex.com), Ruhr-Universität Bochum

### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.;

### You should have received a copy of the GNU General Public License
### along with this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use POSIX ":sys_wait_h";

my $f = 0.001;
my @fq;
my @pipe;

if (!@ARGV) {
	print STDERR "Error: no input specified\ntry fastqtester.pl -h\n";
	exit;
}

if ($ARGV[0] eq '-h' || $ARGV[0] eq '-help' || $ARGV[0] eq '--help' || $ARGV[0] eq '-v' || $ARGV[0] eq '-version' || $ARGV[0] eq '--version') {
	print STDERR "Use this program as follows: fastqtester.pl <input files>
	the input files must belong to the same run and be in fastq format, optionally gzip'ed or bzip2'ed\n";
	exit;
}

if    ($ARGV[0] =~ /.gz$/) {
	@pipe = ("gzip", "-cd");
}
elsif ($ARGV[0] =~ /.bz2?$/) {
	@pipe = ("bzip2", "-cd");
}
elsif ($#ARGV > 0) {
	@pipe = ("cat");
}

if (@pipe) {
	open DATA, "-|", @pipe, @ARGV or die "Error: failed to open pipe: $!\n";
}
else {
	open DATA, "<", $ARGV[0] or die "Error: failed to open file: $!\n";;
}

$fq[$_] = <DATA> foreach (0..3);
my @parts = split(":", $fq[0]);
my $cam   = length($parts[4]) == 5;

my $instr  = substr($parts[0], 1);

printf "Instrument:      %s\n", $instr;
printf "Run:             %s\n", $parts[1];
printf "Flowcell ID:     %s\n", $parts[2];

my $lanemax  = 0;
my $tilemax  = 0;
my $swathmax = 0;
my @fields;
my $bot       = 0;
my $lenmax    = 0;
my $nreads    = 0;
my $cameramax = 1;
my $ntiles    = 0;
my $nparts = 3 + $cam;

my ($surface, $lane, $swath, $camera, $tilenum, $slen);

my $tile = 0;

while (!eof(DATA)) {
	$fq[$_] = <DATA> foreach (0..3);
	$nreads++;
	if (rand() < $f) {
		chomp $fq[1];
		$slen   = length($fq[1]);
		$lenmax = $slen if $lenmax < $slen;
	}
	@parts = split(":", $fq[0], 6);
	next if $tile == $parts[4];
	$ntiles++;

	$lane    = $parts[3];
	$tile    = $parts[4];
	@fields  = split("", $tile, $nparts);

	$surface = $fields[0];
	$swath   = $fields[1];
	$tilenum = $fields[-1];

	$lanemax   = $lane    if $lanemax  < $lane;
	$tilemax   = $tilenum if $tilemax  < $tilenum;
	$swathmax  = $swath   if $swathmax < $swath;
	$bot       = 1        if $surface == 2;
	if ($cam) {
		$camera = $fields[2];
		$cameramax = $camera if $cameramax < $camera;
	}
}

my $kit       = "unrecognized";
my $ntilesmax = 10000;

if    ($swathmax == 3 && $cameramax == 1 && $tilemax == 16 && $lanemax < 9) {
	$kit = "hiseq2000/2500:high-output";
	$ntilesmax = 768;
}
elsif ($swathmax == 2 && $cameramax == 1 && $tilemax == 16 && $lanemax < 9) {
	$kit = "hiseq2000/2500:rapid-run";
	$ntilesmax = 512;
}
elsif ($swathmax == 2 && $cameramax == 1 && $tilemax == 28 && $lanemax < 9) {
	$kit = "hiseq3000/4000";
	$ntilesmax = 896;
}
elsif ($swathmax == 2 && $cameramax == 1 && $tilemax == 24 && $lanemax < 9) {
	$kit = "hiseqX";
	$ntilesmax = 768;
}
elsif ($swathmax == 1 && $cameramax == 1 && $tilemax == 14 && $lanemax < 2) {
	$kit = "miseq";
	$ntilesmax = 28;
}
elsif ($swathmax == 1 && $cameramax == 1 && $tilemax == 4  && $lanemax < 2) {
	$kit = "miseq:micro";
	$ntilesmax = 8;
}
elsif ($swathmax == 1 && $cameramax == 1 && $tilemax == 2  && $lanemax < 2 && !$bot) {
	$kit = "miseq:nano";
	$ntilesmax = 2;
}
elsif ($swathmax == 3 && $cameramax % 3 == 0 && $tilemax == 12 && $lanemax < 5) {
	$kit = "nextseq:high-output";
	$ntilesmax = 864;
}
elsif ($swathmax == 1 && $cameramax % 3 == 0 && $tilemax == 12 && $lanemax < 5) {
	$kit = "nextseq:mid-output";
	$ntilesmax = 288;
}

printf "Flowcell:        %s\n", $kit;
printf "Lanes:           %d\n", $lanemax;
printf "Cameras:         %d\n", $cameramax;
printf "Swaths:          %d\n", $swathmax;
printf "Tiles per block: %d\n", $tilemax;
printf "Max read length: %d\n", $lenmax;
printf "Bottom surface:  %s\n", $bot ? "yes" : "no";
printf "Number of reads: %d\n", $nreads;
printf "Number of tiles: %d\n", $ntiles;
printf "The sequences are not sorted!\n" if $ntiles > $ntilesmax;
