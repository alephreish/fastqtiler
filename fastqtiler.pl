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

my $verion = "1.1b";

my $help = "fastqtiler.pl version $verion
Use as:
  fastqtiler.pl [parameters] [files]

  Mandatory parameters:
    -k, --kit : sequencing kit/flowcell type (see below)

  Optional parameters:
    -n, --name      : name for output files            [default: tiles]
    -o, --outdir    : output directory                 [default: .]
    -t, --threads   : max number of threads            [default: none]
    -l, --len       : read length                      [default: taken from the first read]
    -f, --frac      : fraction of the reads to analyze [default: 0.01]
    -s, --scale     : scale factor (pixel size in FC units) [default: 1000]
    -g, --goodcolor : color for the high-quality bases [default: magenta]
    -b, --badcolor  : color for the low-quality bases  [default: green]
    -i, --lanemin   : lane number reads start at       [default: 1]
    -j, --lanemax   : lane number reads end   at       [default: kit-specific]
    -h/-v, --help/--version : print this short manual

  Files:
	either gzipped, bzipped or raw fastq, but not a mix of different types

  Supported kit types (case-insensitive):
    hiseq2000/2500:high-output = hiseq2000:high-output, hiseq2500:high-output, hs2000:ho, hs2500:ho
    hiseq2000/2500:rapid-run   = hiseq2000:rapid-run, hiseq2500:rapid-run, hs2000:rr, hs2500:rr
    hiseq3000/4000             = hiseq3000, hiseq4000, hs3000, hs4000
    hiseqX                     = hsX
    miseq                      = ms
    miseq:micro                = ms:micro
    miseq:nano                 = ms:nano
    nextseq500/550:high-output = nextseq500:high-output, nextseq550:high-output, ns500:ho, ns550:ho
    nextseq500/550:mid-output  = nextseq500:mid-output, nextseq550:mid-output, ns500:mo, ns550:mo
";
print(STDERR $help) && exit if !@ARGV;

use Getopt::Long;

my $kit;
my $name      = "tiles";
my $outdir    = ".";
my $lanemin   = 1;
my $lanemax;
my $threads   = 2;
my $len;
my $frac      = 0.01;
my $goodcolor = "magenta";
my $badcolor  = "green";
my $scale     = 1000;
my $h = 0;

GetOptions (
	"k=s" => \$kit,       "kit=s"       => \$kit,
	"n=s" => \$name,      "name=s"      => \$name,
	"o=s" => \$outdir,    "outdir=i"    => \$outdir,
	"t=i" => \$threads,   "threads=i"   => \$threads,
	"l=i" => \$len,       "len=i"       => \$len,
	"s=i" => \$scale,     "scale=i"     => \$scale,
	"f=f" => \$frac,      "frac=f"      => \$frac,
	"g=s" => \$goodcolor, "goodcolor=s" => \$goodcolor,
	"b=s" => \$badcolor,  "badcolor=s"  => \$badcolor,
	"i=i" => \$lanemax,   "lanemax=i"   => \$lanemax,
	"j=i" => \$lanemin,   "lanemin=i"   => \$lanemin,
	"h" => \$h, "v" => \$h, "version" => \$h, "help" => \$h,
) || exit 0;

print(STDERR $help) && exit if $h || !@ARGV;
die("Error: kit not specified\n")               if !defined($kit);
die("Error: incorrect output name\n")           if $name !~ /^\w+$/;
die("Error: fraction must be within (1e-08; 1]\n") if $frac > 1 || $frac < 1e-08;
die("Error: scale must be a positive number\n") if $scale < 1;
die("Error: no input files specified\n")        if !@ARGV;

use POSIX ":sys_wait_h";
use Fcntl qw(:flock SEEK_END);
use File::Path qw/make_path/;

my ($bot, $nswaths, $nlanes, $ncameras, $ntiles, $max, $may);

$kit = lc($kit);

my %hs2000_ho = (
	'hiseq2000/2500:high-output' => 1,
	'hiseq2000:high-output' => 1, 'hiseq2500:high-output' => 1,
	'hs2000:ho'             => 1, 'hs2500:ho'             => 1);
my %hs2000_rr = (
	'hiseq2000/2500:rapid-run' => 1,
	'hiseq2000:rapid-run' => 1,   'hiseq2500:rapid-run' => 1,
	'hs2000:rr'           => 1,   'hs2500:rr'           => 1);
my %hs3000    = (
	'hiseq3000/4000' => 1,
	'hiseq3000'   => 1, 'hiseq4000' => 1,
	'hs3000'      => 1, 'hs4000'    => 1);
my %hsX       = (
	'hiseqx'      => 1, 'hsx'       => 1);
my %ms        = (
	'miseq'       => 1, 'ms'        => 1);
my %ms_micro  = (
	'miseq:micro' => 1, 'ms:micro'  => 1);
my %ms_nano   = (
	'miseq:nano'  => 1, 'ms:nano'   => 1);
my %ns500_ho  = (
	'nextseq500/550:high-output' => 1,
	'nextseq500:high-output' => 1, 'nextseq550:high-output' => 1,
	'ns500:ho'               => 1, 'ns550:ho'               => 1);
my %ns500_mo  = (
	'nextseq500/550:mid-output' => 1,
	'nextseq500:mid-output'  => 1, 'nextseq550:mid-output'  => 1,
	'ns500:mo'               => 1, 'ns550:mo'               => 1);
	
if (exists $hs2000_ho{$kit}) {
	$bot = 1; $nswaths = 3; $nlanes = 8; $ncameras = 1; $ntiles = 16;
	$max = 30000; $may = 110000;
}
elsif (exists $hs2000_rr{$kit}) {
	$bot = 1; $nswaths = 2; $nlanes = 8; $ncameras = 1; $ntiles = 16;
	$max = 30000; $may = 110000;
}
elsif (exists $hs3000{$kit}) {
	$bot = 1; $nswaths = 2; $nlanes = 8; $ncameras = 1; $ntiles = 28;
	$max = 30000; $may = 110000; # not verified
}
elsif (exists $hsX{$kit}) {
	$bot = 1; $nswaths = 2; $nlanes = 8; $ncameras = 1; $ntiles = 24;
	$max = 30000; $may = 110000; # not verified
}
elsif (exists $ms{$kit}) {
	$bot = 1; $nswaths = 1; $nlanes = 1; $ncameras = 1; $ntiles = 14;
	$max = 30000; $may = 30000;
}
elsif (exists $ms_micro{$kit}) {
	$bot = 1; $nswaths = 1; $nlanes = 1; $ncameras = 1; $ntiles = 4;
	$max = 30000; $may = 30000;
}
elsif (exists $ms_nano{$kit}) {
	$bot = 0; $nswaths = 1; $nlanes = 1; $ncameras = 1; $ntiles = 2;
	$max = 30000; $may = 30000;
}
elsif (exists $ns500_ho{$kit}) {
	$bot = 1; $nswaths = 3; $nlanes = 4; $ncameras = 3; $ntiles = 12;
	$max = 30000; $may = 30000;
}
elsif (exists $ns500_mo{$kit}) {
	$bot = 1; $nswaths = 1; $nlanes = 4; $ncameras = 3; $ntiles = 12;
	$max = 30000; $may = 30000;
}
else {
	die "Error: kit '$kit' unrecognized\n";
}

if (!defined($lanemax)) {
	$lanemax = $nlanes;
}

die("Error: incorrect ending lane number: $lanemax > $nlanes\n")    if $lanemax > $nlanes;
die("Error: incorrect starting lane number: $lanemin > $lanemax\n") if $lanemin > $lanemax;
die("Couldn't create the '$outdir/png_$name' folder: $!\n") if ! -d "$outdir/png_$name" && !make_path("$outdir/png_$name");
die("Couldn't create the '$outdir/txt_$name' folder: $!\n") if ! -d "$outdir/txt_$name" && !make_path("$outdir/txt_$name");

my @colors;
get_colors();

my @pipe;

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

my %pids;

my $min = 0;
my $start   = $min / $scale;
my $tilepix = $max / $scale - $start;
my $tilepiy = $may / $scale - $start;

my $width  = $tilepix * $ntiles  * $ncameras;
my $height = $tilepiy * $nswaths * ($lanemax - $lanemin + 1);
my $cam = $ncameras > 1;
my $nf = 3 + $cam;

my @fq;
my @parts;
my $lane;
my $tile = -1;
my $eof = 0;

$fq[$_] = <DATA> foreach (0..3);
chomp ($fq[1], $fq[2]);
$len = length($fq[1]) if !defined($len);
die "Error: not a valid fastq file\n" if $fq[2] ne "+" || !$fq[0] || substr($fq[0], 0, 1) ne "@";
@parts = split(":", $fq[0]);
die "Error: unrecognized fastq header\n" if $#parts < 6;

my $head = sprintf("# ImageMagick pixel enumeration: %d,%d,255,srgb\n", $width, $height);

for (my $i = 0; $i < $len; $i++) {
	open(FP, ">", fname($i, 1, "txt")) || die "Error: failed to open output file: $!\n";
	print FP $head;
	close(FP);
	if ($bot) {
		open(FP, ">", fname($i, 2, "txt")) || die "Error: failed to open output file: $!\n";
		print FP $head;
		close(FP);
	}
}

while (!eof(DATA)) {
	$fq[$_] = <DATA> foreach (0..3);
	next if rand() > $frac;
	@parts = split(":", $fq[0], 6);
	$lane = $parts[3];
	next if $lane < $lanemin || $lane > $lanemax;
	if ($tile != $parts[4]) {
		close KID if $tile > 0;
		$tile = $parts[4];
		die "Error: incorrect fastq header: $fq[0]\n" if length($tile) < 4 || $tile < 0;
		check_kids();
		printf "\rStarting with tile $tile";
		my $pid = pipe_to_fork();
		do_job() if !$pid;
		$pids{$pid} = 1;
		printf KID "%d:%d\n", $lane, $tile;
	}
	chomp $parts[5];
	printf KID "%s:%s", $parts[5], $fq[3];
}
close(KID);
close(DATA);

foreach my $pid (keys %pids) {
	sleep 0.1 while (waitpid($pid, WNOHANG) >= 0);
}
%pids = ();
wait;

print "\r                                        ";

for (my $i = 0; $i < $len; $i++) {
#	check_kids();
#	my $pid = fork();
	printf "\rStarting conversion of pos %d", $i + 1;
#	if ($pid) {
#		$pids{$pid} = 1;
#		next;
#	}
	system("convert", fname($i, 1, "txt"), fname($i, 1, "png"));
	system("convert", fname($i, 2, "txt"), fname($i, 2, "png")) if $bot;
#	exit;
}

wait;

open HTML, ">", "$outdir/$name.html";
print HTML html_top();
print HTML table_head();
print HTML table_body();
print HTML html_bottom();
close HTML;

print "\rDone                                    \n";

sub check_kids {
	while (1) {
		foreach my $pid (keys %pids) {
			if (waitpid($pid, WNOHANG) != 0) {
				delete $pids{$pid};
				last;
			}
		}
		last if scalar keys %pids < $threads;
		sleep 1;
	}
}

sub fname {
	my $cycle   = shift;
	my $surface = shift;
	my $dir     = shift;
	return sprintf("%s/%s_%s/%d.%03d.%s", $outdir, $dir, $name, $surface, $cycle + 1, $dir);
}

sub do_job {
	my %quals;
	my @quala;
	my ($px, $py);

	my $header = <STDIN>;
	chomp $header;
	my @p = split(":", $header, 2);
	my $mylane = $p[0] - $lanemin;
	my $mytile = $p[1];
	my @f = split("", $mytile, $nf);

	my $surface = $f[0]  - 1;
	my $swath   = $f[1]  - 1;
	my $tilenum = $f[-1] - 1;
	my $camera = 0;

	if ($cam) {
		$camera = $f[2] - 1 ;
		$camera = 5 - $camera if $camera > 2;
	}

	my $offsetx = $tilepix * ($camera * $ntiles  + $tilenum) - $start;
	my $offsety = $tilepiy * ($mylane * $nswaths + $swath  ) - $start;

	while (<STDIN>) {
		chomp;
		@p = split(/[: ]/);
		$px = sprintf('%d', $offsetx + $p[0] / $scale);
		$py = sprintf('%d', $offsety + $p[1] / $scale);
		next if $px > $width || $py > $height || $px < 1 || $py < 1;
		push(@{$quals{$px}{$py}}, $p[-1]);
	}
	close(STDIN);
	my @sum;
	my @num;
	my $qual;
	foreach $px (keys %quals) {
		foreach $py (keys %{$quals{$px}}) {
			foreach $qual (@{$quals{$px}{$py}}) {
				my $i = 0;
				foreach my $q (split("", $qual)) {
					last if $i == $len;
					$sum[$i]{$px}{$py} += ord($q);
					$num[$i]{$px}{$py}++;
					$i++;
				}
			}
		}
	}
	my %mysum;
	my %mynum;
	my $avg;

	for (my $i = 0; $i <= $#sum; $i++) {
		%mysum = %{$sum[$i]};
		%mynum = %{$num[$i]};
		open(FP, ">>", fname($i, $surface + 1, "txt")) || die "Error: failed to open output file: $!\n";
		flock(FP, LOCK_EX)    || die "Error: cannot lock: $!\n";
		seek(FP, 0, SEEK_END) || die "Error: cannot seek: $!\n";
		foreach $px (keys %mysum) {
			foreach $py (keys %{$mysum{$px}}) {
				$avg = sprintf('%d', ($mysum{$px}{$py} / $mynum{$px}{$py} - 33) / 41 * 255);
				die "Error: unexpected average quality: $avg\n" if $avg > 255 || $avg < 0;
				printf FP "%d,%d:(%s)\n", $px, $py, $colors[$avg];
			}
		}
		flock(FP, LOCK_UN) || die "Cannot unlock: $!\n";
		close(FP);
	}
	exit;
}

sub pipe_to_fork {
	pipe(my $child, KID) || die;
	my $pid = fork();
	die "fork() failed: $!\n" unless defined $pid;
	if ($pid) {
		close $child;
	}
	else {
		close KID;
		open(STDIN, "<&=" . fileno($child)) || die;
	}
	return $pid;
}

sub get_colors {
	open(COL, "-|", "convert -size 1x256 'gradient:$badcolor-$goodcolor' -depth 8 -colorspace RGB txt:");
	my @parts;
	my $q;
	my $rgb;
	my @p;
	while (<COL>) {
		@p = split(/:\s*\(|\)/);
		next if substr($p[0], 0, 1) eq "#";
		$q   = substr($p[0], 2);
		$rgb = $p[1];
		$rgb =~ s/\s+//g;
		$colors[$q] = $rgb;
	}
	close(COL);
	die "Error: convert failed. Check the colors and that the ImageMagick's convert tool is installed and available on the \$PATH\n" if $? > 0;
}

sub table_body {
	my $tbody = "<tbody>";
	my $td = '<td><img src="png_%s/%d.%03d.png" width="%dpx" /></td>';
	my $w = $width < 1000 ? $width : 1000;
	for (my $i = 1; $i <= $len; $i++) {
		my $img  = sprintf($td, $name, 1, $i, $w);
		$img    .= sprintf($td, $name, 2, $i, $w) if $bot;
		$tbody  .= sprintf("<tr><td>%d</td>%s</tr>\n", $i, $img);
	}
	$tbody .= "</tbody>";
	return $tbody;
}

sub table_head {
	my $th1  = '<th>Top tiles</th>';
	$th1 .= '<th>Bottom tiles</th>' if $bot;
	my $tiles = "";
	$tiles .= sprintf("%02d ", $_) foreach 1..$ntiles;
	$tiles = substr($tiles x $ncameras, 0, -1);
	my $th2 = "";
	my $th = '<th class="tilenum">&nbsp;%s&nbsp;</th>';
	$th2 .= sprintf($th, $tiles);
	$th2 .= sprintf($th, $tiles) if $bot;
	return sprintf('<thead>
<tr>
<th rowspan="2">Cycle</th>%s</tr>
<tr>
%s</tr>
</thead>', $th1, $th2);
}

sub html_top {
	return <<'HEAD';
<!DOCTYPE html>
<html>
<head>

<title>Tiles</title>
<meta http-equiv="content-type" content="text/html; charset=UTF-8">
<script type='text/javascript' src='http://code.jquery.com/jquery-1.9.1.js'></script>

<style type='text/css'>
table {
	border: 0;
	padding: 0;
	margin: 0 0 20px 0;
	border-collapse: collapse;
}
th {
	padding: 5px;
	text-align: center;
	font-weight:bold;
	line-height: 1em;
	color: #FFF;
	background-color: #555;
}
tbody td {
	padding: 5px;
	line-height: 18px;
	border-top: 1px solid #E0E0E0;
}
tbody tr:hover {
	background-color: #EEEEEE;
}
.tilenum {
	text-align:justify;
	padding:5px 5px 0px 5px;
	font-size:10px;
}
.tilenum:after {
  content: "";
  display: inline-block;
  width: 100%;
}
</style>

</head><body><table>
HEAD
}

sub html_bottom {
	return <<'BOT';
	</table>

<script type='text/javascript'>//<![CDATA[

$(function () {
	$("table").stickyTableHeaders();
});

/*! Copyright (c) 2011 by Jonas Mosbech - https://github.com/jmosbech/StickyTableHeaders
	MIT license info: https://github.com/jmosbech/StickyTableHeaders/blob/master/license.txt */

;
(function ($, window, undefined) {
	'use strict';

	var name = 'stickyTableHeaders',
	    id = 0,
	    defaults = {
	        fixedOffset: 0,
	        leftOffset: 0,
	        marginTop: 0,
	        scrollableArea: window
	    };

	function Plugin(el, options) {
	    // To avoid scope issues, use 'base' instead of 'this'
	    // to reference this class from internal events and functions.
	    var base = this;

	    // Access to jQuery and DOM versions of element
	    base.$el = $(el);
	    base.el = el;
	    base.id = id++;
	    base.$window = $(window);
	    base.$document = $(document);

	    // Listen for destroyed, call teardown
	    base.$el.bind('destroyed',
	    $.proxy(base.teardown, base));

	    // Cache DOM refs for performance reasons
	    base.$clonedHeader = null;
	    base.$originalHeader = null;

	    // Keep track of state
	    base.isSticky = false;
	    base.hasBeenSticky = false;
	    base.leftOffset = null;
	    base.topOffset = null;

	    base.init = function () {
	        base.$el.each(function () {
	            var $this = $(this);

	            // remove padding on <table> to fix issue #7
	            $this.css('padding', 0);

	            base.$originalHeader = $('thead:first', this);
	            base.$clonedHeader = base.$originalHeader.clone();
	            $this.trigger('clonedHeader.' + name, [base.$clonedHeader]);

	            base.$clonedHeader.addClass('tableFloatingHeader');
	            base.$clonedHeader.css('display', 'none');

	            base.$originalHeader.addClass('tableFloatingHeaderOriginal');

	            base.$originalHeader.after(base.$clonedHeader);

	            base.$printStyle = $('<style type="text/css" media="print">' +
      '.tableFloatingHeader{display:none !important;}' +
      '.tableFloatingHeaderOriginal{position:static !important;}' +
      '</style>');
	            $('head').append(base.$printStyle);
	        });

	        base.setOptions(options);
	        base.updateWidth();
	        base.toggleHeaders();
	        base.bind();
	    };

	    base.destroy = function () {
	        base.$el.unbind('destroyed', base.teardown);
	        base.teardown();
	    };

	    base.teardown = function () {
	        if (base.isSticky) {
	            base.$originalHeader.css('position', 'static');
	        }
	        $.removeData(base.el, 'plugin_' + name);
	        base.unbind();

	        base.$clonedHeader.remove();
	        base.$originalHeader.removeClass('tableFloatingHeaderOriginal');
	        base.$originalHeader.css('visibility', 'visible');
	        base.$printStyle.remove();

	        base.el = null;
	        base.$el = null;
	    };

	    base.bind = function () {
	        base.$scrollableArea.on('scroll.' + name, base.toggleHeaders);
	        if (!base.isWindowScrolling) {
	            base.$window.on('scroll.' + name + base.id, base.setPositionValues);
	            base.$window.on('resize.' + name + base.id, base.toggleHeaders);
	        }
	        base.$scrollableArea.on('resize.' + name, base.toggleHeaders);
	        base.$scrollableArea.on('resize.' + name, base.updateWidth);
	    };

	    base.unbind = function () {
	        // unbind window events by specifying handle so we don't remove too much
	        base.$scrollableArea.off('.' + name, base.toggleHeaders);
	        if (!base.isWindowScrolling) {
	            base.$window.off('.' + name + base.id, base.setPositionValues);
	            base.$window.off('.' + name + base.id, base.toggleHeaders);
	        }
	        base.$scrollableArea.off('.' + name, base.updateWidth);
	    };

	    base.toggleHeaders = function () {
	        if (base.$el) {
	            base.$el.each(function () {
      var $this = $(this),
          newLeft,
          newTopOffset = base.isWindowScrolling ? (
          isNaN(base.options.fixedOffset) ? base.options.fixedOffset.outerHeight() : base.options.fixedOffset) : base.$scrollableArea.offset().top + (!isNaN(base.options.fixedOffset) ? base.options.fixedOffset : 0),
          offset = $this.offset(),

          scrollTop = base.$scrollableArea.scrollTop() + newTopOffset,
          scrollLeft = base.$scrollableArea.scrollLeft(),

          scrolledPastTop = base.isWindowScrolling ? scrollTop > offset.top : newTopOffset > offset.top,
          notScrolledPastBottom = (base.isWindowScrolling ? scrollTop : 0) < (offset.top + $this.height() - base.$clonedHeader.height() - (base.isWindowScrolling ? 0 : newTopOffset));

      if (scrolledPastTop && notScrolledPastBottom) {
          newLeft = offset.left - scrollLeft + base.options.leftOffset;
          base.$originalHeader.css({
              'position': 'fixed',
                  'margin-top': base.options.marginTop,
                  'left': newLeft,
                  'z-index': 3 // #18: opacity bug
          });
          base.leftOffset = newLeft;
          base.topOffset = newTopOffset;
          base.$clonedHeader.css('display', '');
          if (!base.isSticky) {
              base.isSticky = true;
              // make sure the width is correct: the user might have resized the browser while in static mode
              base.updateWidth();
          }
          base.setPositionValues();
      } else if (base.isSticky) {
          base.$originalHeader.css('position', 'static');
          base.$clonedHeader.css('display', 'none');
          base.isSticky = false;
          base.resetWidth($('td,th', base.$clonedHeader), $('td,th', base.$originalHeader));
      }
	            });
	        }
	    };

	    base.setPositionValues = function () {
	        var winScrollTop = base.$window.scrollTop(),
	            winScrollLeft = base.$window.scrollLeft();
	        if (!base.isSticky || winScrollTop < 0 || winScrollTop + base.$window.height() > base.$document.height() || winScrollLeft < 0 || winScrollLeft + base.$window.width() > base.$document.width()) {
	            return;
	        }
	        base.$originalHeader.css({
	            'top': base.topOffset - (base.isWindowScrolling ? 0 : winScrollTop),
      'left': base.leftOffset - (base.isWindowScrolling ? 0 : winScrollLeft)
	        });
	    };

	    base.updateWidth = function () {
	        if (!base.isSticky) {
	            return;
	        }
	        // Copy cell widths from clone
	        if (!base.$originalHeaderCells) {
	            base.$originalHeaderCells = $('th,td', base.$originalHeader);
	        }
	        if (!base.$clonedHeaderCells) {
	            base.$clonedHeaderCells = $('th,td', base.$clonedHeader);
	        }
	        var cellWidths = base.getWidth(base.$clonedHeaderCells);
	        base.setWidth(cellWidths, base.$clonedHeaderCells, base.$originalHeaderCells);

	        // Copy row width from whole table
	        base.$originalHeader.css('width', base.$clonedHeader.width());
	    };

	    base.getWidth = function ($clonedHeaders) {
	        var widths = [];
	        $clonedHeaders.each(function (index) {
	            var width, $this = $(this);

	            if ($this.css('box-sizing') === 'border-box') {
      width = $this[0].getBoundingClientRect().width; // #39: border-box bug
	            } else {
      var $origTh = $('th', base.$originalHeader);
      if ($origTh.css('border-collapse') === 'collapse') {
          if (window.getComputedStyle) {
              width = parseFloat(window.getComputedStyle(this, null).width);
          } else {
              // ie8 only
              var leftPadding = parseFloat($this.css('padding-left'));
              var rightPadding = parseFloat($this.css('padding-right'));
              // Needs more investigation - this is assuming constant border around this cell and it's neighbours.
              var border = parseFloat($this.css('border-width'));
              width = $this.outerWidth() - leftPadding - rightPadding - border;
          }
      } else {
          width = $this.width();
      }
	            }

	            widths[index] = width;
	        });
	        return widths;
	    };

	    base.setWidth = function (widths, $clonedHeaders, $origHeaders) {
	        $clonedHeaders.each(function (index) {
	            var width = widths[index];
	            $origHeaders.eq(index).css({
      'min-width': width,
          'max-width': width
	            });
	        });
	    };

	    base.resetWidth = function ($clonedHeaders, $origHeaders) {
	        $clonedHeaders.each(function (index) {
	            var $this = $(this);
	            $origHeaders.eq(index).css({
      'min-width': $this.css('min-width'),
          'max-width': $this.css('max-width')
	            });
	        });
	    };

	    base.setOptions = function (options) {
	        base.options = $.extend({}, defaults, options);
	        base.$scrollableArea = $(base.options.scrollableArea);
	        base.isWindowScrolling = base.$scrollableArea[0] === window;
	    };

	    base.updateOptions = function (options) {
	        base.setOptions(options);
	        // scrollableArea might have changed
	        base.unbind();
	        base.bind();
	        base.updateWidth();
	        base.toggleHeaders();
	    };

	    // Run initializer
	    base.init();
	}

	// A plugin wrapper around the constructor,
	// preventing against multiple instantiations
	$.fn[name] = function (options) {
	    return this.each(function () {
	        var instance = $.data(this, 'plugin_' + name);
	        if (instance) {
	            if (typeof options === 'string') {
      instance[options].apply(instance);
	            } else {
      instance.updateOptions(options);
	            }
	        } else if (options !== 'destroy') {
	            $.data(this, 'plugin_' + name, new Plugin(this, options));
	        }
	    });
	};

})(jQuery, window);
//]]> 

</script>


</body>
</html>
BOT
}
