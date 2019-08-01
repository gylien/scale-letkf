#!/usr/bin/perl -w
#
# Ex.
#	./ltestgen.pl ../testdata/JITDTLOG.20170510 ../testdata/data.tar /tmp/ 2
#   In this case,
#     something like the following commands are executed each interval time:
#	cd /tmp/;
#	tar --to-stdout -x -f /home/ishikawa/jit-dt/src/../testdata/data.tar \
#	       kobe_20170510170000_A08_pawr_vr.dat.gz \
#	    | gunzip > kobe_20170510170000_A08_pawr_vr.dat
#
use Cwd;
use Time::HiRes qw(setitimer ITIMER_REAL);

if (@ARGV < 3) {
    die("ltestgen.pl <log file> <tar file> <out directory> " .
	            "[inteval] [count] \n");
}
open(INPUT, $ARGV[0]) or die("Can't open file: $ARGV[0]");
$cwd = Cwd::getcwd();
$targzfile = $ARGV[1];
if ($targzfile !~ "^/" && $targzfile !~ "^~") {
    $targzfile = $cwd . "/" . $targzfile;
}
$dir = $ARGV[2];
if (@ARGV >= 4) {
    $watch_interval = int($ARGV[3]);
} else {
    $watch_interval = 30;
}
$first_sleep = 1;
$count = -1;
if (@ARGV == 5) {
    $count = int($ARGV[4]);
}
$retval = 0;
printf("Interval: %d sec, Count: %d\n", $watch_interval, $count);

$SIG{ALRM} = sub {
#    printf("ALRAM: %s\n", $cmd);
    if (eof(INPUT)) {
	printf("Finish.\n");
	exit 0;
    }
#    $nfiles = 2;
#    Modified to 3 by T.Honda (07/31/2019)
    $nfiles = 3;
    while ($nfiles-- > 0) {
	do {
	    $line = <INPUT>;
	} while ($line !~ /^[\d]/);
	($_, $fname,$_) = split(',', $line);
	$fname =~ /\s+(.+).dat/;
	$fgzname = $fname . ".gz";
	$cmd = "cd " . $dir . "; tar --to-stdout -x -f " . $targzfile . $fgzname
	    . " | gunzip >" . $fname;
	$retval = system $cmd;
	$time = localtime(time);
	if ($retval == 0) {
	    printf("[%s] %s\n", $time, $cmd);
	} else {
	    printf("[%s] exit code(%d) %s\n", $time, $retval, $cmd);
	    exit 0;
	}
    }
    if ($count > 0) {
	$count--;
	if ($count == 0) {
	    printf("Finish.\n");
	    exit 0;
	}
    }
};

setitimer(ITIMER_REAL, $first_sleep, $watch_interval);

while (<>) {
    sleep(1000);
}

#
#while ($line = <INPUT>) {
#    if ($line !~ /^[\d]/) { next; }
#    ($_, $fname,$_) = split(',', $line);
#    $fname =~ /\s+(.+).dat/;
#    $cmd = "tar gzf " . $targzfile . $fname;
#    printf("%s\n", $cmd);
#    $wpath = $dir . $1 . ".dat";
#    printf("%s\n", $wpath);
#    open(OUT, "> $wpath");
#    printf OUT "$wpath\n";
#    close(OUT);
#    sleep(1000);
#}
