#!/usr/bin/perl

# Run a few standard cases and compare the results to the reference
# output stored in maint/selftest_ref.  For each test, the numbers in
# the new output are compared to the reference numbers, and the number
# of significant digits to which they agree is reported.
#
# Usage:  maint/selftest.pl [-s] [-m MINDIGITS] [TEST ...]
#
#   -s      Store the current output as the new reference.
#   -m      Minimum number of significant digits needed to pass (default 5).
#   TEST    Run only the named tests.
#
# Two numbers a (reference) and b (new) are compared using
#
#     |a-b| / max(|a|, |b|, scale)
#
# where scale is the natural scale of quantities that can be close to zero:
#  - in data lines  "lam kabs ksca g",          scale = 1 for g
#  - in data lines  "F11 F12 F22 F33 F34 F44",  scale = F11 for F12...F44
# For all other numbers scale = 0, so they are compared relative to their
# own value.  Without the scale, the cancellation noise in, for example,
# g = 1e-6 in the Rayleigh limit, or F34 near 180 degrees, would dominate
# the result.  The agreement in digits is -log10 of the largest relative
# difference, reported as "printed precision" when it exceeds the number
# of digits that optool prints.  Any change in text, or in the number of
# lines or numbers, fails a test.

use Getopt::Std;
use File::Basename;
use File::Path qw(make_path remove_tree);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip;
getopts('sm:');
$mindigits = defined $opt_m ? $opt_m : 5;
$dir = 'selftest_optool';
$refdir = dirname($0) . '/selftest_ref';
%only = map { $_ => 1 } @ARGV;

my $tests = [
  {name => 'noscat',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'full-diana-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'sparse-file',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -na 20 -sparse 2.2 -sp 300 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'lognormal-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -a 0.001 4.9 0.1:1.0 30 -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'diana-plus-ice-mantle',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -na 20 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'div-in-10-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -d 3 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa_*.dat"},

  {name => 'large-grain',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.696 -c c-z 0.104 -m h2o-w 0.2 -p 0.25 -a 1000 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 10 -na 20 -s 720 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'chop-peak',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -chop 4 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'cde-rayleigh',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -a 0.01 0.1 3.5 10 -l 10 60 30 -s -cde -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'mmf-scat',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool pyr 1 -q -a 10 30 -mmf 0.3 .005 -l 30 3000 20 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'mmf-opac',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr 1 -q -a 10 30 3.5 2 -mmf 0.3 2.2 -nlam 3 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-astrosil',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool astrosil 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-gra',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-gra 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-nano',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-nano 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-org',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-org 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-p',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-p 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-z',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-z 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch3oh-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch3oh-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch3oh-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch3oh-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch4-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch4-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch4-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch4-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co2-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co2-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co2-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co2-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co2-w',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co2-w 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-cor-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool cor-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-fe-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool fe-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-fes',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool fes 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-h2o-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool h2o-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-h2o-w',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool h2o-w 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-nh3-m',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool nh3-m 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-c-mg00',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-c-mg00 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-c-mg100',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-c-mg100 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-c-mg95',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-c-mg95 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-mg40',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-mg40 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-mg50',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-mg50 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-c-mg96',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-c-mg96 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg100',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg100 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg40',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg40 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg50',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg50 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg60',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg60 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg70',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg70 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg80',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg80 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg95',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg95 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-sic',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool sic 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-sio2',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool sio2 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => 'fb1a85451e94e2d4d8a5d05317ebde9b6fc9e48e'}
  ];


die "./optool file not present" if not -e "./optool";
die "./optool is not executable" if not -x "./optool";
for (keys %only) {
  my $n = $_;
  die "No such test: $n\n" unless grep { $_->{name} eq $n } @$tests;
}

$npassed=0;$nfailed=0;$nskipped=0;$nstored=0;
foreach $test (@$tests) {
  $name = $test->{name};
  next if %only and not $only{$name};
  print "\nRunning test: $name....\n";
  system($test->{prepare}) if $test->{prepare};
  system($test->{cmd});
  $getsha = $test->{getsha};
  if ($getsha eq "accuracy") {
    @files = glob($test->{glob});
    if ($opt_s) {
      &store_reference($name,@files);
      $result = "Test $name: reference output stored\n";
      $r = sprintf "Test %-30s reference stored",$name;
      $nstored++;
    } else {
      ($status,$digits,$where) = &compare_to_reference($name,@files);
      if ($status eq "identical") {
        $result = "Test $name passed, output is numerically identical\n";
        $r = sprintf "Test %-30s passed",$name;
        $npassed++;
      } elsif ($status eq "printed") {
        $result = "Test $name passed, output agrees to the printed precision"
          . " of $digits digits, largest difference at\n  $where\n";
        $r = sprintf "Test %-30s OK to printed precision (%d digits)",$name,$digits;
        $npassed++;
      } elsif ($status eq "different") {
        if ($digits >= $mindigits) {
          $r = sprintf "Test %-30s OK to %4.1f significant digits",$name,$digits;
          $npassed++;
        } else {
          $r = sprintf "Test %-30s FAILED: only %4.1f significant digits",$name,$digits;
          $nfailed++;
        }
        $result = sprintf "Test %s: agreement %.1f significant digits, worst at\n  %s\n",
          $name,$digits,$where;
      } else {
        $result = "Test $name failed: $digits\n";
        $r = sprintf "Test %-30s FAILED: %s",$name,$status;
        $nfailed++;
      }
    }
  } else {
    $sha = qx($getsha);
    $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
    if ($sha eq $test->{sha}) {
      $result = "Test $test->{name} passed.\n";
      $r = sprintf "Test %-30s passed",$test->{name};
      $npassed++;
    } elsif (($name eq 'fits-output') and (`./optool -feature fits` =~/False/)) {
      $result = "Test $test->{name} skipped (fits support not implemented).\n";
      $r = sprintf "Test %-30s skipped",$test->{name};
      $nskipped++;
    } else {
      $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
      $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
      $nfailed++;
    }
  }
  print "$result";
  push @results,$r;
}
print "\n";
print "===========================================\n";
print "     Summary of optool test results\n";
print "===========================================\n";
for (@results) {print "$_\n"}
print "-------------------------------------------\n";
print "$nstored references stored; " if $opt_s;
print "$npassed passed; $nfailed failed; $nskipped skipped\n";
print "Tests need at least $mindigits significant digits to pass\n" unless $opt_s;
exit($nfailed ? 1 : 0);

sub store_reference {
  # Store compressed copies of the output files of test NAME as the new
  # reference.  Minimal gzip headers contain no file name or time stamp,
  # so unchanged output gives unchanged reference files.
  my ($name,@files) = @_;
  my $d = "$refdir/$name";
  die "Test $name produced no output files\n" unless @files;
  remove_tree($d);
  make_path($d);
  for my $file (@files) {
    my $ref = "$d/" . basename($file) . ".gz";
    gzip($file => $ref, Minimal => 1, -Level => 9)
      or die "gzip failed for $file: $GzipError\n";
  }
}

sub read_lines {
  # Return a reference to the list of lines in FILE, which may be gzipped.
  my $file = shift;
  my $fh;
  if ($file =~ /\.gz$/) {
    $fh = IO::Uncompress::Gunzip->new($file) or die "Can't open $file\n";
  } else {
    open $fh, '<', $file or die "Can't open $file: $!\n";
  }
  my @lines = <$fh>;
  close $fh;
  chomp @lines;
  return \@lines;
}

sub as_number {
  # Return the numerical value of TOKEN, or undef if it is not a number.
  # Accepts Fortran exponents with D, and the Fortran form 1.234-100
  # used when the exponent needs three digits.
  my $t = shift;
  return $t + 0 if $t =~ /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/;
  return "$1e$2" + 0 if $t =~ /^([-+]?(?:\d+\.?\d*|\.\d+))[dD]([-+]?\d+)$/;
  return "$1e$2" + 0 if $t =~ /^([-+]?\d*\.\d+)([-+]\d+)$/;
  return undef;
}

sub compare_to_reference {
  # Compare the output FILES of test NAME to the reference output.
  # Returns ("identical"), ("different", DIGITS, WHERE),
  # ("printed", DIGITS, WHERE) when the agreement reaches the printed
  # precision of DIGITS significant digits, or
  # (STATUS, MESSAGE) when the files cannot be compared number by number.
  my ($name,@files) = @_;
  my $d = "$refdir/$name";
  my @refs = sort glob("$d/*.gz");
  return ("no reference","no reference output in $d, store it with -s")
    unless @refs;
  my %new = map { basename($_) => $_ } @files;
  my $nref = @refs; my $nnew = @files;
  return ("file count","$nnew output files, but $nref reference files")
    unless $nnew == $nref;
  my $maxrel = 0;
  my $where = "";
  my $prec = 0;      # most significant digits printed in a differing number
  for my $ref (@refs) {
    (my $base = basename($ref)) =~ s/\.gz$//;
    return ("missing file","output file $base not produced") unless $new{$base};
    my $rlines = &read_lines($ref);
    my $nlines = &read_lines($new{$base});
    my ($nr,$nn) = (scalar(@$rlines),scalar(@$nlines));
    return ("structure","$base has $nn lines, reference has $nr") unless $nr == $nn;

    for my $i (0..$nr-1) {
      my @rt = split ' ', $rlines->[$i];
      my @nt = split ' ', $nlines->[$i];
      my $line = $i+1;
      return ("structure","$base line $line: different number of tokens")
        unless @rt == @nt;
      my (@a,@b);
      for my $j (0..$#rt) {
        $a[$j] = &as_number($rt[$j]);
        $b[$j] = &as_number($nt[$j]);
        if (not defined $a[$j] and not defined $b[$j]) {
          return ("text differs","$base line $line: '$rt[$j]' became '$nt[$j]'")
            unless $rt[$j] eq $nt[$j];
        } elsif (not defined $a[$j] or not defined $b[$j]) {
          return ("structure","$base line $line: '$rt[$j]' became '$nt[$j]'");
        }
      }
      my $data = not grep { not defined } @a;
      for my $j (0..$#rt) {
        my ($a,$b) = ($a[$j],$b[$j]);
        next if not defined $a or $a == $b;
        my $scale = 0;
        $scale = 1          if $data and @a == 4 and $j == 3;  # g
        $scale = abs($a[0]) if $data and @a == 6 and $j >= 1;  # F12...F44
        my $norm = abs($a) > abs($b) ? abs($a) : abs($b);
        $norm = $scale if $scale > $norm;
        my $rel = abs($a-$b)/$norm;
        my $p = &printed_digits($rt[$j]);
        $prec = $p if $p > $prec;
        if ($rel > $maxrel) {
          $maxrel = $rel;
          $where = sprintf "%s line %d number %d: reference %s, now %s",
            $base,$line,$j+1,$rt[$j],$nt[$j];
        }
      }
    }
  }
  return ("identical") if $maxrel == 0;
  my $digits = -log($maxrel)/log(10);
  # Agreement beyond the printed precision is not meaningful
  return ("printed", $prec, $where) if $digits >= $prec;
  return ("different", $digits, $where);
}

sub printed_digits {
  # Return the number of significant digits printed in number TOKEN.
  my $t = shift;
  $t =~ s/[eEdD].*$//;                   # remove exponent
  $t =~ s/^([-+]?\d*\.\d+)[-+]\d+$/$1/;  # Fortran form 1.234-100
  $t =~ s/[^0-9]//g;
  $t =~ s/^0+//;
  return length($t);
}
