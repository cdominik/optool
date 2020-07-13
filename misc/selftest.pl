#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

use Digest::SHA1 qw/ sha1_hex /;

$dir = 'selftest_optool';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'df939aa359c7530c3167fb6eeba7e97f750f24b3'},

  {name => 'full-diana-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'c0b89c3f564c441c795c3311b815486b37ad2ff5'},

  {name => 'diana-plus-ice-mantle',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -na 20 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '48c6b0cae1287916a3f27202c76f29b098b31894'},

  {name => 'radmc-in-10-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -s -d 3 -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_*.inp|shasum",
   sha => 'dded4f829be86c6b343d53c01f5702d8d3bf0259'},

  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 10 -na 20 -s 720 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '823bea94e4a3d79fd155ea8fa41907cbac02f164'},

  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '3732c66467f2e7503f0223920131a62369d25908'},

  {name => 'chop-peak',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -chop 4 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '598aa09b0da59e3f790827fd7649c1c92fd86925'},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => '5db1e1afef60ffdfd94b9e4cb7a9626790d7a501'}
  ];

foreach $test (@$tests) {
  print "\nRunning test: $test->{name}....\n";
  system($test->{prepare}) if $test->{prepare};
  system($test->{cmd});
  $file = $test->{file};
  $getsha = $test->{getsha};
  $sha = qx($getsha);
  $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
  if ($sha eq $test->{sha}) {
    $result = "Test $test->{name} passed.\n";
    $r = sprintf "Test %-30s passed",$test->{name};
    # system('say Great!')
  } else {
    $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
    $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
    # make selftest
    # system('say Ooops!')
  }
  print "$result";
  push @results,$r;
}
print "\n";
print "==========================================\n";
print "     Summary of optool test results\n";
print "==========================================\n";
for (@results) {print "$_\n"}



