#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

$dir = 'selftest_optool';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -na 10 -nl 30 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => "eef3f6260dd4a242ef70d042317786d994a5f0e7"},
  {name => 'full_diana',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'd4940b842854a53f47af4c9d7e25ffb6f62083a2'},
  {name => 'radmc-in-42-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -na 42 -nl 100 -s -d -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_*.inp|shasum",
   sha => 'dcfdcc9c8d6ba971bdac3d3d93fe9fcdea581367'}
  ];

foreach $test (@$tests) {
  print "Running test: $test->{name}....\n";
  system($test->{prepare}) if $test->{prepare};
  system($test->{cmd});
  $file = $test->{file};
  $getsha = $test->{getsha};
  $sha = qx($getsha);
  $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
  if ($sha eq $test->{sha}) {
    $result = "Test $test->{name} passed.\n";
    $r = sprintf "Test %-30s passed",$test->{name}
  } else {
    $result = "Test $test->{name} failed: $sha should have been $test->{sha}\n";
    $r = sprintf "Test %-30s FAILED",$test->{name}
  }
  print "$result";
  push @results,$r;
}
print "\n";
print "==========================================\n";
print "     Summary of optool test results\n";
print "==========================================\n";
for (@results) {print "$_\n"}

