#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

$dir = 'selftest_optool';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -na 10 -nl 30 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => "83cf962afca0d894fca831b01daa5fe504dd161c"},
  {name => 'full_diana',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'e6f6693e7c04288600be6034e92bf2cb42ee44e4'},
  {name => 'radmc-in-42-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -na 42 -nl 100 -s -d -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_*.inp|shasum",
   sha => '4a4bfe26ec344e4048dcf7d4f3ef8abb1f1999bc'}
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

