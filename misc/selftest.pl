#!/usr/bin/perl

# Take all the link files on the command line and turn them into
# subroutines of a fortaran 90 code spit out to STDOUT. The names of
# the lnk files must follow narrow conventions - see the User guide
# for more information.

# Usage: ingestlnk path/to/lnk/file.lnk [...]

#unless (@ARGV) {
#  &usage();
#  exit(1);
#}

use Data::Dumper;

$dir = 'op_selftest';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat"
   cmd => "./optool -na 10 -nl 30 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => "83cf962afca0d894fca831b01daa5fe504dd161c"},
  {name => 'full_diana',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'e6f6693e7c04288600be6034e92bf2cb42ee44e4'}
  {name => 'radmc-in-42-parts',
   prepare => "rm -f $dir/dustkapscatmat_0{01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42}.inp",
   cmd  => "./optool -na 42 -nl 100 -s -d -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_0{01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42}.inp|shasum",
   sha => 'a9160658c8c77c8f17628913765b9b7cbec242e7'}
  ];

foreach $test (@$tests) {
  print "Running test: $test->{name}....\n";
  system($test->{cmd});
  $file = $test->{file};
  $getsha = $test->{getsha};
  print "getsha is $getsha\n";
  $sha = qx($getsha);
  print "output is: $sha\n";
  $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
  #print "Extracted sha is: $sha\n";
  if ($sha eq $test->{sha}) {
    $result = "Test $test->{name} passed: $sha was indeed equal to $test->{sha}\n";
    $r = sprintf "Test %-30s passed",$test->{name}
  } else {
    $result = "Test $test->{name} failed: $sha should have been $test->{sha}\n";
    $r = sprintf "Test %-30s FAILED",$test->{name}
  }
  print "$result";
  push @results,$r;
}
print "\nSummary of test results:\n";
for (@results) {
  print "$_\n";
}
