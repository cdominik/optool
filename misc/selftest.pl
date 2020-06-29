#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

$dir = 'selftest_optool';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => "445ebda49c10b282b29ea2b11a2a86629b0ae9cb"},
  {name => 'full-diana-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '8d41f0c1d2603cfaee0b98a0867adbf90b5a88c6'},
  {name => 'diana-plus-ice-mantle',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '30b4e73dc48318a1ee3ff55c78f827d88c1f7a3e'},
  {name => 'radmc-in-21-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 21 -nl 100 -s -d 3 -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_*.inp|shasum",
   sha => 'f10c409d2ce2c9ee84f03a53904668932753869c'},
  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 100 -na 20 -s 720 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'e435cd1f2c792a78f1302e7be1317acb92c9a014'},
  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'f04843fa301db4c3e046710b2a85cafdf7c30d95'},
  {name => 'fits-output',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => '57c3ef30de85d66cd46abdad87ec31cca3a4e879'}
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
    $r = sprintf "Test %-30s passed",$test->{name}
  } else {
    $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
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

