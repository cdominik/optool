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
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -na 20 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '8d49dd45b2b5f0f8c46d451b4c9f2885ee1f07e4'},

  {name => 'radmc-in-10-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -s -d 3 -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_*.inp|shasum",
   sha => '79f75953c67805209e1aaeb31174f8705b3a7b78'},

  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 10 -na 20 -s 720 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '790bc492a96d4d2f60a41d1582f09b561977d9c5'},

  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'f04843fa301db4c3e046710b2a85cafdf7c30d95'}
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

