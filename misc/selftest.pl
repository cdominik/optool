#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

$dir = 'selftest_optool';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '211b9a615d40404164581643b06615f62a957245'},

  {name => 'full-diana-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => 'f7e71ff8cbe9051e107731381169f53ac544d893'},

  {name => 'diana-plus-ice-mantle',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -na 20 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '93e748b947e1ef97575af514cf966871d6448794'},

  {name => 'radmc-in-10-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -s -d 3 -radmc -o $dir",
   getsha => "cat $dir/dustkapscatmat_*.inp|shasum",
   sha => '6b175091013b73195bff49aa4889b11482173a5a'},

  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 10 -na 20 -s 720 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '22edd2f6602fb6db1f7fce5a89e5423ac0a0bd09'},

  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '5a6cb272cd522de512c2b4e8bedeb13b320c87bc'},

  {name => 'chop-peak',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -chop 4 -o $dir",
   getsha => "shasum $dir/dustkapscatmat.dat",
   sha => '29402bbe1277722f938055d3764ad6bdc20f22e5'},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => '05cea43000033be968b1e2c0dbe6081760344a0b'}
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

