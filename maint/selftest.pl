#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.
# When run with -s, the hashes at the end of the file are set.  You may
# want to clean out the old velues by hand.

use Digest::SHA qw/ sha1_hex /;
use Getopt::Std;
getopts('s');
$dir = 'selftest_optool';

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

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => '392ccf20ccb2a74cc60305b2b12c5da6deac70d0'}
  ];

# Read the hashes
while (<DATA>) {
  next if /^\s*#/;
  next if /^\s*$/;
  ($name,$hashes) = split(/:::/,$_);
  $hashes{$name} = $hashes;
}

foreach $test (@$tests) {
  $name = $test->{name};
  #next unless $test->{name} eq "noscat";
  print "\nRunning test: $name....\n";
  system($test->{prepare}) if $test->{prepare};
  system($test->{cmd});
  $getsha = $test->{getsha};
  if ($getsha eq "accuracy") {
    $glob = $test->{glob};
    @files = glob($glob);
    #print "files are: @files\n"; 
    if ($opt_s) {
      $hashes = &accuracy_hashes(@files);
      $allhashes .= "${name}:::${hashes}\n";
    } else {
      $acc = &find_accuracy($hashes{$name},@files);
      #print "accuracy found is $acc\n";
      if ($acc == 1000) {
        $result = "Test $test->{name} perfect pass\n";
        $r = sprintf "Test %-30s passed",$test->{name};
      } elsif ($acc == 0) {
        $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
        $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
      } else {
        $acc = $1+1 if $acc =~ /([0-9]+):/;
        $result = "Test $test->{name} passed with $acc significant digits\n";
        $r = sprintf("Test %-30s OK to $acc significant digits",$test->{name})
      }
    }
  } else {
    $sha = qx($getsha);
    $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
    if ($sha eq $test->{sha}) {
      $result = "Test $test->{name} passed.\n";
      $r = sprintf "Test %-30s passed",$test->{name};
    } else {
      $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
      $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
    }
  }
  print "$result";
  push @results,$r;
}
if ($opt_s) {
  $date = `date`;
  open SELF,">> $0" or die "Unable to open SELF\n";
  print SELF "\n\n# HASHES set on $date\n",$allhashes,"\n";
} else {
  print "\n";
  print "==========================================\n";
  print "     Summary of optool test results\n";
  print "==========================================\n";
  for (@results) {print "$_\n"}
}

sub accuracy_hashes {
  # Create a series of hashes of the concatenation of all
  # files in @_, where the accuracy of number with exponential format
  # is limited to N digits after the comma.
  # The return value is a string with 10 different hashes,
  # accuracy 10 to 1. The hashes have the format ACC:HASH where ACC
  # is an integer and HASH is the first 10 digits of the sha1 hash of
  # the concatenated files with whitespace compressed and accuracy limited.
  # The different hashes are joined with "::" into a single string.
  my @files = sort @_;
  my $all = "";
  my @all = ();
  my @numbers = ();
  my $file,$content,$fh;
  my $nh = 10;            # Number of chars from sha1 hash to use
  my $debug = 0;
  my $maxacc = 10;        # Maximum accuracy for which to produce a hash
  my $perfect = 1000;     # "accuracy: key for perfect, i.e. unmodified
  my $scatlim = 1e-2;
  
  foreach $file (@files) {
    # print "Doing file $file\n";
    open $fh, '<', $file or die "Can't open file $!";
    $content = do { local $/; <$fh> };

    #$content =~ s/[-+][0-9]\.[0-9]+[eE]-[0-9]+/1.234567E-01/g;
    $content =~ s/[\n\t ]+/ /g; # turn multiple whitespace into single space
    $content =~ s/\s+$//;
    $content =~ s/^\s+//;
    print "c3",$content if $debug;
    push @all,$content;
  }

  $all = join("\n",@all);
  print "all is ",$all if $debug;

  @digits = (1..$maxacc);
  push @digits,$perfect;
  @hashes = ();
  for $i (reverse @digits) {
    $re = "([-+]?[0-9]+[.][0-9]{$i})([0-9]+)([eE][-+][0-9]+)";
    #print "re is $re\n";
    $all =~ s/$re/"$1" . "0" x length($2) . "$3"/ge if $i<$perfect;
    #print $i,":::\n",$all,"\n\n";
    push @hashes,"$i:" . substr(&sha1_hex($all),0,$nh);
  }
  return join("::",@hashes);
}

sub find_accuracy {
  # Find the highest accuracy for which the files are identical
  # to the files used to produce the hashes given
  # in the first argument
  #print "@_";
  my $goalhashes = shift @_;
  my @files = sort(@_);
  my $hashes = &accuracy_hashes(@files);
  $hashes = "::" . $hashes . "::";
  @goalhashes = split('::',$goalhashes);
  for $goal (@goalhashes) {
    return $goal if $hashes =~ m/::${goal}::/;
  }
  return 0;
}

__DATA__

# HASHES set on Thu Nov 10 19:25:17 CET 2022

noscat:::1000:18843f03d4::10:18843f03d4::9:18843f03d4::8:18843f03d4::7:18843f03d4::6:18843f03d4::5:c395b90667::4:b8b7c8aaec::3:0439e7d021::2:be78c20493::1:735d6b4bfc
quick:::1000:22fb28ab58::10:22fb28ab58::9:22fb28ab58::8:22fb28ab58::7:22fb28ab58::6:22fb28ab58::5:18b8099297::4:0bbcd1dba4::3:d3b9253c73::2:4ca4320855::1:19b3b43e22
full-diana-lowres:::1000:040f64207e::10:040f64207e::9:040f64207e::8:040f64207e::7:040f64207e::6:040f64207e::5:8da5320161::4:a535432bf4::3:548805a7a7::2:f8466c6a3e::1:eb52685017
sparse-file:::1000:19c7554cb6::10:19c7554cb6::9:19c7554cb6::8:19c7554cb6::7:19c7554cb6::6:19c7554cb6::5:c125475e08::4:a1b2f353e8::3:701d3d1cca::2:1376d37107::1:6eb3b3d969
lognormal-lowres:::1000:0657684c37::10:0657684c37::9:0657684c37::8:0657684c37::7:0657684c37::6:0657684c37::5:d95a419365::4:e5839a581e::3:1ebfb49e26::2:feb2f63812::1:62ce857d5f
diana-plus-ice-mantle:::1000:55a5f43520::10:55a5f43520::9:55a5f43520::8:55a5f43520::7:55a5f43520::6:55a5f43520::5:96ec0c5335::4:77f3f4c09d::3:7bcf3b947f::2:32a1adb475::1:2d7fb84780
div-in-10-parts:::1000:c8f8d2f52b::10:c8f8d2f52b::9:c8f8d2f52b::8:c8f8d2f52b::7:c8f8d2f52b::6:c8f8d2f52b::5:89dac49573::4:cd0c9bd56c::3:c010cba240::2:f507b33e09::1:cf6e9a43cb
large-grain:::1000:a72be7b67b::10:a72be7b67b::9:a72be7b67b::8:a72be7b67b::7:a72be7b67b::6:a72be7b67b::5:2a84b0dc11::4:2b38359836::3:e6cac4175e::2:7a1549eaa6::1:2c11eb3264
high-angular-resolution:::1000:9a9abff575::10:9a9abff575::9:9a9abff575::8:9a9abff575::7:9a9abff575::6:9a9abff575::5:eeb1b5ec7e::4:d15d774d4a::3:27b0cb074b::2:0e616238e6::1:837cb8ba33
pure-ice-grain:::1000:1f1dc3247e::10:1f1dc3247e::9:1f1dc3247e::8:1f1dc3247e::7:1f1dc3247e::6:1f1dc3247e::5:4cd021606f::4:2a2008a9cd::3:bf7aa7396f::2:81724ae0fd::1:4446b78e80
chop-peak:::1000:85e7873eb4::10:85e7873eb4::9:85e7873eb4::8:85e7873eb4::7:85e7873eb4::6:85e7873eb4::5:a24f8e0479::4:aa681ac7fb::3:584f6878e2::2:56cba6bf0f::1:cb4d2b5829
cde-rayleigh:::1000:503e063b60::10:503e063b60::9:503e063b60::8:503e063b60::7:503e063b60::6:503e063b60::5:de3c2965c1::4:950b7e9122::3:1a286b7ca5::2:1e44bf7ab9::1:a8f7f57541
mmf-scat:::1000:a919f27210::10:a919f27210::9:a919f27210::8:a919f27210::7:a919f27210::6:a919f27210::5:cacdd11446::4:aa09fa0c4e::3:c6930562d4::2:f34ba2591e::1:9cb68b3ded
mmf-opac:::1000:c1637cf821::10:c1637cf821::9:c1637cf821::8:c1637cf821::7:c1637cf821::6:c1637cf821::5:2f79f583e7::4:7f738c463f::3:e03c264ce5::2:3a97e5e2a4::1:5e4701113e

