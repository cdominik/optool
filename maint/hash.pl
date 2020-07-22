#!/usr/bin/perl

use Digest::SHA1 qw/ sha1_hex /;

&accuracy_hashes(@ARGV);

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

  foreach $file (@files) {
    print "Doing file $file\n" if $debug;
    open $fh, '<', $file or die "Can't open file $!";
    $content = do { local $/; <$fh> };
    print "c1",$content if $debug;
    $content =~ s/^\s*$//mg; $content =~ s/\n+/\n/mg; # remove empty lines
    print "c2",$content if $debug;
    $content =~ s/\s\s+/ /g; # turn multiple whitespace into single space
    print "c3",$content if $debug;
    push @all,$content;
  }

  $all = join("\n",@all);
  print "all is ",$all if $debug;

  @digits = (1..$maxacc);
  push @digits,$perfect;
  @hashes = ();
  for $i (reverse @digits) {
    print "doing i=$i\n";
    $re = "([-+]?[0-9]+\.[0-9]{$i})([0-9]+)([eE][-+][0-9]+)";
    $all =~ s/$re/"$1" . "0" x length($2) . "$3"/ge if $i<$perfect;
    print $i,":::\n",$all,"\n\n";
    push @hashes,"$i:" . substr(&sha1_hex($all),0,$nh);
  }
  return join("::",@hashes);
}

sub find_accucary {
  # Find the highest accuracy for which the files are identical
  # to the files used to produce the hashes given
  # in the first argument
  my $goalhashes = unshift @_;
  my @files = sort @_;
  my $hashes = &accuracy_hashes(@files);
  my $hashes = "::" . $hashes . "::";
  @goalhashes = split('::',$goalhashes);
  for $goal (@goalhashes) {
    return $goal if $hashes =~ /::${goal}::/;
  }
  return 0;
}
