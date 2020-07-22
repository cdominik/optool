#!/usr/bin/perl

@files = @ARGV;

foreach $file (@files) {
  $old = $file;
  $new = $old;
  $new =~ s/\.inp$/.dat/;
  $cmd = "mv $old $new";
  print "$cmd\n";
  system($cmd);
}
    
