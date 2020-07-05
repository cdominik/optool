#!/usr/bin/perl

$cmd = "/Applications/Emacs.app/Contents/MacOS/Emacs UserGuide.org --batch -f org-ascii-export-to-ascii --kill";
system($cmd);

$manual = do {
  local $/ = undef;
  open my $fh, "<", "UserGuide.txt" or die "could not open UserGuide.txt: $!";
  <$fh>;
};


@lines = split(/^/,$manual);
$manf90 = "subroutine manual()\n";
foreach (@lines) {
  if ($skip and /^$/) {
    $skip = 0;
    next;
  } elsif ($skip) {
    next;
  } elsif (/^[[]/) {
    $skip = 1;
    next;
  } elsif (/^Footnotes/) {
    last;
  } else {
    chomp;
    if (length($_) >= 110) {
      $_ = substr($_,2,65) . substr($_,91,14);
      s/\*/ /g;
      s/{//g;
      s/}/  /g;
    }
    s/'/ /g;
    s/`/ /g;
    $manf90 .= "  write(*,'(\"$_\")')\n"
  }
}
$manf90 .= "end subroutine manual\n";

print "$manf90\n";
