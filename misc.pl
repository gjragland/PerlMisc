#misc. useful code blocks

#process piped output from external command line by line
open(my $fh, '-|', 'powercfg -l') or die $!;
while (my $line = <$fh>) {
    # Do stuff with each $line.
}
