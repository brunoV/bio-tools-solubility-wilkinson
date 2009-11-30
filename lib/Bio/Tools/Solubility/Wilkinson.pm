package Bio::Tools::Solubility::Wilkinson;

# ABSTRACT: Calculate the probability of a protein to be soluble using the Wilkinson-Harrison model

use strict;
use warnings;
use Sub::Exporter -setup => {
    exports => ['solubility'],
    groups  => { default => ['solubility'] },
};

use constant {
    L1  =>  15.43,
    L2  => -29.56,
    CVp =>   1.71,
    A   =>   0.4934,
    B   =>   0.276,
    C   =>  -0.0392,
};

sub solubility {
    my $protein = shift // die "No protein argument";

    my $CV      = _CV($protein);
    my $CV_norm = abs($CV - CVp);

    my $probability = A + B * $CV_norm + C * ($CV_norm**2);

    return $probability;
}

sub _CV {
    my $protein = shift;

    my %n = map { $_ => _aa_count($protein, $_) } qw(N G P S R K D E);
    my $n = length $protein;

    my $cv =
        L1 * (   ($n{N} + $n{G} + $n{P} + $n{S}) / $n       )
      + L2 * abs(($n{R} + $n{K} - $n{D} - $n{E}) / $n - 0.03);

    return $cv;

}

sub _aa_count {
    my ($protein, $aa) = @_;

    my @occurences = $protein =~ /$aa/ig;

    return scalar @occurences;

}

1;
