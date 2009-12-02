package Bio::Tools::Solubility::Wilkinson;

# ABSTRACT: Calculate the probability of a protein to be soluble using the Wilkinson-Harrison model

use strict;
use warnings;
use Carp 'croak';

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

=method solubility

Returns the probability of the sequence C<$seq> being soluble. C<$seq>
should be a correct and validated protein sequence written in one-letter
aminoacid code. Any whitespace (spaces, newlines, tabs, etc) will be
ignored.

    my $probability = solubility($seq);

This function is exported by default.

=cut

sub solubility {
    my $protein = shift or croak "No protein argument";

    $protein =~ s/\s+//g;

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

__END__

=head1 SYNOPSIS

    use Bio::Tools::Solubility::Wilkinson;

    my $seq = 'MMAEELLVIKP...';

    my $s = solubility($seq);

=cut

=head1 DESCRIPTION

This module implements a simple method for the prediction of protein
solubility, as described by Wilkinson, D.L. and Harrison [1]. It only
takes the primary sequence of the protein as a parameter, and returns
the calculated probability that the protein will be expressed in soluble
form in E. Coli.

It uses aminoacid composition to compute the two main parameters found
to have strong correlation with solubility: charge average and turn
forming residue fraction. For a detailed description of the formula
used, check the article by Harrison [2] in the L<References> section.

=cut

=head1 Customize your imports

If you'd like to rename the only imported subroutine for some reason,
you can do:

    use Bio::Tools::Solubility::Wilkinson
        solubility => { -as => 'solubility_wilkinson' };

    use Some::Other::Solubility::Module 'solubility_foo';

    my $p1 = solubility_wilkinson($seq);
    my $p2 = solubility_foo($seq);

This (and other goodies) are possible thanks to L<Sub::Exporter>, which
is worth checking out.

=cut

=head1 References

=over 4

=item * [1] Wilkinson, D.L. and Harrison, R.G. (1991) BioTechnology 9, 443-448
L<http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=1367308&dopt=Abstract>

=item * [2] R.G. Harrison. 2000. Expression of soluble heterologous proteins via fusion with NusA protein. inNovations. 11:4-7.
L<http://www.biotech.ou.edu/innovations.pdf>

=back

=cut

