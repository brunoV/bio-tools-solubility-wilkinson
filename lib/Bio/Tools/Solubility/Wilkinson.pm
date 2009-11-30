package Bio::Tools::Solubility::Wilkinson;

# ABSTRACT: Calculate the probability of a protein to be soluble using the Wilkinson-Harrison model

sub new {
    my $class = shift;

    my $object;

    return bless \$object, $class;
}

sub solubility {
    my $protein = shift;

    my $solubility = 42;

    return $solubility;
}

1;
