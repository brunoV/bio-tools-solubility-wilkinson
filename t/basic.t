use Test::More;

use_ok 'Bio::Tools::Solubility::Wilkinson';

my $s = Bio::Tools::Solubility::Wilkinson->new;

isa_ok $s, 'Bio::Tools::Solubility::Wilkinson';

can_ok $s, 'solubility';

done_testing();
