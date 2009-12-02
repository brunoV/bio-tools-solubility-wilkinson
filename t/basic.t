use Test::More;
use Test::Exception;
use strict;
use warnings;

BEGIN {
    use_ok q{ Bio::Tools::Solubility::Wilkinson };
}

ok solubility('MAELLKKVIKP'), 'subroutine returns a value';
dies_ok { solubility() }  '... and dies with no arguments';

my $prots = read_sequences('DATA');

cmp_ok
    abs(solubility($_->{seq}) - $_->{solubility}),
    '<', 0.02, "Solubility of $_->{name}" for @$prots;

cmp_ok
    abs(solubility("E" x 15) -
        solubility("E" x 15 . " " x 100 . "\n" x 100)), "<", 0.02,
    "Whitespace doesn't affect the result";

done_testing();

sub read_sequences {
    my $fh = shift;

    my @entries;

    {
        local $/ = ">";
        @entries = <$fh>;
        shift @entries; # discard first empty record
    }

    my @sequences;

    foreach my $entry (@entries) {

        $entry =~ s/[\n|>]//g;
        my ($name, $seq, $solubility) = split /:/, $entry;

        push @sequences, {
            name => $name, seq => $seq, solubility => $solubility
        };
    }

    return \@sequences;
}

__DATA__
>gi|242378714|emb|CAQ33504.1| nusA [Escherichia coli BL21(DE3)]:
MNKEILAVVEAVSNEKALPREKIFEALESALATATKKKYEQEIDVRVQIDRKSGDFDTFRRWLVVDEVTQ
PTKEITLEAARYEDESLNLGDYVEDQIESVTFDRITTQTAKQVIVQKVREAERAMVVDQFREHEGEIITG
VVKKVNRDNISLDLGNNAEAVILREDMLPRENFRPGDRVRGVLYSVRPEARGAQLFVTRSKPEMLIELFR
IEVPEIGEEVIEIKAAARDPGSRAKIAVKTNDKRIDPVGACVGMRGARVQAVSTELGGERIDIVLWDDNP
AQFVINAMAPADVASIVVDEDKHTMDIAVEAGNLAQAIGRNGQNVRLASQLSGWELNVMTVDDLQAKHQA
EAHAAIDTFTKYLDIDEDFATVLVEEGFSTLEELAYVPMKELLEIEGLDEPTVEALRERAKNALATIAQA
QEESLGDNKPADDLLNLEGVDRDLAFKLAARGVCTLEDLAEQGIDDLADIEGLTDEKAGALIMAARNICW
FGDEA:0.95
>gi|169757120|gb|ACA79819.1| thioredoxin [Escherichia coli ATCC 8739]:
MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLNIDQNPGTAPK
YGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA:0.73
>gi|345025|emb|CAA00428.1| Bovine growth hormone [synthetic construct]:
AFPAMSLSGLFANAVLRAQHLHQLAADTFKEFERTYIPEGQRYSIQNTQVAFCFSETMPAPTGKNEAQQK
SDLELLRISLLLIQSWLGPLQFLSRVFTNSLVFGTSDRVYEKLKDLEEGILALMRELEDGTPRRGQILKQ
TYDKFDTNMRSDDALLKNYGLLSCFRKDLHKTETYLRVMKCRRFGEASCAF:0.63
>gi|242378863|emb|CAQ33655.1| bfr [Escherichia coli BL21(DE3)]:
MKGDTKVINYLNKLLGNELVAINQYFLHARMFKNWGLKRLNDVEYHESIDEMKHADRYIERILFLEGLPN
LQDLGKLNIGEDVEEMLRSDLALELDGAKNLREAIGYADSVHDYVSRDMMIEILRDEEGHIDWLETELDL
IQKMGLQNYLQAQIREEG:0.95

