use strict;
use warnings;
use Test::More tests => 6;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::ANOVA') };

my $aov = Statistics::ANOVA->new();
isa_ok($aov, 'Statistics::ANOVA');

# Example from Matthew & Delaney (1990) p. 209ff.

my @g1 = (116, 179, 182, 143, 156, 174);
my @g2 = (177, 172, 137, 196, 145, 168);
my @g3 = (170, 156, 188, 212, 164, 184);
my @g4 = (181, 190, 210, 173, 172, 187);
my @g5 = (177, 186, 199, 202, 204, 198);

eval {
    $aov->load_data({1 => \@g1, 2 => \@g2, 3 => \@g3, 4 => \@g4, 5 => \@g5 });
};
ok(!$@, $@);

my %ref_vals = (
	f_value => 3.47,
    sk => 8.86,
);

eval {$aov->anova(independent => 1, parametric => 1, ordinal =>0);};
ok(!$@, $@);

ok( about_equal($aov->{'_stat'}->{'f_value'}, $ref_vals{'f_value'}), "Parametric Clustering:f_value: $aov->{'_stat'}->{'f_value'} = $ref_vals{'f_value'}" );

my $vars;

eval {$vars = $aov->cluster(independent => 1, parametric => 1, ordinal =>0);};
ok(!$@, $@);

sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}


