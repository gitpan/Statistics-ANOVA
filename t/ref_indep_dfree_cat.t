use strict;
use warnings;
use Test::More tests => 5;
use constant EPS => 1e-9;

BEGIN { use_ok('Statistics::ANOVA') };

my $aov = Statistics::ANOVA->new();
isa_ok($aov, 'Statistics::ANOVA');

# Example for Siegal (1956) p. 190.

my @g1 = (2.0, 2.8, 3.3, 3.2, 4.4, 3.6, 1.9, 3.3, 2.8, 1.1);
my @g2 = (3.5, 2.8, 3.2, 3.5, 2.3, 2.4, 2.0, 1.6);
my @g3 = (3.3, 3.6, 2.6, 3.1, 3.2, 3.3, 2.9, 3.4, 3.2, 3.2);
my @g4 = (3.2, 3.3, 3.2, 2.9, 3.3, 2.5, 2.6, 2.8);
my @g5 = (2.6, 2.6, 2.9, 2.0, 2.0, 2.1);
my @g6 = (3.1, 2.9, 3.1, 2.5);
my @g7 = (2.6, 2.2, 2.2, 2.5, 1.2, 1.2);
my @g8 = (2.5, 2.4, 3.0, 1.5);

eval {
    $aov->load_data({g1 => \@g1, g2 => \@g2, g3 => \@g3, g4 => \@g4, g5 => \@g5, g6 => \@g6, g7 => \@g7, g8 => \@g8});
};
ok(!$@, $@);

my %ref_vals = (
	h_value => 18.5654138029782,
	p_value => 0.00966343126440618,
);

$aov->anova(independent => 1, parametric => 0, ordinal =>0);

ok( about_equal($aov->{'_stat'}->{'h_value'}, $ref_vals{'h_value'}), "Kruskal-Wallis H-value $aov->{'_stat'}->{'h_value'} = $ref_vals{'h_value'}" );

ok( about_equal($aov->{'_stat'}->{'p_value'}, $ref_vals{'p_value'}), "Kruskal-Wallis p-value $aov->{'_stat'}->{'p_value'} = $ref_vals{'p_value'}" );
 
sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}


