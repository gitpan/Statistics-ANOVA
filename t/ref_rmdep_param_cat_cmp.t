use strict;
use warnings;
use Test::More tests => 4; #9;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::ANOVA') };

my $aov = Statistics::ANOVA->new();
isa_ok($aov, 'Statistics::ANOVA');

# Example 1: Maxwell & Delaney (1990), pp. 555ff.
my @y1 = (2, 4, 6, 8, 10, 3, 6, 9);
my @y2 = (3, 7, 8, 9, 13, 4, 9, 11);
my @y3 = (5, 9, 8, 8, 15, 9, 8, 10);

eval {
    $aov->load_data({y1 => \@y1, y2 => \@y2, y3 => \@y3 });
};
ok(!$@, $@);

my %ref_vals = (
    f_12 => .32,
    f_12_adj_e => .23,
    ms_w => 67.375,
);

eval {
    $aov->anova(independent => 0, parametric => 1);
};
ok(!$@, $@);

sub anon { # unimplemented methods set aside for now

ok( about_equal($aov->{'_stat'}->{'ms_w'}, $ref_vals{'ms_w'}), "F-test pair comparison: h1,h2: $aov->{'_stat'}->{'ms_w'} = $ref_vals{'ms_w'}" );

my $pair_dat;
eval {$pair_dat = $aov->compare(independent => 0, parametric => 1, ordinal => 0, adjust_p => 0, adjust_e => 2, use_t => 0);};
ok(!$@, $@);

ok( about_equal($pair_dat->{"h1,h2"}->{'t_value'}, $ref_vals{'f_12_adj_e'}), "F-test pair comparison (variance-adjusted denom.): g1,g2: $pair_dat->{'h1,h2'}->{'t_value'} = $ref_vals{'f_12_adj_e'}" );

eval {$pair_dat = $aov->compare(independent => 1, parametric => 1, ordinal => 0, adjust_p => 0, adjust_e => 0, use_t => 0);};
ok(!$@, $@);

ok( about_equal($pair_dat->{"h1,h2"}->{'t_value'}, $ref_vals{'f_12'}), "F-test pair comparison (unadjusted denom.): g1,g2: $pair_dat->{'h1,h2'}->{'t_value'} = $ref_vals{'f_12'}" );

sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}

} # end sub anon

__END__
