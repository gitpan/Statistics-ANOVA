use strict;
use warnings;
use Test::More tests => 9;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::ANOVA') };

my $aov = Statistics::ANOVA->new();
isa_ok($aov, 'Statistics::ANOVA');

# Example 1: Maxwell & Delaney (), pp. 134-5.
my @h1 = (84, 95, 93, 104);
my @h2 = (81, 84, 92, 101, 80, 108);
my @h3 = (98, 95, 86, 87, 94);
my @h4 = (91, 78, 85, 80, 81);

eval {
    $aov->load_data({h1 => \@h1, h2 => \@h2, h3 => \@h3, h4 => \@h4 });
};
ok(!$@, $@);

my %ref_vals = (
	f_12 => .32,
    f_12_adj_e => .23,
    ms_w => 67.375,
);

eval {
    $aov->anova(independent => 1, parametric => 1);
};
ok(!$@, $@);
ok( about_equal($aov->{'_stat'}->{'ms_w'}, $ref_vals{'ms_w'}), "F-test pair comparison: h1,h2: $aov->{'_stat'}->{'ms_w'} = $ref_vals{'ms_w'}" );

my $pair_dat;
eval {$pair_dat = $aov->compare(independent => 1, parametric => 1, ordinal => 0, adjust_p => 0, adjust_e => 2, use_t => 0);};
ok(!$@, $@);

ok( about_equal($pair_dat->{"h1,h2"}->{'t_value'}, $ref_vals{'f_12_adj_e'}), "F-test pair comparison (variance-adjusted denom.): g1,g2: $pair_dat->{'h1,h2'}->{'t_value'} = $ref_vals{'f_12_adj_e'}" );

eval {$pair_dat = $aov->compare(independent => 1, parametric => 1, ordinal => 0, adjust_p => 0, adjust_e => 0, use_t => 0);};
ok(!$@, $@);

ok( about_equal($pair_dat->{"h1,h2"}->{'t_value'}, $ref_vals{'f_12'}), "F-test pair comparison (unadjusted denom.): g1,g2: $pair_dat->{'h1,h2'}->{'t_value'} = $ref_vals{'f_12'}" );

sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}

__END__

# Example 1: from Gardner (2000) p. 81ff.

my @g1 = (50, 54, 55, 68, 66, 76, 75, 63, 62, 61, 75, 63, 62, 77, 68);
my @g2 = (72, 75, 65, 60, 69, 64, 56, 57, 55, 63, 54, 63, 51, 72, 78);
my @g3 = (60, 60, 70, 72, 64, 79, 73, 75, 72, 77, 67, 83, 60, 77, 70);
my @g4 = (68, 65, 76, 74, 75, 83, 82, 76, 85, 75, 63, 64, 77, 62, 83);


%ref_vals = (
	md_12 => 1.400,
    sig_12 => 1.00,
    md_13 => -5.600,
    sig_13 => .344,
    md_14 => -8.8667,
    sig_14 => .020,
    md_23 => -7.000,
    sig_23 => .111,
    md_24 => -10.2667,
    sig_24 => .005,
    md_34 => -3.2667,
    sig_34 => 1.000,
    flag => ["g1,g4", "g2,g4"],
);

eval {
    $aov->load_data({g1 => \@g1, g2 => \@g2, g3 => \@g3, g4 => \@g4 });
};
ok(!$@, $@);

eval {$pair_dat = $aov->compare(independent => 1, parametric => 1, ordinal => 0, adjust_p => 1, use_t => 0);};
ok(!$@, $@);

ok( about_equal($pair_dat->{"g1,g2"}->{'p_value'}, $ref_vals{'sig_12'}), "T-test comparison: g1,g2: $pair_dat->{'g1,g2'}->{'p_value'} = $ref_vals{'sig_12'}" );

ok( about_equal($pair_dat->{"g1,g3"}->{'p_value'}, $ref_vals{'sig_13'}), "T-test comparison: g1,g3: $pair_dat->{'g1,g3'}->{'p_value'} = $ref_vals{'sig_13'}" );

ok( about_equal($pair_dat->{"g1,g4"}->{'p_value'}, $ref_vals{'sig_14'}), "T-test comparison: g1,g4: $pair_dat->{'g1,g4'}->{'p_value'} = $ref_vals{'sig_14'}" );

ok( about_equal($pair_dat->{"g2,g3"}->{'p_value'}, $ref_vals{'sig_23'}), "T-test comparison: g2,g3: $pair_dat->{'g2,g3'}->{'p_value'} = $ref_vals{'sig_23'}" );

ok( about_equal($pair_dat->{"g2,g4"}->{'p_value'}, $ref_vals{'sig_24'}), "T-test comparison: g2,g4:  $pair_dat->{'g2,g4'}->{'p_value'} = $ref_vals{'sig_24'}" );

ok( about_equal($pair_dat->{"g3,g4"}->{'p_value'}, $ref_vals{'sig_34'}), "T-test comparison g3,g4: $pair_dat->{'g3,g4'}->{'p_value'} = $ref_vals{'sig_34'}" );





