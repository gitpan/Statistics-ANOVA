 use strict;
 use Statistics::ANOVA 0.03;
 my $varo = Statistics::ANOVA->new();

 # Some data:
 my @gp1 = (qw/8 7 11 14 9/);
 my @gp2 = (qw/11 9 8 11 13/);
 my @gp3 = (qw/7 13 12 8 10/);

 # Load the data, names can be arbitrary
 $varo->load_data({gp1 => \@gp1, gp2 => \@gp2});
 # Oh, forgot one:
 $varo->add_data(gp3 => \@gp3);

 # If they are independent data, test equality of variances, difference between them, and means:
 $varo->obrien_test()->dump(title => 'O\'Brien\'s test of equality of variances');
 $varo->levene_test()->dump(title => 'Levene\'s test of equality of variances');
 $varo->anova_indep()->dump(title => 'Independent groups ANOVA', eta_squared => 1, omega_squared => 1);
 $varo->comparisons_indep();

 # or if they are repeated measures:
 $varo->anova_dep()->dump(title => 'Dependent groups ANOVA');
 $varo->comparisons_dep();
 # or:
 $varo->anova_friedman()->dump(title => 'Friedman test');
 # or:
 $varo->anova_friedman(f_equiv => 1)->dump(title => 'Friedman test');
 
 __END__

##Should print out this:

O'Brien's test of equality of variances
F(2, 12) = 0.386222473178995, p = 0.687770018342956
Levene's test of equality of variances
F(2, 12) = 0.388785046728972, p = 0.686116469325902
Independent groups ANOVA
F(2, 12) = 0.0777777777777777, p = 0.925632513123226, eta-squared = 0.0127970749542962, omega-squared = -0.140202702702703
gp2 - gp1: t(7.17482707174827) = 0.395628284037472 2-p = 0.70416
gp2 - gp3: t(7.48562356676663) = 0.278693205716647 2-p = 0.78854
gp1 - gp3: t(7.94327358676384) = 0.118678165819385 2-p = 0.90886
Bonferroni-adjusted alpha = 0.0166666666666667
Dependent groups ANOVA
F(2, 8) = 0.060475161987041, p = 0.941743285288526
gp2 - gp1: t(4) = 0.399114063142644, 2-p = 0.7102
gp2 - gp3: t(4) = 0.221539510248685, 2-p = 0.83552
gp1 - gp3: t(4) = -0.103417537999004, 2-p = 0.9226
Bonferroni-adjusted alpha = 0.0166666666666667
Friedman test
chi^2(3) = 0.400000000000006, p = 0.81873075307798O'Brien's test of equality of variances
