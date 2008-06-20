package Statistics::ANOVA;

use 5.008008;
use strict;
use warnings;
use Carp qw/carp croak/;
use Statistics::Descriptive;
use Algorithm::Combinatorics qw(combinations);
use Math::Cephes qw(:dists);

use vars qw($VERSION);

our $VERSION = '0.01';

#-----------------------------------------------------------------------
sub new {
#-----------------------------------------------------------------------
                
        my $proto = shift;
        my $class = ref($proto) || $proto;
    
        my $self= {};
        
        foreach (qw/dump/) {
                $self->{$_} = 0;
        }
        ##$self->{$_} = '' foreach qw/df_t df_e f_value chi_value p_value ss_t ss_e title/;
        bless($self, $class);
        return $self;
}

#-----------------------------------------------------------------------
sub load {
#-----------------------------------------------------------------------        
    my $self = shift;
    
    $self->unload();
    
    if (ref $_[0] eq 'HASH') {
      while (my ($sample_name, $sample_data) = each %{$_[0]}) {
         if (ref $sample_data) {
              $self->{'data'}->{$sample_name} = Statistics::Descriptive::Full->new();
              $self->{'data'}->{$sample_name}->add_data(@{$sample_data});
         } 
      }
    }
    else {
       my $sample_name = shift;
       my $sample_data = ref $_[0] eq 'ARRAY' ? $_[0] : scalar (@_) ? \@_ : croak 'No list of data';
       $self->{'data'}->{$sample_name} = Statistics::Descriptive::Full->new();
       $self->{'data'}->{$sample_name}->add_data(@{$sample_data});
    }
}

#-----------------------------------------------------------------------        
sub unload {
#-----------------------------------------------------------------------        
    my $self = shift;
    $self->{'data'} = {};
    $self->{$_} = undef foreach qw/df_t df_e f_value chi_value p_value ss_t ss_e ms_t ms_e/;
}

#-----------------------------------------------------------------------        
sub anova_indep {
#-----------------------------------------------------------------------        
    my ($self, $data) = @_;

    my %data = ref $data eq 'HASH' ? %{$data} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for performing ANOVA';
    
    my $k = scalar(keys(%data));
    if (! $k  || $k == 1) {
        carp 'No groups in the data for performing ANOVA';
        return; 
    }
        
    my ($ss_t, $df_e, $ss_e, $means, @s, @dat, $mean, $count) = ();

    foreach (keys %data) {
        @dat = $data{$_}->get_data;
        $mean = $data{$_}->mean;
        $count = $data{$_}->count;
        croak 'Empty data sent to ANOVA' if !defined $mean || !defined $count;

        $means += $mean;

        # Accumulate within-groups SS, and df:
        foreach (@dat) {
            croak 'Empty data sent to ANOVA' if !defined $_;
            $ss_e += ($_ - $mean)**2;
        }

        $df_e += ($count - 1);
     
        push @s, [$count, $mean]; # store for calculating between SS
    }
    ##croak 'No within-groups for performing ANOVA' if ! $ss_e;
    if (!$ss_e || !$df_e) { carp 'No within-groups for performing ANOVA'; return;}
    # Calc. between groups SS:
    # 1st need the grand mean:
    my $grand_mean = $means / $k;
    foreach (@s) {
        $ss_t += $_->[0] * ($_->[1] - $grand_mean)**2;
    }

    # Calc F, and assoc'd probability:
    my $df_t = $k - 1;
    my $ms_t = $ss_t / $df_t;
    my $ms_e = $ss_e / $df_e;
    my $f = $ms_t / $ms_e;
    my $f_prob = fdtrc($df_t, $df_e, $f);
    #my $f_prob = Statistics::Distributions::fprob($df_t, $df_e, $f);

    $self->{'f_value'} = $f;
    $self->{'p_value'} = $f_prob;
    $self->{'df_t'} = $df_t;
    $self->{'df_e'} = $df_e;
    $self->{'ss_t'} = $ss_t;
	$self->{'ss_e'} = $ss_e;
    $self->{'ms_t'} = $ms_t;
	$self->{'ms_e'} = $ms_e;
        
    return $self;

}

#-----------------------------------------------------------------------        
sub anova_dep {
#-----------------------------------------------------------------------        
    my ($self, $data) = @_;

    my %data = ref $data eq 'HASH' ? %{$data} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for performing ANOVA';
    
    my $k = scalar(keys(%data));
    if (! $k  || $k == 1) {
        carp 'No groups in the data for performing ANOVA';
        return; 
    }

    # Check counts:
    my $count = _check_counts(\%data);
    
    my ($i, $means, @s, @dat, $mean, @p_means, %p_data) = ();
    
    for ($i = 0; $i < $count; $i++) {
        $p_data{$i} = Statistics::Descriptive::Full->new();
        foreach (keys %data) {
            $p_data{$i}->add_data( ($data{$_}->get_data)[$i] );
        }
        $p_means[$i] = $p_data{$i}->mean();
    }
    
    my $grand_mean;
    $grand_mean += $_ foreach @p_means;     
    $grand_mean /= scalar(@p_means);

    # Between groups variance:
    ##my $ss_s;
    #foreach (@p_means) {
    #    $ss_s += ( $_ - $grand_mean)**2;
    #}
    #$ss_s *= $k;
    my $df_b = ( $count - 1);
    #my $ms_b = $ss_s / $df_b;
    
    my $ss_t;
    my %j_means = ();
    foreach (keys %data) {
        $j_means{$_} = $data{$_}->mean();
        $ss_t += ($j_means{$_} - $grand_mean)**2;
    }
    $ss_t *= $count;
    my $df_t = ($k - 1);
    my $ms_t = $ss_t / $df_t;
    
    my $ss_e;
    foreach (keys %data) {
        for ($i = 0; $i < $count; $i++) {
            my $o = ($data{$_}->get_data)[$i];
            $ss_e += ($o - $p_means[$i] - $j_means{$_} + $grand_mean)**2;
        }
    }
    my $df_e = $df_t * $df_b;
    my $ms_e = $ss_e / $df_e;
    
    my $f = $ms_t / $ms_e;
    my $f_prob = fdtrc($df_t, $df_e, $f);
    ##my $f_prob = Statistics::Distributions::fprob($df_t, $df_e, $f);
    
    ##print "Between Ss\t$ss_s\t$df_b\t$ms_b\n";
    #print "Within Ss\n B\t$ss_t $df_t\t$ms_t\t$f\t$f_prob\n BS\t$ss_e\t$df_e\t$ms_e\n";
    ##    $self->{'p_value'} = $self->{'p_precision'} ? sprintf('%.' . $self->{'p_precision'} . 'f', $f_prob) : $f_prob;
    $self->{'f_value'} = $f;
    $self->{'p_value'} = $f_prob;
    $self->{'df_t'} = $df_t;
    $self->{'df_e'} = $df_e;
    $self->{'ss_t'} = $ss_t;
	$self->{'ss_e'} = $ss_e;
    $self->{'ms_t'} = $ms_t;
    $self->{'ms_e'} = $ms_e;
 
    return $self;
}

#-----------------------------------------------------------------------        
sub anova_friedman {
#-----------------------------------------------------------------------        
    my ($self, $data) = @_;
    my %data = ref $data eq 'HASH' ? %{$data} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for performing ANOVA';
    
    my $k = scalar(keys(%data));
    if (! $k  || $k == 1) {
        carp 'No groups in the data for performing ANOVA';
        return; 
    }
    
    # Check counts - should be equal:
    my $count = _check_counts(\%data);
        
    # Get column sum of ranks of row values:
    my ($i, $j, $col, @row_sorted, %ranks, %row_values) = ();
    for ($i = 0; $i < $count; $i++) {
        foreach (keys %data) {
            push @{ $row_values{ ($data{$_}->get_data)[$i] } }, $_;
        }
        @row_sorted = sort {$a <=> $b} keys %row_values;
        for ($j = 0; $j < scalar @row_sorted; $j++) {
            foreach $col( @{$row_values{ $row_sorted[$j]} } ) {
                $ranks{$col} += $j + 1;
            }
        }
        %row_values = ();
    }
    my $sum_squared_ranks;
    foreach (keys %ranks) {
        $sum_squared_ranks += $ranks{$_}**2;
    }
    my $chi = ( ( 12 / ($count * $k * ($k + 1) ) ) * $sum_squared_ranks ) - ( 3 * $count * ($k + 1) );
    my $df = $k - 1;
    my $chi_prob = chdtrc($df, $chi);
    ##use Statistics::Distributions;
    ##my $chi_prob = Statistics::Distributions::chisqrprob ($df, $chi);
    $self->{'chi_value'} = $chi;
    $self->{'p_value'} = $chi_prob;
    $self->{'df_t'} = $k;
    return $self;
}

#-----------------------------------------------------------------------        
sub eta_squared {
#-----------------------------------------------------------------------        
    my ($self) = @_;
    croak 'Need to run ANOVA to obtain requested statistic' if !defined $self->{'df_t'} || !defined $self->{'f_value'};
    $self->{'eta_sq'} = ( $self->{'df_t'} * $self->{'f_value'} ) / ($self->{'df_t'} * $self->{'f_value'} + $self->{'df_e'} );
    ##$ss_t / ( $ss_t / $ss_e );
    return $self->{'eta_sq'};
}

#-----------------------------------------------------------------------        
sub omega_squared {
#-----------------------------------------------------------------------        
    my ($self) = @_;
    my $eta = $self->eta_squared();
    my $a = $eta * ($self->{'df_t'} + $self->{'df_e'}) - $self->{'df_t'};
    my $n;
    foreach (keys %{$self->{'data'}}) {
        $n += $self->{'data'}->{$_}->count;
    } 
    $self->{'omega_sq'} =  $a / ( $a + $n * ( 1 - $eta ) ); 
    return $self->{'omega_sq'};
}

#-----------------------------------------------------------------------        
sub obrien_test {
#-----------------------------------------------------------------------        
    my ($self, $data) = @_;

    my %data = ref $data eq 'HASH' ? %{$data} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for assessing equality of variance';
        
    my ($m, $v, $n, @r) = ();       
    $self->{'obrien'} = {};
    
    # Traverse each sample of data:
    while (my($sample_name, $sample_data) = each %data){
                
        # For each group, compute the sample mean and the unbiased sample variance:
        $m = $sample_data->mean;
        $v = $sample_data->variance;
        $n = $sample_data->count;
    
        # For each observation, compute a transformed score, r:
        my @data = $sample_data->get_data;
        foreach (@data) {
            push @r,
            # Calculate the transformed scores:
                (
                    # Numerator:
                    (
                        (
                           ($n - 1.5 ) * $n * (($_ - $m)**2)
                        )
                        -
                        (
                           .5 * $v * ($n - 1)
                        )
                    )
                    /
                    # Denominator:
                    (
                        ($n - 1) * ($n - 2)
                    )
                );
            }
            $self->{'obrien'}->{$sample_name} = Statistics::Descriptive::Full->new();
            $self->{'obrien'}->{$sample_name}->add_data(@r);
            @r = ();

            # Check that each group mean of the O'Briens are equal to the variance of the original data:
            if (sprintf('%.2f', $self->{'obrien'}->{$sample_name}->mean) != sprintf('%.2f', $v)) {
                croak "mean for sample $sample_name does not equal variance";
            }
        }
        
    ##$self->{'title'} = 'O\'Brien\'s Test for equality of variances:';
    #if ($self->{'p_value'} < .05) {print "Unequal variances by Obrien Test\n"; $self->{'equal_variances'} = 0; } else {$self->{'equal_variances'} = 1;}
    # Perform an ANOVA using the O'Briens as the DV:
    $self->anova_indep($self->{'obrien'});
    
    return $self;
}

#-----------------------------------------------------------------------        
sub levene_test {
#-----------------------------------------------------------------------        
    my ($self, $data) = @_;

    my %data = ref $data eq 'HASH' ? %{$data} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for assessing equality of variance';
        
    my ($m, $v, $n, @d) = ();       
    $self->{'levene'} = {};
    
    # Traverse each sample of data:
    while(my($sample_name, $sample_data) = each %data){
                
        # For each group, compute the sample mean and the unbiased sample variance:
        $m = $sample_data->mean;
        $v = $sample_data->variance;
        $n = $sample_data->count;
                
        # For each observation, compute the absolute deviation:
        my @data = $sample_data->get_data;
        my $m = $sample_data->mean;
        push @d, abs($_ - $m) foreach @data;
        $self->{'levene'}->{$sample_name} = Statistics::Descriptive::Full->new();
        $self->{'levene'}->{$sample_name}->add_data(@d);
        @d = ();
    }

    #$self->{'title'} = 'Levene\'s Test for equality of variances:';   
       
    # Perform an ANOVA using the abs. deviations as the DV:
    $self->anova_indep($self->{'levene'});
    # if ($self->{'p_value'} < .05) {print "Unequal variances by Levene Test\n"; $self->{'equal_variances'} = 0; } else {$self->{'equal_variances'} = 1;}

    return $self;
}

#-----------------------------------------------------------------------        
sub fmax_test {
#-----------------------------------------------------------------------        
    my ($self, $data) = @_;
return;
    my %data = ref $data eq 'HASH' ? %{$data} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for assessing equality of variance';
        
    my (@n, @v) = ();       

    # Traverse each sample of data:
    while(my($sample_name, $sample_data) = each %data){
        push @v, $sample_data->variance;
        push @n, $sample_data->count;
    }
        
    my $vars = Statistics::Descriptive::Sparse->new();
    $vars->add_data(@v);
        
    my $f = ($vars->max / $vars->min);
        
    my $df_t = ($n[$vars->maxdex] - 1);
    my $df_e = ($n[$vars->mindex] - 1);
        
    my $f_prob = fdtrc($df_t, $df_e, $f);
    #my $fprob = Statistics::Distributions::fprob($df_t, $df_e, $f);
    
    $self->{$_} = '' foreach qw/ss_t ss_e ms_t ms_e/;
    
    $self->{'f_value'} = $f;
    $self->{'p_value'} = $f_prob;
    $self->{'df_t'} = $df_t;
    $self->{'df_e'} = $df_e;
    
    return $self;
    #$self->{'title'} = 'F-max test for equality of variances:';
}

#-----------------------------------------------------------------------        
sub comparisons_indep {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    require Statistics::TTest;
    my $ttest = new Statistics::TTest;  
    $ttest->set_significance(95);
    my @all_pairs = combinations([keys %{$self->{'data'}}], 2);
    my $alpha = .05 / scalar(@all_pairs);
    foreach my $pairs (@all_pairs) {
        print "$pairs->[0] - $pairs->[1]: ";
        $ttest->load_data([$self->{'data'}->{$pairs->[0]}->get_data()], [$self->{'data'}->{$pairs->[1]}->get_data]);
        my $p_value = $ttest->{'t_prob'}; # returns the 2-tailed p_value
        my $p_str;
        if ($args{'tails'} && $args{'tails'} == 1) {
            $p_value /= 2;
            $p_str = '1-p';
        }
        else {
            $p_str = '2-p';
        }
        $p_value .= $args{'flag'} ? $p_value < $alpha ? ' *' : '' : '';
        print 't(', $ttest->df, ') = ', $ttest->t_statistic, ' ', $p_str, ' = ', $p_value, "\n";
    }
    print "Bonferroni-adjusted alpha = $alpha\n";
}

#-----------------------------------------------------------------------        
sub comparisons_dep {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    require Statistics::DependantTTest;
    require Statistics::Distributions;
    my $ttest = new Statistics::DependantTTest; 
    my @all_pairs = combinations([keys %{$self->{'data'}}], 2);
    my $alpha = .05 / scalar(@all_pairs);
    foreach my $pairs (@all_pairs) {
        print "$pairs->[0] - $pairs->[1]: ";
        $ttest->load_data($pairs->[0], $self->{'data'}->{$pairs->[0]}->get_data());
        $ttest->load_data($pairs->[1], $self->{'data'}->{$pairs->[1]}->get_data());
        my ($t_value, $deg_freedom) = $ttest->perform_t_test($pairs->[0], $pairs->[1]);
        my $p_value = Statistics::Distributions::tprob($deg_freedom, abs($t_value)); # returns the 1-tailed p_value
        my $p_str;
        if ($args{'tails'} and $args{'tails'} != 2) {
            $p_str = '1-p';
        }
        else {
            $p_value *= 2;
            $p_str = '2-p';
        }
        $p_value .= $args{'flag'} ? $p_value < $alpha ? ' *' : '' : '';
        #$p_value = 1 if $p_value > 1;
        print "t($deg_freedom) = $t_value, $p_str = $p_value\n";
    }
    print "Bonferroni-adjusted alpha = $alpha\n";
}

#-----------------------------------------------------------------------        
sub string {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my $str;
    if (defined $self->{'f_value'}) {
        $str = "F($self->{'df_t'}, $self->{'df_e'}) = $self->{'f_value'}, p = $self->{'p_value'},";
        $str .= ' MSe = ' . $self->{'ms_e'} . ',' if $args{'mse'};
        $str .= ' eta-squared = ' . $self->eta_squared() . ',' if $args{'eta_squared'};
        $str .= ' omega-squared = ' . $self->omega_squared() . ',' if $args{'omega_squared'};
        chop($str);
    }
    elsif (defined $self->{'chi_value'}) {
        $str = "chi^2($self->{'df_t'}) = $self->{'chi_value'}, p = $self->{'p_value'}";
    }
    else {
        croak 'Need to run ANOVA to obtain results string'
    } 
    return $str;
}

#-----------------------------------------------------------------------        
sub dump {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    print "$args{'title'}\n" if $args{'title'};
    print $self->string(%args), "\n";
}

sub _check_counts {
    my $data = shift;
    my $count;
    foreach (keys %{$data}) {
        if (defined $count) {
            my $tmp_count = $data->{$_}->count;
            if ($tmp_count != $count) {
                carp 'Number of observations per group need to be equal for repeated measures ANOVA';
                return;
            }
            else {
                $count = $tmp_count;
            }
        }
        else {
            $count = $data->{$_}->count;
        }
    }
    return $count;
}

# Aliases:
*load_data = \&load;
*anova_rm = \&anova_dep;
*friedman_test = \&anova_friedman;

1;

__END__

=head1 NAME

Statistics::ANOVA - Perform oneway analyses of variance

=head1 SYNOPSIS

 use Statistics::ANOVA;

 my $varo = Statistics::ANOVA->new();

my @data1 = (qw/50	54	55	68	66	76	75	63	62	61	75	63	62	77	68/);
my @data2 = (qw/72	75	65	60	69	64	56	57	55	63	54	63	51	72	78/);
my @data3 = (qw/60	60	70	72	64	79	73	75	72	77	67	83	60	77	70/);
my @data4 = (qw/68	65	76	74	75	83	82	76	85	75	63	64	77	62	83/);

 $varo->load({dat1 => \@data1, dat2 => \@data2, dat3 => \@data3});

 $varo->obrien_test()->dump();
 $varo->levene_test()->dump();
 $varo->fmax_test()->dump();
 $varo->anova_indep()->dump(eta_squared => 1, omega_squared => 1, title => 'Independent groups ANOVA');
 $varo->comparisons_indep();

=head1 DESCRIPTION

Performs oneway between groups and repeated measures ANOVAs, with estimates of proportion of variance acounted for (eta-squared) and effect-size (omega-squared), plus  pairwise comparisons by the relevant t-tests. Also performs equality of variances tests (O'Brien's, Levene's, F-max).

=head1 METHODS

Probabilities for all F-tests are computed using the C<fdtrc> function in the L<Math::Cephes|Math::Cephes> module, rather than the L<Statistics::Distributions|Statistics::Distributions> module, as the former appears to be more accurate for indicating higher-level significances than the latter.

=head2 new

Create a new Statistics::ANOVA object

=head2 load

 $varo->load('data1', @data1)
 $varo->load('data1', \@data1)
 $varo->load({'data1' => \@data1, 'data2' => \@data2})

Alias: C<load_data>

Provided for comparability with other Perl Statistics modules - but you can just send the data you want to test with the call to C<anova_indep> or C<anova_dep> if you like. 

Accepts a single C<name value =E<gt> pair> of a sample name, and a list (referenced or not) of data; or a hash reference of named array references of data. The data are loaded into the class object by name, within a hash called C<data>, as L<Statistics::Descriptive::Full|Statistics::Descriptive> objects. So you could get at the data again, for instance, by going $varo->{'data'}->{'data1'}->get_data(). You can keep updating the data this way, without unloading earlier loads.

Returns the Statistics::ANOVA object.

=head2 unload

 $varo->unload();

Empties all cached data and calculations upon them, ensuring these will not be used for testing. This will not be automatically called with each new load or test. So it should be used whenever switching from one dataset to another. Alternatively, if you send a hash of data directly to C<anova_indep> or C<anova_dep>, that's what will be tested; i.e., you don't have to specifically load and unload the data.

=head2 anova_indep

 $varo->anova_indep()

An implementation of a one-way between-groups analysis of variance. Feeds the class object C<$varo> as follows:

 $varo->{'f_value'}
 $varo->{'df_t'} : the "treatment" or numerator or between-groups degree(s) of freedom
 $varo->{'df_e'} : the "error" or denominator or within-groups degree(s) of freedom
 $varo->{'p_value'} : associated with the F-value with the above dfs
 $varo->{'ss_t'} : treatment sums of squares
 $varo->{'ss_e'} : error sums of squares
 $varo->{'ms_t'} : treatment mean squares
 $varo->{'ms_e'} : error mean squares

=head2 anova_dep

 $varo->anova_dep()

Alias: anova_rm

Performs a one-way repeated measures analysis of variance (sphericity assumed). See C<anova_indep> for fed values.

=head2 anova_friedman

 $varo->anova_friedman()

Alias: friedman_test

Performs Friedman's nonparametric analysis of variance - for two or more dependent (matched, related) groups. The statistical attributes now within the class object (see C<anova_indep>) pertain to this test, e.g., $varo->{'chi_value'} gives the chi-square statistic from the Friedman test; and $varo->{'p_value'} gives the asociated p-value (area under the right-side, upper tail). There is now no defined 'f_value'. See some other module for performing nonparametric pairwise comparisons.

=head2 obrien_test

Performs O'Brien's Test for equality of variances. The statistical attributes now within the class object (see C<anova_indep>) pertain to this test, e.g., $varo->{'f_value'} gives the F-statistic for Obrien's Test; and $varo->{'p_value'} gives the p-value associated with the F-statistic for Obrien's Test.

=head2 levene_test

Performs Levene's (1960) Test for equality of variances. The statistical attributes now within the class object (see C<anova_indep>) pertain to this test, e.g., $varo->{'f_value'} gives the F-statistic for Levene's Test; and $varo->{'p_value'} gives the p-value associated with the F-statistic for Levene's Test.

=head2 fmax_test

I<Presently unimplemented; returns nothing>

Performs the standard F-max test for equality of variances. All statistical attributes now within the class object, $varo, pertain to this test: 'f_value', 'p_value', 'df_t' and 'df_e'.

=head2 comparisons_indep

 $varo->comparisons_indep(tails => 1|2, flag => 1|0 )

Performs independent samples t-tests for each pair of the loaded data, using L<Statistics::TTest|Statistics::TTest>. Simply prints the results to STDOUT. The p_value is 2-tailed, by default, unless otherwise specified, as above.  The output strings are appended with an asterisk if the logical value of the optional attribute C<flag> equals 1 and the C<p_value> is less than the Bonferroni-adjusted alpha level. This alpha level, relative to alpha = .05, for the number of paired comparisons, is printed at the end of the list.

=head2 comparisons_dep

 $varo->comparisons_dep(tails => 1|2, flag => 1|0 )

Performs dependent samples t-tests for each pair of the loaded data, using L<Statistics::DependantTTest|Statistics::DependantTTest>. The number of observations must be equal for each of the data-sets tested. Simply prints the results to STDOUT. The p_value is 2-tailed, by default, unless otherwise specified, as above. The output strings are appended with an asterisk if the logical value of the optional attribute C<flag> equals 1 and the C<p_value> is less than the Bonferroni-adjusted alpha level. This alpha level, relative to alpha = .05, for the number of paired comparisons, is printed at the end of the list.

=head2 eta_squared

Returns eta-squared if an ANOVA has been performed; otherwise C<croak>s. Also feeds $varo with the value, named 'eta_sq'. Values range from 0 to 1, 0 indicating no effect, 1 indicating difference between at least two DV means. Generally indicates the proportion of variance in the DV related to an effect.

=head2 omega_squared

Returns the effect size statistic omega-squared if an ANOVA has been performed; otherwise C<croak>s. Also feeds $varo with the value, named 'omega_sq'. Generally, size is small where omega_sq = .01, medium if omega_sq = .059, and strong if omega_sq = .138. 

=head2 string

 $str = $varo->string(mse => 1, eta_squared => 1, omega_squared => 1)

Returns a statement of result, in the form of C<F(df_t, df_e) = f_value, p = p_value>; or, for Friedman test C<chi^2(df_t) = chi_value, p = p_value>. Optionally also get MSe, eta_squared and omega_squared values appended to the string, where relevant.

=head2 dump

 $varo->dump(mse => 1, eta_squared => 1, omega_squared => 1, title => 'ANOVA test')

Prints the string returned by C<string>. Optionally also get MSe, eta_squared and omega_squared values appended to the string. A newline - "\n" - is appended at the end of the C<print>. Above this string, a title can also be printed, by giving a value to the optional C<title> attribute.

=head1 EXPORT

None by default.

=head1 REFERENCES

None yet.

=head1 SEE ALSO

L<Math::Cephes|lib::Math::Cephes>

=head1 BUGS/LIMITATIONS

No computational bugs as yet identfied. Hopefully this will change, given usage over time. Optimisation of code welcomed.

No adjustment for violations of sphericity in repeated measures ANOVA.

Print only of t-test results.

=head1 REVISION HISTORY

=over 4

=item v 0.01

June 2008

Initital release via PAUSE.

=back

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2008 Roderick Garton

rgarton@utas_DOT_edu_DOT_au

This program is free software. It may may be modified, used, copied, and redistributed at your own risk, and under the terms of the Perl Artistic License (see L<http://www.perl.com/perl/misc/Artistic.html>).
Publicly redistributed modified versions must use a different name.

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut
