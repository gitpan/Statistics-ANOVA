package Statistics::ANOVA;

use 5.008008;
use strict;
use warnings;
use Carp qw/croak/;
use Algorithm::Combinatorics qw(combinations);
use Math::Cephes qw(:dists);
use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;
use vars qw($VERSION);

our $VERSION = '0.06';

#-----------------------------------------------------------------------
sub new {
#-----------------------------------------------------------------------
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self= {};
    ##$self->{$_} = '' foreach qw/df_t df_e f_value chi_value p_value ss_t ss_e title/;
    bless($self, $class);
    return $self;
}

#-----------------------------------------------------------------------
sub load {
#-----------------------------------------------------------------------        
    my $self = shift;
    $self->unload();
    $self->add(@_);
}

#-----------------------------------------------------------------------
sub add {
#-----------------------------------------------------------------------        
    my $self = shift;
    
    if (ref $_[0] eq 'HASH') {
      while (my ($sample_name, $sample_data) = each %{$_[0]}) {
         if (ref $sample_data) {
              $self->{'data'}->{$sample_name} = Statistics::Descriptive::Full->new();
              my $dat = $self->_purge($sample_data);
              $self->{'data'}->{$sample_name}->add_data($dat);
         } 
      }
    }
    else {
       my $sample_name = shift;
       my $sample_data = ref $_[0] eq 'ARRAY' ? $_[0] : scalar (@_) ? \@_ : croak 'No list of data';
       $self->{'data'}->{$sample_name} = Statistics::Descriptive::Full->new();
       my $dat = $self->_purge($sample_data);
       $self->{'data'}->{$sample_name}->add_data($dat);
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
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to a hash of data for performing ANOVA';
    
    my $k = scalar(keys(%data));
    croak 'Not enough groups, if any, in the data for performing ANOVA' if ! $k  || $k == 1;
        
    my ($ss_t, $df_e, $ss_e, $sum_means, @s, @dat, $mean, $count) = ();
    my $all = Statistics::Descriptive::Sparse->new();
    foreach (keys %data) {
        @dat = $data{$_}->get_data;
        $mean = $data{$_}->mean;
        $count = $data{$_}->count;
        $all->add_data(@dat);
        $sum_means += $mean;

        # accumulate within-groups SS, and df:
        $ss_e += ($_ - $mean)**2 foreach @dat; 
        $df_e += ($count - 1);

        push @s, [$count, $mean]; # for calculating between SS
    }
    croak 'No within-groups for performing ANOVA' if !$ss_e || !$df_e;

    # Calc. between groups SS:     # 1st need grand mean:
    my $grand_mean = $all->mean();
    foreach (@s) { # - arefs of counts and means per group
        $ss_t += $_->[0] * ($_->[1] - $grand_mean)**2; # M&D Eq. 57
    }

    # Calc F & prob:
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
    $self->{'nparam'} = 0;
    return $self;

}

#-----------------------------------------------------------------------        
sub anova_dep {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to a hash of data for performing ANOVA';

    my $k = scalar(keys(%data));
    if (! $k  || $k == 1) {
        croak 'Not enough groups, if any, in the data for performing ANOVA';
        return; 
    }

    # Check counts:
    my $count = _check_counts(\%data);
    #my $all = Statistics::Descriptive::Sparse->new();

    my ($i, $grand_mean, @i_means, %i_data) = ();

    for ($i = 0; $i < $count; $i++) {
        $i_data{$i} = Statistics::Descriptive::Full->new();
        foreach (keys %data) {
            #$all->add_data(($data{$_}->get_data)[$i]);
            $i_data{$i}->add_data( ($data{$_}->get_data)[$i] );
        }
        $i_means[$i] = $i_data{$i}->mean();
    }

    $grand_mean += $_ foreach @i_means;     
    $grand_mean /= scalar(@i_means);
    #$grand_mean = $all->mean(); # rely on above as s/be equal obs per gp

    # Between groups variance:
    ##my $ss_s;
    #foreach (@i_means) {
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
            $ss_e += ($o - $i_means[$i] - $j_means{$_} + $grand_mean)**2;
        }
    }

    my $df_e = $df_t * $df_b;
    my $ms_e = $ss_e / $df_e;

    my $f = $ms_t / $ms_e;    
    my $f_prob = fdtrc($df_t, $df_e, $f);
    ##my $f_prob = Statistics::Distributions::fprob($df_t, $df_e, $f);

    $self->{'f_value'} = $f;
    $self->{'p_value'} = $f_prob;
    $self->{'df_t'} = $df_t;
    $self->{'df_e'} = $df_e;
    $self->{'ss_t'} = $ss_t;
	$self->{'ss_e'} = $ss_e;
    $self->{'ms_t'} = $ms_t;
    $self->{'ms_e'} = $ms_e;
    $self->{'nparam'} = 0;
    return $self;
}

#-----------------------------------------------------------------------        
sub anova_friedman {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to a hash of data for performing ANOVA';
    
    my $k = scalar(keys(%data));
    if (! $k  || $k == 1) {
        croak 'Not enough groups, if any, in the data for performing ANOVA';
    }
    
    # Check counts - should be equal:
    my $count = _check_counts(\%data);
  
    # Get column sum of ranks of row values:
    my ($i, $col, %ranks, %row_values) = ();
    my ($old, $cur) = (0, 0);

    # Set the averaged ranks.
    my @ranks;
    # - go down each index:
    for ($i = 0; $i < $count; $i++) {
        # - list the values at this index in each data-array (being Statistics-Descriptive objects):
        # - a value might occur in more than one group at this index, so store an array of the groups:
        push @{ $row_values{ ($data{$_}->get_data)[$i] } }, $_ foreach keys %data;
        # This loop adapted from Gene Boggs' "rank" function in Statistics-RankCorrelation:
        for my $x (sort { $a <=> $b } keys %row_values) {
            # Get the number of ties:
            my $ties = scalar(@{ $row_values{$x} });
            $cur += $ties;
            if ($ties > 1) {
                # Average the tied data:
                my $average = $old + ($ties + 1) / 2;
                $ranks{$_} += $average for @{ $row_values{$x} };
            }
            else {
                # Add the single rank to the list of ranks:
                $ranks{ $row_values{$x}[0] } += $cur;
            }
            $old = $cur;
        }

        ($old, $cur) = (0, 0);
        %row_values = ();
    }
    
    my $sum_squared_ranks;
    foreach (keys %ranks) {
        $sum_squared_ranks += $ranks{$_}**2;
    }
    # x=[12/(N k (k+1))]Sum(1 ... k)(Rj^2) - 3 N (k+1) where Rj is the sum of ranks for j (j=1 ... k)
    my $chi = ( 12 / ($count * $k * ($k + 1) ) ) * $sum_squared_ranks - 3 * $count * ($k + 1) ;
    my $df = $k - 1;
    if ($args{'f_equiv'}) {
        my $f = (($count - 1) * $chi) / ($count * ($k - 1) - $chi);
        my $df_t = $k - 1;
        my $df_e = ($count - 1) * ($df_t);
        my $f_prob = fdtrc($df_t, $df_e, $f);
        $self->{'f_value'} = $f;
        $self->{'p_value'} = $f_prob;
        $self->{'df_t'} = $df_t;
        $self->{'df_e'} = $df_e;
        $self->{'nparam'} = 0;
    }
    else {
        my $chi_prob = chdtrc($df, $chi);
        ##use Statistics::Distributions;
        ##my $chi_prob = Statistics::Distributions::chisqrprob ($df, $chi);
        $self->{'chi_value'} = $chi;
        $self->{'p_value'} = $chi_prob;
        $self->{'df_t'} = $k;
        $self->{'nparam'} = 1;
    }
    return $self;
}

#-----------------------------------------------------------------------        
sub obrien_test {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to a hash of data for performing ANOVA';
        
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
                croak "Mean for sample $sample_name does not equal variance";
            }
        }

    # Perform an ANOVA using the O'Briens as the DV:
    $self->anova_indep(data => $self->{'obrien'});
    
    return $self;
}

#-----------------------------------------------------------------------        
sub levene_test {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for performing ANOVA';
       
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

    # Perform an ANOVA using the abs. deviations as the DV:
    $self->anova_indep(data => $self->{'levene'});

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
    print "Adjusted alpha = $alpha\n";
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
    print "Adjusted alpha = $alpha\n";
}

#-----------------------------------------------------------------------        
sub string {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my $str;
    my $p = $args{'p_precision'} ? sprintf('%.' . $args{'p_precision'} . 'f', $self->{'p_value'}) : $self->{'p_value'};
    my $s_precision = $args{'s_precision'} || 3;
    if (defined $self->{'f_value'} && !$self->{'nparam'}) {
        $str .= "F($self->{'df_t'}, $self->{'df_e'}) = ";
        $str .= sprintf('%.' . $s_precision . 'f', $self->{'f_value'});
        $str .= ", p = $p,";
        $str .= ' MSe = ' . sprintf('%.' . $s_precision . 'f', $self->{'ms_e'}) . ',' if $args{'mse'};
        $str .= ' eta^2 = ' . sprintf('%.' . $s_precision . 'f', $self->eta_squared()) . ',' if $args{'eta_squared'};
        $str .= ' omega^2 = ' . sprintf('%.' . $s_precision . 'f', $self->omega_squared()) . ',' if $args{'omega_squared'};
        chop($str);
    }
    elsif (defined $self->{'chi_value'}) {
        $str .= "chi^2($self->{'df_t'}) = ";
        $str .= sprintf('%.' . $s_precision . 'f', $self->{'chi_value'});
        $str .= ", p = $p";
    }
    else {
        croak 'Need to run ANOVA to obtain results string'
    } 
    return $str;
}

#-----------------------------------------------------------------------        
sub table {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    my $tbl;
    my $p = $args{'p_precision'} ? sprintf('%.' . $args{'p_precision'} . 'f', $self->{'p_value'}) : $self->{'p_value'};
    my $s_precision = $args{'s_precision'} || 3;
    # F-table:
    if (defined $self->{'f_value'} && !$self->{'nparam'}) {
        $tbl .= "\t$_" foreach ('Df', 'Sum Sq', 'Mean Sq', 'F value', 'Pr(>F)');
        $tbl .=  "\n";
        $tbl .=  "$_\t" foreach ('Factor', $self->{'df_t'}); 
        $tbl .= sprintf('%.' . $s_precision . 'f', $_) . "\t" foreach ($self->{'ss_t'}, $self->{'ms_t'}, $self->{'f_value'});
        $tbl .= $p;
        $tbl .=  "\n";
        $tbl .=  "$_\t" foreach ('Error', $self->{'df_e'}); 
        $tbl .= sprintf('%.' . $s_precision . 'f', $_) . "\t" foreach ($self->{'ss_e'}, $self->{'ms_e'});
        $tbl .=  "\n";
    }
   
    return $tbl;
}

#-----------------------------------------------------------------------        
sub dump {
#-----------------------------------------------------------------------        
    my ($self, %args) = @_;
    print "$args{'title'}\n" if $args{'title'};
    if ($args{'table'}) {
        print $self->table(%args);
        print $self->string(%args), "\n" if $args{'string'};
    }
    else {
        print $self->string(%args), "\n";
    }
    print "Observations purged as undefined or not-a-number: ". $self->{'purged'} . "\n" if $self->{'purged'} && $args{'verbose'};
}

#-----------
# private functions - do drop in

sub _check_counts {
    my $data = shift;
    my $count;
    foreach (keys %{$data}) {
        if (defined $count) {
            my $tmp_count = $data->{$_}->count;
            if ($tmp_count != $count) {
                croak 'Number of observations per group need to be equal for repeated measures ANOVA';
                ##return;
            }
            else {
                $count = $tmp_count;
            }
        }
        else {
            $count = $data->{$_}->count;
        }
    }
    croak 'One observation per group can not be tested in ANOVA' if $count == 1;
    return $count;
}

sub _purge {
   my ($self, $dat) = @_;
   my @true = grep { looks_like_number($_) } @{$dat};
   $self->{'purged'} += ( scalar(@{$dat}) - scalar(@true) );
   croak "Empty data sent to ANOVA" if !scalar(@true);
   return \@true;
}

# Aliases:
*load_data = \&load;
*add_data = \&add;
*anova_rm = \&anova_dep;
*friedman_test = \&anova_friedman;

1;

__END__

=head1 NAME

Statistics::ANOVA - Perform oneway analyses of variance

=head1 SYNOPSIS

 use Statistics::ANOVA 0.06;
 my $varo = Statistics::ANOVA->new();

 # Some data:
 my @gp1 = (qw/8 7 11 14 9/);
 my @gp2 = (qw/11 9 8 11 13/);

 # Load the data (names can be arbitrary):
 $varo->load_data({gp1 => \@gp1, gp2 => \@gp2});
 # Oh, here comes another one:
 my @gp3 = (qw/7 13 12 8 10/);
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
 # or 
 $varo->anova_friedman(f_equiv => 1)->dump(title => 'Friedman test');

=head1 DESCRIPTION

Performs oneway between groups and repeated measures ANOVAs, with estimates of proportion of variance acounted for (eta-squared) and effect-size (omega-squared), plus  pairwise comparisons by the relevant t-tests. Also performs equality of variances tests (O'Brien's, Levene's).

=head1 METHODS

=head2 new

Create a new Statistics::ANOVA object

=head2 load

 $varo->load('aname', @data1)
 $varo->load('aname', \@data1)
 $varo->load({'aname' => \@data1, 'another_name' => \@data2})

I<Alias>: C<load_data>

Accepts either (1) a single C<name =E<gt> value> pair of a sample name, and a list (referenced or not) of data; or (2) a hash reference of named array references of data. The data are loaded into the class object by name, within a hash named C<data>, as L<Statistics::Descriptive::Full|Statistics::Descriptive> objects. So you can easily get at any descriptives for the groups you've loaded - e.g., $varo->{'data'}->{'aname'}->mean() - or you could get at the data again by going $varo->{'data'}->{'aname'}->get_data(); and so on. The names of the data are up to you.

I<Missing/Invalid values>: Any observations that are undefined or not-a-number are purged prior to being shunted off to L<Statistics::Descriptive::Full|Statistics::Descriptive> (which does not handle missing or invalid values itself). The number of such purged values are cached thus: $varo->{'purged'}. The L<dump|dump> method can also reveal this value. The C<looks_like_number> method in L<Scalar::Util|Scalar::Util/looks_like_number> is used for this purpose.

Each call L<unload|unload>s any previous loads.

Returns the Statistics::ANOVA object.

=head2 add

 $varo->add('another_name', @data2)
 $varo->add('another_name', \@data2)
 $varo->add({'another_name' => \@data2})

I<Alias>: C<add_data>

Same as L<load|load> except that any previous loads are not L<unload|unload>ed.

=head2 unload

 $varo->unload();

Empties all cached data and calculations upon them, ensuring these will not be used for testing. This will be automatically called with each new load, but, to take care of any development, it could be good practice to call it yourself whenever switching from one dataset for testing to another.

=head2 anova_indep

 $varo->anova_indep()

An implementation of a one-way between-groups analysis of variance. Feeds the class object C<$varo> as follows:

 $varo->{'f_value'}
 $varo->{'df_t'} : the treatment or numerator or between-groups degree(s) of freedom
 $varo->{'df_e'} : the error or denominator or within-groups degree(s) of freedom
 $varo->{'p_value'} : associated with the F-value with the above dfs
 $varo->{'ss_t'} : treatment sums of squares
 $varo->{'ss_e'} : error sums of squares
 $varo->{'ms_t'} : treatment mean squares
 $varo->{'ms_e'} : error mean squares

=head2 anova_dep

 $varo->anova_dep()

I<Alias>: anova_rm

Performs a one-way repeated measures analysis of variance (sphericity assumed). See L<anova_indep|anova_indep> for fed values.

=head2 anova_friedman

 $varo->anova_friedman()

I<Alias>: friedman_test

Performs Friedman's nonparametric analysis of variance - for two or more dependent (matched, related) groups. The statistical attributes now within the class object (see L<anova_indep|anova_indep>) pertain to this test, e.g., $varo->{'chi_value'} gives the chi-square statistic from the Friedman test; and $varo->{'p_value'} gives the associated p-value (area under the right-side, upper tail). There is now no defined 'f_value'.

Accepts, however, one argument: If I<f_equiv> => 1, then, instead of the chi_value, and p_value read off the chi-square distribution, you get the F-value equivalent, with the p-value read off the F-distribution.

See some other module for performing nonparametric pairwise comparisons.

=head2 obrien_test

Performs O'Brien's (1981) test for equality of variances within each group: based on transforming each observation in relation to its group variance and its deviation from its group mean; and performing an ANOVA on these transformed scores (for which the group mean is equal to the variance of the original observations). The procedure is recognised to be robust against violations of normality (unlike F-max).

The statistical attributes now within the class object (see L<anova_indep|anova_indep>) pertain to this test, e.g., $varo->{'f_value'} gives the F-statistic for O'Brien's Test; and $varo->{'p_value'} gives the p-value associated with the F-statistic for O'Brien's Test.

=head2 levene_test

Performs Levene's (1960) test for equality of variances within each group: an ANOVA of the absolute deviations, i.e., absolute value of each observation less its group mean.

The statistical attributes now within the class object (see L<anova_indep|anova_indep>) pertain to this test, e.g., $varo->{'f_value'} gives the F-statistic for Levene's Test; and $varo->{'p_value'} gives the p-value associated with the F-statistic for Levene's Test.

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

 $str = $varo->string(mse => 1, eta_squared => 1, omega_squared => 1, p_precision => integer, s_precision => integer)

Returns a statement of result, in the form of C<F(df_t, df_e) = f_value, p = p_value>; or, for Friedman test C<chi^2(df_t) = chi_value, p = p_value> (to the value of I<p_precision>, if any). Optionally also get MSe, eta_squared and omega_squared values appended to the string, where relevant. These and the test statistic are "sprintf"'d to the I<s_precision> specified (default = 3).

=head2 table

 $tble = $varo->table(p_precision => integer, s_precision => integer);

Returns a table listing the degrees of freedom, sums of squares, and mean squares for the tested "factor" and "error" (between/within groups), and the F and p values. The test statistics are "sprintf"'d to the I<s_precision> specified (default = 3); the p value's precision can be specified by I<p_precision>.

Formatting with right-justification where appropriate is left as an exercise for the user.

=head2 dump

 $varo->dump(title => 'ANOVA test', p_precision => integer, s_precision => integer, mse => 1, eta_squared => 1, omega_squared => 1, verbose => 1)

Prints the string returned by L<string|string>, or, if specified with the attribute I<table> => 1, the table returned by L<table|table>; and the string as well if I<string> => 1. A newline - "\n" - is appended at the end of the print of the string. Above this string or table, a title can also be printed, by giving a value to the optional C<title> attribute.

If I<verbose> => 1, then any curiosities arising in the calculations are noted at the end of other dumps. At the moment, this is only the number of observations that might have been purged were they identified as undefined or not-a-number upon loading/adding.

=head1 REFERENCES

Gardner, R. C. (2001). I<Psychological Statistics using SPSS for Windows>. Upper Saddle River, NJ, US: Prentice Hall. : An interesting source for open-source.

Maxwell, S. E., & Delaney, H. D. (1990). I<Designing Experiments and Analyzing Data: A Model Comparison Perspective.> Belmont, CA, US: Wadsworth.

=head1 SEE ALSO

L<Statistics::FisherPitman|lib::Statistics::FisherPitman> For an alternative to independent groups ANOVA when the variances are unequal.

L<Math::Cephes|lib::Math::Cephes> Probabilities for all F-tests are computed using the C<fdtrc> function in this, rather than the more commonly used L<Statistics::Distributions|lib::Statistics::Distributions> module, as the former appears to be more accurate for higher values of F.

L<Statistics::Descriptive|Statistics::Descriptive> Fundamental calculations of means and variances are left up to this old standard; any limitations/idiosyncrasies therein are naturally passed onto the present one; although the present one purges missing and non-numerical values, unlike Statistics::Descriptive.

L<Statistics::Table::F|Statistics::Table::F> Simply returns an F value. Note that it does not handle missing values, treating them as zero and thus returning an erroneous F-value in these cases.

=head1 BUGS/LIMITATIONS

Computational bugs will hopefully be identified with usage over time.

Optimisation welcomed.

No adjustment for violations of sphericity in repeated measures ANOVA.

Print only of t-test results.

=head1 REVISION HISTORY

=over 4

=item v 0.01

June 2008: Initital release via PAUSE. 

=back

See CHANGES in installation distribution for subsequent updates.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2008 Roderick Garton

rgarton AT cpan DOT org

This program is free software. This module is free software. It may be used, redistributed and/or modified under the stame terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut
