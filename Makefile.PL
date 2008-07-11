use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Statistics::ANOVA',
    VERSION_FROM      => 'lib/Statistics/ANOVA.pm', # finds $VERSION
    PREREQ_PM         => {
                         Statistics::TTest => 1.1,
                         Statistics::DependantTTest => 0.03,
                         Algorithm::Combinatorics => 0.23,
                         Statistics::Distributions => 1.02,
                         Statistics::Descriptive => 2.6,
                         Math::Cephes => 0.43,
                         Scalar::Util => 1.18,
                         },
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (
       ABSTRACT_FROM  => 'lib/Statistics/ANOVA.pm', # retrieve abstract from module
       AUTHOR         => 'R. Garton <rgarton AT cpan DOT org>') : ()
      ),
);
