package ALYA::ExecFactory;

#---------------------------------------------------------------------
#--Package that creates the html reports from the xml result----------
#---------------------------------------------------------------------

#usage modules
use strict; use warnings;

    sub instantiate {
        my $classe         = shift;
        my $requested_type = shift;
        my $location       = "ALYA/$requested_type.pm";
        my $class          = "ALYA::$requested_type";

        require $location;

        return $class->new(@_);
    }

1;
