# Copyright (c) 2019, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# This material may be reproduced by or for the U.S. Government 
# pursuant to the copyright license under the clauses at DFARS 
# 252.227-7013 and 252.227-7014.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the copyright holder nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.


use strict;

open(CONFIG, "config.h") || die "cannot open config.h for reading\n";
my @config = <CONFIG>;
close(CONFIG);
my $enabled = grep /ENABLE_NONBLOCKING 1/, @config;
my $old_enabled = $enabled;

my @lines = <>;
my @blines = grep /^1/, @lines;
my @nblines = grep /^0/, @lines;
die "incomplete data for tuning\n" if $#blines != $#nblines;
die "insufficient data for tuning\n" if $#blines < 10;

my @bdata = map { (split)[2] } @blines;
my @nbdata = map { (split)[2] } @nblines;

if (! $enabled) {
    for (my $i = 0; $i < @bdata; $i++) {
        $bdata[$i] = ($bdata[$i] + $nbdata[$i]) / 2;
    }
}
my ($peak, $value) = &findPeak(\@bdata);
if ($enabled) {
    my ($npeak, $nvalue) = &findPeak(\@nbdata);
    if ($nvalue >= 0.95 * $value) {
        $peak = $npeak;
        $value = $nvalue;
        @blines = @nblines;
    }
    else {
        $enabled = 0;
    }
}

my $bufsiz = (split ' ',$blines[$peak])[1];
my $old_bufsiz = (join '', grep /CONVEY_BUFFER_SIZE [0-9]+\s*$/, @config);
$old_bufsiz =~ s/^.*CONVEY_BUFFER_SIZE ([0-9]+)\s$/\1/;

printf "Estimated bandwidth is %.1f MB/sec/PE\n", ($value);
print "Setting ENABLE_NONBLOCKING to $enabled (was $old_enabled)\n";
print "Setting CONVEY_BUFFER_SIZE to $bufsiz (was $old_bufsiz)\n";
print "Rerun 'make' to build these values into the library.\n";

my $config = join '', @config;
$config =~ s/ENABLE_NONBLOCKING [0-1]/ENABLE_NONBLOCKING $enabled/gm;
$config =~ s/CONVEY_BUFFER_SIZE [0-9]+\s*$/CONVEY_BUFFER_SIZE $bufsiz/gm;
open(CONFIG, ">config.h") || die "cannot open config.h for writing\n";
print CONFIG $config;
close(CONFIG);

exit;


sub findPeak {
    my ($array) = @_;
    my $peak = 0;
    my $value = 0;
    for (my $i = 1; $i < @$array - 1; $i++) {
        my $avg = ($$array[$i-1] + $$array[$i] + $$array[$i+1]) / 3;
        if ($avg > $value) {
            $value = $avg;
            $peak = $i;
        }
    }
    return ($peak, $value);
}
