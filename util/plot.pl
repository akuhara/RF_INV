#!/usr/bin/perl
use strict;
use warnings;
use SAC;
`gmt gmtset FONT_ANNOT_PRIMARY 10p,Helvetica FONT_ANNOT_SECONDARY 10p,Helvetica FONT_LABEL 11p,Helvetica FONT_TITLE 12p,Helvetica`;

# Get arguments
my $param_file = "./params.in";
if (defined $ARGV[0]) {
    $param_file = $ARGV[0];
}
else {
    print STDERR "\n";
    print STDERR "USAGE: perl plot.pl [param.in] (test)\n";
    print STDERR "\n";
    print STDERR " NOTE: This program is assumed to be called at" .
	" working directory of inversion\n";
    print STDERR "\n";
    exit;
}
my $test_mode = 0;
if (defined $ARGV[1] and $ARGV[1] eq "test") {
    $test_mode = 1;
}

# Read parameter file
my %param = %{get_param($param_file)};

# Input file names
my $out_dir = $param{out_dir};
my $nk_file    = "$out_dir/num_interface.ppd";
my $sig_file   = "$out_dir/sigma.ppd";
my $syn_file   = "$out_dir/syn_trace.ppd";
my $vs_file    = "$out_dir/vs_z.ppd";
my $vp_file    = "$out_dir/vp_z.ppd";




# Layout
#
my $width_nk = "7c";
my $hight_nk = "4c";
my $ymax_nk  = 0.3;
my $ymin_nk  = 0.0;
my $xtic_nk  = "5f1";
my $ytic_nk  = "0.1";
my $xlabel_nk = "# of layer interface";
my $ylabel_nk = "Probability";
my $xshift_nk = "2c";
my $yshift_nk = "24c";
#
my $width_sig = "7c";
my $hight_sig = "4c";
my $ymax_sig  = 0.2;
my $ymin_sig  = 0.0;
my $xtic_sig  = "0.01";
my $ytic_sig  = "0.05";
my $xlabel_sig = "Noise \@~s\@~";;
my $ylabel_sig = "Probability";
my $xshift_sig = "10c";
my $yshift_sig = "0c";
# Synthetic 
my $width_syn = "10c";
my $hight_syn = "4c";
my $ymin_syn  = -0.4;
my $ymax_syn  = 0.4;
my $zmin_syn  = 0.0;
my $zmax_syn  = 0.2;
my $xtic_syn  = "2f1";
my $ytic_syn  = "0.5f0.1";
my $xlabel_syn = "Time after P (s)";
my $ylabel_syn = "RF amp.";
my $xshift_syn = "-$xshift_sig";
my $yshift_syn = "-7c";
# Vs profile
my $width_vs = "5c";
my $hight_vs = "-9c";
my $ymin_vs  = -1.2;
my $ymax_vs  = 1.2;
my $zmin_vs  = 0.0;
my $zmax_vs  = 0.1;
my $xtic_vs  = "1f0.2";
my $ytic_vs  = "10f5";
my $xlabel_vs = "Vs (km/s)";
my $ylabel_vs = "Depth (km)";
my $xshift_vs = "0c";
my $yshift_vs = "-12c";  
# Vp
my $width_vp = "5c";
my $hight_vp = "-9c";
my $ymin_vp  = -1.2;
my $ymax_vp  = 1.2;
my $zmin_vp  = 0.0;
my $zmax_vp  = 0.1;
my $xtic_vp  = "2f1";
my $ytic_vp  = "10f5";
my $ztic_vp  = 0.1;
my $xlabel_vp = "Vp (km/s)";
my $ylabel_vp = "Depth (km)";
my $xshift_vp = "7c";
my $yshift_vp = "0c";  


#---------------------------------------------------------------------

foreach my $itrc (1..$param{ntrc}) {

    # output file name
    my $out = sprintf "$out_dir/plot.%02d.ps", $itrc;

    # Number of layer interfaces
    system "gmt psxy $nk_file -Sb1u -W1.0 -Ggray "
	. "-JX$width_nk/$hight_nk "
	. "-R$param{k_min}/$param{k_max}/$ymin_nk/$ymax_nk "
	. "-B$xtic_nk:\"$xlabel_nk\":/${ytic_nk}:\"$ylabel_nk\":WSne "
	. "-X$xshift_nk -Y$yshift_nk "
	. "-K -P > $out";
    
    # Noise sigma
    my $dbin_sig = ($param{sig_max} - $param{sig_min}) / $param{nbin_sig};
    system "gawk '\$3=='$itrc'{print \$1, \$2}' $sig_file "
	. "| gmt psxy -Sb${dbin_sig}u -W1.0 -Ggray "
	. "-JX$width_sig/$hight_sig "
	. "-R$param{sig_min}/$param{sig_max}/$ymin_sig/$ymax_sig "
	. "-B$xtic_sig:\"$xlabel_sig\":/${ytic_sig}:\"$ylabel_sig\":WSne "
	. "-X$xshift_sig -Y$yshift_sig "
	. "-O -K -P >> $out";

    # Synthetic trace
    my $dbin_amp = ($param{amp_max} - $param{amp_min}) / $param{nbin_amp};
    my $delta = get_header($param{obs_files}->[$itrc-1], 'delta');
    my $dz_syn = ($zmax_syn - $zmin_syn) / 10.0;
    $delta = sprintf "%.5f", $delta;
    $dbin_amp = sprintf "%.5f", $dbin_amp;
    my $amp_min = $param{amp_min} + 0.5 * $dbin_amp;
    my $amp_max = $param{amp_max} - 0.5 * $dbin_amp;
    system "gawk '\$4=='$itrc'{print \$1, \$2, \$3}' $syn_file "
	. "| gmt xyz2grd -G/tmp/syn.grd -I$delta/$dbin_amp "
	. "-R$param{t_start}/$param{t_end}/$amp_min/$amp_max";
    system "gmt makecpt -Cocean -I -Z -T$zmin_syn/$zmax_syn/$dz_syn -Do > /tmp/syn.cpt";
    system "gmt grdimage /tmp/syn.grd -C/tmp/syn.cpt "
	. "-JX$width_syn/$hight_syn "
	. "-R$param{t_start}/$param{t_end}/$ymin_syn/$ymax_syn "
	. "-B$xtic_syn:\"$xlabel_syn\":/${ytic_syn}:\"$ylabel_syn\":WSne "
	. "-X$xshift_syn -Y$yshift_syn "
	. "-O -K -P >> $out";
    my $t_beg = get_header($param{obs_files}->[$itrc-1], 'b');
    my $it1 = int($param{t_start} - $t_beg + 0.5) / $delta;
    my $it2 = int($param{t_end} - $t_beg + 0.5) / $delta;
    my @obs = @{get_data($param{obs_files}->[$itrc-1])}[$it1..$it2];
    open my $OBS, "| gmt psxy -W1.0,red -J -R -O -K -P >> $out" or die;
    foreach my $i (0..$#obs) {
	my $t = $i * $delta + $param{t_start};
	print {$OBS} "$t $obs[$i]\n";
    }
    close $OBS or die;

    # Vs profile
    my $dbin_z = sprintf "%.5f", $param{z_max} / $param{nbin_z};
    my $dbin_vs = sprintf "%.5f", ($param{vs_max} - $param{vs_min}) / $param{nbin_vs};
    my $dz_vs = ($zmax_vs - $zmin_vs) / 10.0;
    my $vs_min = $param{vs_min} + 0.5 * $dbin_vs;
    my $vs_max = $param{vs_max} - 0.5 * $dbin_vs;
    my $z_min  = 0.5 * $dbin_z;
    my $z_max  = $param{z_max} - 0.5 * $dbin_z;
    system "gawk '{print \$1, \$2, \$3}' $vs_file "
	. "| gmt xyz2grd -G/tmp/vs.grd -I$dbin_vs/$dbin_z "
	. "-R$vs_min/$vs_max/$z_min/$z_max";
    system "gmt makecpt -Cocean -I -Z -T$zmin_vs/$zmax_vs/$dz_vs -Do > /tmp/vs.cpt";
    system "gmt grdimage /tmp/vs.grd -C/tmp/vs.cpt "
	. "-JX$width_vs/$hight_vs "
	. "-R0/$param{vs_max}/0/$param{z_max} "
	. "-B$xtic_vs:\"$xlabel_vs\":/$ytic_vs:\"$ylabel_vs\":WSne "
	. "-X$xshift_vs -Y$yshift_vs "
	. "-O -K -P >> $out";
    open my $V_REF, "<", $param{vel_ref_file} or die;
    my (@z_ref, @vp_ref, @vs_ref);
    my (@vp_c, @vs_c);
    while (my $line = <$V_REF>) {
	chomp $line;
	my @item = split q{ }, $line;
	my ($z_ref, $vp_ref, $vs_ref) = @item;
	my $vs_c = $vs_ref;
	my $vp_c = $vp_ref;
	if ($vs_ref + $param{dvs_min} < $param{vs_min}) {
	    $vs_c = $param{vs_min} - $param{dvs_min} 
	}
	elsif ($vs_ref + $param{dvs_max} > $param{vs_max}) {
	    $vs_c = $param{vs_max} - $param{dvs_max};
	}
	if ($vp_ref + $param{dvp_min} < $param{vp_min}) {
	    $vp_c = $param{vp_min} - $param{dvp_min} 
	}
	elsif ($vp_ref + $param{dvp_max} > $param{vp_max}) {
	    $vp_c = $param{vp_max} - $param{dvp_max};
	}
	push @z_ref, $z_ref;
	push @vp_ref, $vp_ref;
	push @vs_ref, $vs_ref;
	push @vp_c, $vp_c;
	push @vs_c, $vs_c;
    }
    close $V_REF or die;
    open my $VS_REF, "| gmt psxy -W1.0,red -J -R -O -K >> $out" or die;
    foreach my $i (0..$#z_ref) {
	print {$VS_REF} "$vs_ref[$i] $z_ref[$i]\n";
    }
    close $VS_REF or die;
    open my $VS_REF_LOW, "| gmt psxy -W1.0,red,- -J -R -O -K >> $out" or die;
    foreach my $i (0..$#z_ref) {
	print {$VS_REF_LOW} $vs_c[$i] + $param{dvs_min}, q{ }, $z_ref[$i],"\n";
    }
    close $VS_REF_LOW or die;
    open my $VS_REF_UP, "| gmt psxy -W1.0,red,- -J -R -O -K >> $out" or die;
    foreach my $i (0..$#z_ref) {
	print {$VS_REF_UP} $vs_c[$i] + $param{dvs_max}, q{ }, $z_ref[$i],"\n";
    }
    close $VS_REF_UP or die;
    
    if ($test_mode) {
	open my $TEST_VEL, "<", "test_vel" or die;
	open my $TEST_VS, "| gmt psxy -J -R -O -K -W1.0,pink ". 
	    " >> $out" or die;
	my $tmp_dep = 0.0;
	while (my $line = <$TEST_VEL>) {
	    chomp $line;
	    my @item = split q{ }, $line;
	    my ($vs, $h) = @item[1,3];
	    print {$TEST_VS} "$vs $tmp_dep\n";
	    $tmp_dep += $h;
	    print {$TEST_VS} "$vs $tmp_dep\n";
	}
	close $TEST_VEL or die;
	close $TEST_VS or die;
    }


    # Vp profile
    my $dbin_vp = sprintf "%.5f", ($param{vp_max} - $param{vp_min}) / $param{nbin_vp};
    my $dz_vp = ($zmax_vp - $zmin_vp) / 10.0;
    my $vp_min = $param{vp_min} + 0.5 * $dbin_vp;
    my $vp_max = $param{vp_max} - 0.5 * $dbin_vp;
    system "gawk '{print \$1, \$2, \$3}' $vp_file "
	. "| gmt xyz2grd -G/tmp/vp.grd -I$dbin_vp/$dbin_z "
	. "-R$vp_min/$vp_max/$z_min/$z_max";
    system "gmt makecpt -Cocean -I -Z -T$zmin_vp/$zmax_vp/$dz_vp -Do > /tmp/vp.cpt";
    system "gmt grdimage /tmp/vp.grd -C/tmp/vp.cpt "
	. "-JX$width_vp/$hight_vp "
	. "-R0/$param{vp_max}/0/$param{z_max} "
	. "-B$xtic_vp:\"$xlabel_vp\":/$ytic_vp:\"$ylabel_vp\":WSne "
	. "-X$xshift_vp -Y$yshift_vp "
	. "-O -P -K >> $out";
    open my $VP_REF, "| gmt psxy -W1.0,red -J -R -O -K >> $out" or die;
    foreach my $i (0..$#z_ref) {
	print {$VP_REF} "$vp_ref[$i] $z_ref[$i]\n";
    }
    close $VP_REF or die;
    open my $VP_REF_LOW, "| gmt psxy -W1.0,red,- -J -R -O -K >> $out" or die;
    foreach my $i (0..$#z_ref) {
	print {$VP_REF_LOW} $vp_c[$i] + $param{dvp_min}, q{ }, $z_ref[$i],"\n";
    }
    close $VP_REF_LOW or die;
    open my $VP_REF_UP, "| gmt psxy -W1.0,red,- -J -R -O -K >> $out" or die;
    foreach my $i (0..$#z_ref) {
	print {$VP_REF_UP} $vp_c[$i] + $param{dvp_max}, q{ }, $z_ref[$i],"\n";
    }
    close $VP_REF_UP or die;
    if ($test_mode) {
	open my $TEST_VEL, "<", "test_vel" or die;
	open my $TEST_VP, "| gmt psxy -J -R -O -K -W1.0,pink ". 
	    " >> $out" or die;
	my $tmp_dep = 0.0;
	while (my $line = <$TEST_VEL>) {
	    chomp $line;
	    my @item = split q{ }, $line;
	    my ($vp, $h) = @item[0,3];
	    print {$TEST_VP} "$vp $tmp_dep\n";
	    $tmp_dep += $h;
	    print {$TEST_VP} "$vp $tmp_dep\n";
	}
	close $TEST_VEL or die;
	close $TEST_VP or die;
    }
    
    system "gmt psscale -DJCB+ef+o0c/2.0c -C/tmp/vp.cpt " .
	"-B$ztic_vp:\"Probability\": " .
	"-R -J -O >> $out";
}


#---------------------------------------------------------------------

sub get_param {
    my $param_file = $_[0];

    my $i = -1;
    my %param;
    open my $PARAM, "<", $param_file or die;
    while (my $line = <$PARAM>) {
	chomp $line;
	my @item = split q{ }, $line;
	next if ($item[0] =~/^#/);
	if ($i == -1) {
	    $param{out_dir} = $item[0];
	}
	elsif ($i == 0) {
	    $param{nburn} = $item[0];
	}
	elsif ($i == 1) {
	    $param{niter} = $item[0];
	}
	elsif ($i == 2) {
	    $param{ncorr} = $item[0];
	}
	elsif ($i == 3) {
	    $param{nchains} = $item[0];
	}
	elsif ($i == 4) {
	    $param{ncool} = $item[0];
	}
	elsif ($i == 5) {
	    $param{t_high} = $item[0];
	}
	elsif ($i == 6) {
	    $param{iseed} = $item[0];
	}
	elsif ($i == 7) {
	    $param{ntrc} = $item[0];
	}
	elsif ($i == 8) {
	    push @{$param{rayps}}, $item[0];
	    next if (@{$param{rayps}} < $param{ntrc});
	    
	}
	elsif ($i == 9) {
	    push @{$param{gauss}}, $item[0];
	    next if (@{$param{gauss}} < $param{ntrc});
	}
	elsif ($i == 10) {
	    $param{nfft} = $item[0];
	}
	elsif ($i == 11) {
	    push @{$param{obs_files}}, $item[0];
	    next if (@{$param{obs_files}} < $param{ntrc});
	}
	elsif ($i == 12) {
	    ($param{t_start}, $param{t_end}) = @item[0,1];
	}
	elsif ($i == 13) {
	    $param{sdep} = $item[0];
	}
	elsif ($i == 14) {
	    $param{vel_ref_file} = $item[0];
	}
	elsif ($i == 15) {
	    $param{vp_mode} = $item[0];
	}
	elsif ($i == 16) {
	    ($param{k_min}, $param{k_max}) = @item[0,1];
	}
	elsif ($i == 17) {
	    ($param{z_min}, $param{z_max}) = @item[0,1];
	}
	elsif ($i == 18) {
	    ($param{dvs_min}, $param{dvs_max}) = @item[0,1];
	}
	elsif ($i == 19) {
	    ($param{dvp_min}, $param{dvp_max}) = @item[0,1];
	}
	elsif ($i == 20) {
	    ($param{sig_min}, $param{sig_max}) = @item[0,1];
	}
	elsif ($i == 21) {
	    $param{dev_z} = $item[0];
	}
	elsif ($i == 22) {
	    $param{dev_dvs} = $item[0];
	}
	elsif ($i == 23) {
	    $param{dev_dvp} = $item[0];
	}
	elsif ($i == 24) {
	    $param{dev_sig} = $item[0];
	}
	elsif ($i == 25) {
	    $param{nbin_z} = $item[0];
	}
	elsif ($i == 26) {
	    $param{nbin_vs} = $item[0];
	}
	elsif ($i == 27) {
	    $param{nbin_vp} = $item[0];
	}
	elsif ($i == 28) {
	    $param{nbin_sig} = $item[0];
	}
	elsif ($i == 29) {
	    $param{nbin_amp} = $item[0];
	}
	elsif ($i == 30) {
	    ($param{amp_min}, $param{amp_max}) = @item[0,1];
	}
	elsif ($i == 31) {
	    ($param{vp_min}, $param{vp_max}) = @item[0,1];
	}
	elsif ($i == 32) {
	    ($param{vs_min}, $param{vs_max}) = @item[0,1];
	}
	$i ++;
    }
    close $PARAM or die;
    return \%param;
    
}


