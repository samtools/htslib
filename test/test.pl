#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
#use IPC::Open2;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

my $opts = parse_params();

test_tabix($opts,in=>'merge.a',reg=>'2:3199812-3199812',out=>'tabix.2.3199812.out');
test_vcf_check($opts,in=>'check',out=>'check.chk');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.out',args=>'-n =2');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.both.out',args=>'-n =2 -c both');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.any.out',args=>'-n =2 -c any');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.C.out',args=>'-C -c any');
test_vcf_isec2($opts,vcf_in=>['isec.a'],tab_in=>'isec',out=>'isec.tab.out',args=>'');
test_vcf_merge($opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.out');

print "\nNumber of tests:\n";
printf "    total   .. %d\n", $$opts{nok}+$$opts{nfailed};
printf "    passed  .. %d\n", $$opts{nok};
printf "    failed  .. %d\n", $$opts{nfailed};
print "\n";

exit;

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "About: htslib consistency test script\n",
        "Usage: test.pl [OPTIONS]\n",
        "Options:\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = { keep_files=>0, nok=>0, nfailed=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            't|temp-dir:s' => \$$opts{keep_files}, 
            'r|redo-outputs' => \$$opts{redo_outputs}, 
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp} = $$opts{keep_files} ? $$opts{keep_files} : tempdir(CLEANUP=>1);
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    return $opts;
}
sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid) 
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    } 
    else 
    {      
        # child
        exec('/bin/bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [/bin/sh -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed: $cmd\n", $out); }
    return $out;
}
sub test_cmd
{
    my ($opts,%args) = @_;
    if ( !exists($args{out}) )
    {
        if ( !exists($args{in}) ) { error("FIXME: expected out or in key\n"); }
        $args{out} = "$args{in}.out";
    }
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n";
    print "\t$args{cmd}\n";

    my ($ret,$out) = _cmd($args{cmd});
    if ( $ret ) { failed($opts,$test); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
        print $fh $out;
        close($fh);
        my ($ret,$out) = _cmd("diff -q $$opts{path}/$args{out} $$opts{path}/$args{out}.old");
        if ( !$ret && $out eq '' ) { unlink("$$opts{path}/$args{out}.old"); }
        else
        {
            print "\tthe expected output changed, saving:\n";
            print "\t  old .. $$opts{path}/$args{out}.old\n";
            print "\t  new .. $$opts{path}/$args{out}\n";
        }
    }
    my $exp = '';
    if ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
        my @exp = <$fh>;
        $exp = join('',@exp);
        close($fh);
    }
    elsif ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}: $!"); return; }

    if ( $exp ne $out ) 
    { 
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
        print $fh $out;
        close($fh);
        failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new"); 
        return; 
    }
    passed($opts,$test);
}
sub failed
{
    my ($opts,$test,$reason) = @_;
    $$opts{nfailed}++;
    if ( defined $reason ) { print "\n\t$reason"; }
    print "\n.. failed ...\n\n";
}
sub passed
{
    my ($opts,$test) = @_;
    $$opts{nok}++;
    print ".. ok\n\n";
}
sub is_file_newer
{
    my ($afile,$bfile) = @_;
    my (@astat) = stat($afile) or return 0;
    my (@bstat) = stat($bfile) or return 0;
    if ( $astat[9]>$bstat[9] ) { return 1 }
    return 0;
}
sub bgzip_tabix
{
    my ($opts,%args) = @_;
    my $file = "$args{file}.$args{suffix}";
    if ( $$opts{redo_outputs} or !-e "$$opts{tmp}/$file.gz" or is_file_newer("$$opts{path}/$file","$$opts{tmp}/$file.gz") )
    {
        cmd("cat $$opts{path}/$file | bgzip -c > $$opts{tmp}/$file.gz");
    }
    if ( $$opts{redo_outputs} or !-e "$$opts{tmp}/$file.gz.tbi" or is_file_newer("$$opts{tmp}/$file.gz","$$opts{tmp}/$file.gz.tbi") )
    {
        cmd("$$opts{bin}/htscmd tabix -f $args{args} $$opts{tmp}/$file.gz");
    }
}
sub bgzip_tabix_vcf
{
    my ($opts,$file) = @_;
    bgzip_tabix($opts,file=>$file,suffix=>'vcf',args=>'-p vcf');
}


# The tests --------------------------

sub test_tabix
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/htscmd tabix $$opts{tmp}/$args{in}.vcf.gz $args{reg}");
}
sub test_vcf_check
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/htscmd vcfcheck -s - $$opts{tmp}/$args{in}.vcf.gz | grep -v '^# The command' | grep -v '^ID\t'");
}
sub test_vcf_merge
{
    my ($opts,%args) = @_;
    my @files;
    for my $file (@{$args{in}})
    {
        bgzip_tabix_vcf($opts,$file);
        push @files, "$$opts{tmp}/$file.vcf.gz";
    }
    my $files = join(' ',@files);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/htscmd vcfmerge $files");
}
sub test_vcf_isec
{
    my ($opts,%args) = @_;
    my @files;
    for my $file (@{$args{in}})
    {
        bgzip_tabix_vcf($opts,$file);
        push @files, "$$opts{tmp}/$file.vcf.gz";
    }
    my $files = join(' ',@files);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/htscmd vcfisec $args{args} $files");
}
sub test_vcf_isec2
{
    my ($opts,%args) = @_;
    my @files;
    for my $file (@{$args{vcf_in}})
    {
        bgzip_tabix_vcf($opts,$file);
        push @files, "$$opts{tmp}/$file.vcf.gz";
    }
    my $files = join(' ',@files);
    bgzip_tabix($opts,file=>$args{tab_in},suffix=>'tab',args=>'-s 1 -b 2 -e 3');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/htscmd vcfisec $args{args} -s $$opts{tmp}/$args{tab_in}.tab.gz $files");
}

