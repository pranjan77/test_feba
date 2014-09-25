#!/usr/bin/perl -w
use strict;
use POSIX qw(sys_wait_h);
my $wait = 2;
my $n = undef;
my $verbose = 0;
my $usage = "Usage: submitter.pl [-shell rsh|ssh] [-host X,Y,Z] [-wait $wait] [-n nCPUs] [-cd dir]"
    . "    [-log name] [-env file] commandsFile\n"
    . "host is the list of machine(s) to run the jobs on (or this machine if none)\n\n"
    . "shell (ignored if no host) is ssh by default\n\n"
    . "cd is the directory to change to before running each job\n\n"
    . "log is the prefix of the filename (relative to the cd dir) to use for each log\n"
    . "(it defaults to commandFile)\n\n"
    . "env should be a list of lines of the form VARNAME=VALUE (lines starting with ! or # ignored)\n\n"
    . "the commands in commandsFile will be executed by either\n"
    . "ssh host '(cd dir; var1=value1...; time LINE) >& dir/log-CMDNO' (if host(s) are specified)\n"
    . "or by (cd dir; var1=value1...; time LINE) >& dir/log-CMDNO (if no host is specified)\n"
    . "(lines starting with ! or # also ignored)\n\n";

die $usage if (scalar @ARGV) % 2 == 0;

my $maxload = (`egrep -c '^processor' /proc/cpuinfo`) + 1;
$n = int(0.6 * $maxload) if !defined $n;
$n = 1 if $n < 1;
my @hosts = ("");
my $cd = "";
my $shell="ssh";
my $logPrefix = "";
my $envFile = "";
my $commandsFile = "";
my @commands = ();
my @envNames = ();
my @envValues = ();
my $debug = 0;

while (scalar @ARGV > 2) {
    my $option = shift @ARGV;
    my $value = shift @ARGV;
    if ($option eq "-host") {
	@hosts = split /,/, $value;
    } elsif ($option eq "-wait") {
	$wait = $value;
    } elsif ($option eq "-n" || $option eq "-load") {
	$n = $value;
    } elsif ($option eq "-cd") {
	$cd = $value;
    } elsif ($option eq "-log") {
	$logPrefix = $value;
    } elsif ($option eq "-env") {
	$envFile = $value;
    } elsif ($option eq "-debug") {
	$debug = $value;
    } elsif ($option eq "-shell") {
	$shell = $value;
	$shell = "rsh -n" if $shell eq "rsh"; # so that we don't pause if in the background
    } elsif ($option eq "-verbose") {
	$verbose = $value;
    } else {
	die "Unrecognized option $option\n\n$usage";
    }
}


$commandsFile = shift @ARGV;
die $usage if scalar @ARGV != 0;
$logPrefix = $commandsFile if ($logPrefix eq "");

open(COMMANDS, "<", $commandsFile) || die "Cannot open commands file $commandsFile";
while(<COMMANDS>) {
    $_ =~ s/[\r\n]+$//;
    if ($_ ne "" && !m/^[!#]/ && !m/^[ \t]+$/) {
       push(@commands, $_);
    }
}
close(COMMANDS) || die "Cannot close commands file $commandsFile";

my $envCommand = "";

if ($envFile ne "") {
    open(ENVFILE, "<", $envFile) || die "Cannot open environment file $envFile";
    while(<ENVFILE>) {
	$_ =~ s/[\r\n]+$//;
	if ($_ ne "" && !m/^[!#]/) {
	    die "Cannot parse environment line $_\n" unless m/^(\S+)=(.*)$/;
            push(@envNames, $1);
	    push(@envValues, $2);
        }
    }
    close(ENVFILE) || die "Cannot close environment file $envFile";
    $envCommand = join("; ", map("$envNames[$_]=$envValues[$_]",0..$#envNames)) . ";";
}


    my $cdPrefix = ($cd eq "") ? "" : "$cd/";
    my $cdCmd = ($cd eq "") ? "" : "cd $cd;";

    if ($verbose) {
	print "HOSTS=",join(",",@hosts),"\n";
	print "Wait=$wait\nJobs=$n\nCD=$cd\nmax_load=$maxload\n";
	print "ENV: $envCommand\n";
	print "CMDS:\n",join("\n",@commands),"\n";
    }
    
    my $ncmd = 0;
    my %jobids = ();
    foreach my $host (@hosts) {
	my @list = map(-1, 0..($n-1));
        $jobids{$host} = \@list;
    }
    my @cmdJobs = map(-1, 0..$#commands);
    my %pidToNCmd = ();

    my $failed = 0;

    while ($ncmd < scalar @commands) {
	foreach my $host (@hosts) {
	    my @joblist = @{ $jobids{$host} };
	    my $shellString = $host eq "" ? "" : "$shell $host ";
	    foreach my $slot (0..($n-1)) {
		# See if job in this slot is done yet
		if ($joblist[$slot] != -1) {
		    my $child = waitpid($joblist[$slot], WNOHANG);
		    if ($child == $joblist[$slot]) {
			my $ncmdOld = $pidToNCmd{$child};
			if ($? != 0) {
			    print "Failed child $child job $ncmdOld status $? : $commands[$ncmdOld]\n";
			    $failed++;
			} else {
			    print "Completed child $child job $ncmdOld : $commands[$ncmdOld]\n" if $verbose > 0;
			}
			$joblist[$slot] = -1;
		    }
		}
		if ($joblist[$slot] == -1 && $ncmd < scalar @commands) {
		    my $upstring = `$shellString uptime`;
		    if (!$upstring) {
			die "Cannot run uptime on $host\n" if (scalar @hosts < 2);
			#else
			print "Cannot run uptime on $host\n";
		    }
		    my $up;
		    if ($upstring) {
			my @uparray = split(' ',$upstring);
			$up = $uparray[$#uparray-2];
			$up =~ s/,//;
			print "host $host slot $slot up $up ncmd $ncmd\n" if $verbose;
		    }
		    if ($upstring && $up < $maxload) {
			my $logstart = "$cdPrefix$logPrefix";
			$logstart = $logPrefix if $cdPrefix =~ m!/$! && $logPrefix =~ m!^/!;
			my $cmd = "($cdCmd $envCommand time $commands[$ncmd]) >& $logstart-$ncmd.log";
			my $execute = $host eq "" ? $cmd : "$shellString '$cmd'";
			if ($debug) {
			    print "Would execute $execute\n";
			} else {
			    my $child = fork();
			    if ($child == 0) {
				# we are the child;
				system($execute);
				print "job $ncmd FAILED $? : $commands[$ncmd]\n" if $? != 0;
				exit(0);
			    } else {
				$joblist[$slot] = $child;
				$cmdJobs[$ncmd] = $child;
				$pidToNCmd{$child} = $ncmd;
				print "Spawned $child for job $ncmd : $execute\n" if $verbose > 0;
			    }
			}
			$ncmd++;
		    } else {
			sleep($wait);
		    }
		}
	    }
	    $jobids{$host} = \@joblist;
	}
    }

    # And now wait for all the children to finish
    foreach my $host (@hosts) {
	my @joblist = @{ $jobids{$host} };
	foreach my $slot (0..($n-1)) {
	    if ($joblist[$slot] != -1) {
		my $ncmdOld = $pidToNCmd{$joblist[$slot]};
		my $child = waitpid($joblist[$slot], 0);
		if ($child ne $joblist[$slot]) {
		    print "Ack!!! waitpid for child $joblist[$slot] failed\n";
		} elsif ($? != 0) {
		    print "Failed child $child job $ncmdOld status $? : $commands[$ncmdOld]\n";
		    $failed++;
		} else {
		    print "Completed child $child job $ncmdOld : $commands[$ncmdOld]\n" if $verbose > 0;
		}
	    }
	}
    }
    if (!$debug) {
	if ($failed == 0) {
	    print "Done: all " . scalar(@commands) . " jobs completed successfully\n";
	    
	} else {
	    print "Done: $failed of " . scalar(@commands) . " jobs failed\n";
	}
	sleep(1);
    }
    exit($failed == 0 ? 0 : 1);
