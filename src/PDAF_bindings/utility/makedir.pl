#!/usr/bin/perl -w
# Usage: ./makedir.pl total_ens_number
# example: ./makedir.pl 8
#
# Create dirs etc for ensemble runs
if (@ARGV != 1)
{ print "USAGE: $0 <n_tasks>\n"; exit(1); }
$ntasks=$ARGV[0];

for($i=1;$i<=$ntasks; $i++)
{
  $j=sprintf('%03d',$i);
  system "rm -rf schism_$j; mkdir schism_$j";
  system "cd schism_$j; ln -sf ../* .; mkdir outputs; rm -f schism* run_* *.sh *.pl";
#  open(RE,">schism_$j.cfg");
#  print RE "simulationDirectory: ens$j\n";
#  close(RE);
}#for i
