universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1;
arguments = $(process) 
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 50

