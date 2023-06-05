universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 2 5 5 10 2 1 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10



universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 2 5 10 5 2 1 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10


universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 2 5 5 10 2 4 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10



universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 2 5 10 5 2 4 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10






universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 1 5 5 10 2 1 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10



universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 1 5 10 5 2 1 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10


universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 1 5 5 10 2 4 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10



universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 1 5 10 5 2 4 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10




universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 3 5 5 10 2 1 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10



universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 3 5 10 5 2 1 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10


universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 3 5 5 10 2 4 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10



universe = vanilla
getenv = true
executable = test.sh
#/opt/matlab/latest/bin/matlab
#inst=$1; eg=$2; distance_gp=$3; ini_size=$4; seq_size =$5; design_algo = $6; True_theta = $7;design_criterion = $8
arguments = $(process) 3 5 10 5 2 4 2
#-nodesktop -nosplash -r "test($(process))"
log = logfile/$(Cluster).log
output = logfile/$(Cluster).$(process).out
error = logfile/$(cluster).$(Process).error
notification = Never
request_memory=2048
request_cpus=1
queue 10