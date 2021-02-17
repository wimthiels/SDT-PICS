#!/bin/bash
# will invoke testinfra.sh so that each testcase runs in parallel
# testinfra_parallel_execution layout = runIC runID groupID suffix_statemod  (eg. : EXE-;C001a;group01;_test)
# 1) run : if EXEC than testcase will be run with ./testinfra.sh command
# 2) runID : testcase label
# 3) groupID : optional, to be used with -g 
# 4) suffix_statemod : optional , will be concatenated to INPUT_CSV (testinfra.xls) (if you want to redirect to another statemodifier file)
# useage
# ./testinfra_parallel.sh             => all EXEC testcases will be run
# ./testinfra_parallel.sh restart 5   => all EXEC testcases will be run, all restarted from fifth module
# ./testinfra_parallel.sh C001f       => C001f testcase will be run
# ./testinfra_parallel.sh C001f 5     => C001f testcase will be run, restart from fifth module
# ./testinfra_parallel.sh -g paperTRA => all testcases with grouplabel "paperTRA" will be run
# ./testinfra_parallel.sh -g paperTRA => all testcases with grouplabel "paperTRA" will be run

# naming testcases
#  1) modulenotation = M012b : module test of module 012, run b (C012b = chain starting from module 012
#  2) free notation = D020 or T001
#  3) group notation = 7cell0304 : testcase belonging to group '7cell' third testcase (=replicate), fourth run_nb = state_nb = separate tab in excel
#################################################TESTINFRASTRUCTURE MAIN DRIVER#############################################

start=$(date +%s)

#Parse command line arguments--------------------------------------------------------------------
# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
verbose=false
ic_local="?"
ic_dockerenv=false

while getopts "h?bvg" opt; do
	case "$opt" in
	h?)
		show_help
		exit 0
		;;
	v)  verbose=true
		;;
	g)  grouping=true
		;;
	esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift


#Set variables --------------------------------------------------------------------
declare -x flag_use_exec=false
declare -x arg_runID=""
restart_module=0
if [ $# -eq 0 ]
then
	flag_use_exec=true
else	
	arg_runID=$1

	if [ $arg_runID == "restart" ];then
		flag_use_exec=true
	fi


	if [ $# -eq 2 ]
	then
		restart_module=$2
		echo "restart will take place starting from module= ${restart_module}"
	else	
		restart_module=0
	fi
fi

echo "$HOSTNAME"
if [ $HOSTNAME == "DESKTOP-USJGL6G" ];then
	source ~/anaconda3/etc/profile.d/conda.sh #get access to conda functions
	declare -r DATA_FOLDER="/mnt/f/Downloads/ROET/testinfra"
	declare -r SCRIPT_FOLDER="/mnt/c/Users/wimth/OneDrive/SYNC"
	declare -r PARAM_FOLDER=$SCRIPT_FOLDER
	ic_local=true
	
elif [ $HOSTNAME == "wth-HP-ZBook-Studio-x360-G5" ];then
	source ~/anaconda3/etc/profile.d/conda.sh #get access to conda functions
	declare -r DATA_FOLDER="/home/wth/Downloads/testinfra"
	declare -r SCRIPT_FOLDER="/home/wth/Downloads/SYNC"
	declare -r PARAM_FOLDER=$SCRIPT_FOLDER
	ic_local=true

elif [ -f /.dockerenv ]; then
	echo "info : testinfra.sh run inside a docker container";
	declare -r DATA_FOLDER="/home/docker/DATA"
	declare -r SCRIPT_FOLDER="/home/docker/SCRIPTS"  #= relative to workdir of docker 
	declare -r PARAM_FOLDER="/home/docker/DATA/INPUT"
	ic_local=false
	ic_dockerenv=true

else
	source "${VSC_DATA}/anaconda3/etc/profile.d/conda.sh"
	declare -r DATA_FOLDER="${VSC_DATA}/testinfra"
	declare -r SCRIPT_FOLDER="${VSC_DATA}/SYNC"
	declare -r PARAM_FOLDER=$SCRIPT_FOLDER
	ic_local=false
fi



echo "testinfra_parallel.sh parms => verbose=$verbose, testcase: $1, Leftover: $2, flag_use_exec=$flag_use_exec, ic_local= $ic_local, grouping=$grouping"

#Set locations----------------------------------------------------------------------------------------------------------------------------------------
declare -r INFILE="${PARAM_FOLDER}/testinfra_execution.txt"  #example record EXEC;M10001;group01;_paper_data;

#EXECUTE TESTCASES========================================================================================================================================

filelines=`cat $INFILE`
for line in $filelines ; do				#looping over testcases in testinfra_execution.txt
	IFS=';' read runIC runID groupID suffix_statemod<<<$line  #parse line
	if  [[ ($flag_use_exec == true   &&   $runIC == "EXEC" ) ||  ( $flag_use_exec == false   &&  $runID == $arg_runID ) ||  ( $grouping == true   &&  ${groupID} == ${arg_runID}* ) ]]
	then
		bash testinfra.sh -b $runID $restart_module &>/dev/null &  
		echo "testcase ${runID} is submitted in background.  Restart at module=" $restart_module
	fi
	
done 


if [[ "$ic_dockerenv" != "true" ]] ; then
	$SHELL  #keep shell open
fi