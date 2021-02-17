#!/bin/bash
# testinfra_execution layout = runIC runID groupID suffix_statemod  (eg. : EXE-;C001a;group01;_test)
# 1) run : if EXEC than testcase will be run with ./testinfra.sh command
# 2) runID : testcase label
# 3) groupID : optional, to be used with -g 
# 4) suffix_statemod : optional , will be concatenated to INPUT_CSV (testinfra.xls) (if you want to redirect to another statemodifier file)
# useage
# ./testinfra.sh                  => all EXEC testcases will be run
# ./testinfra.sh restart 5        => all EXEC testcases will be run, all restarted from fifth module
# ./testinfra.sh C001f            => C001f testcase will be run
# ./testinfra.sh C001f 5          => C001f testcase will be run, restart from fifth module
# ./testinfra.sh -g paperTRA      => all testcases with grouplabel "paperTRA" will be run
# ./testinfra.sh -g -b paperTRA 5 => all testcases with grouplabel "paperTRA" will be run, all restarted from fifth module, with no debugging

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
	b)  batch=true  #in batch docker is called non-interactively, debugging is bypassed
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

#source ~/anaconda3/etc/profile.d/conda.sh #get access to conda functions (no longer needed)
echo "$HOSTNAME"
source testinfra_config  # set paths and docker setup for the current environment
SCRIPT_FOLDER="${SCRIPT_FOLDER/./$PWD}" #resolve relative paths
DATA_FOLDER="${DATA_FOLDER/./$PWD}" #resolve relative paths

if [[ "$batch" == true ]]; then
	sed -i 's/#*disable/disable/' .pdbrc
fi


echo "testinfra.sh parms => verbose=$verbose, testcase: $1, Leftover: $2, flag_use_exec=$flag_use_exec"
echo "script folder = ${SCRIPT_FOLDER} // data folder = ${DATA_FOLDER} // param folder = ${PARAM_FOLDER} "

#Set derived locations----------------------------------------------------------------------------------------------------------------------------------------
declare -r INFILE="${PARAM_FOLDER}/testinfra_execution.txt"  #example record EXEC;M10001;group01;_paper_data;
declare -r XML_OUTPUT_FOLDER="${DATA_FOLDER}/XML"  #this also points to the outputfolder
declare -r XML_INPUT="${PARAM_FOLDER}/PARAMS.xml"  
declare -r BACKUP_FOLDER="${DATA_FOLDER}/BACKUP"
  # docker equivalents
declare -r ROOT_DOCKER="/home/docker"
declare -r SCRIPT_FOLDER_DOCKER="${ROOT_DOCKER}/SCRIPT"  # 
declare -r DATA_FOLDER_DOCKER="${ROOT_DOCKER}/DATA"  # 
#declare -r PARAM_FOLDER_DOCKER=$SCRIPT_FOLDER_DOCKER # #moved to config
declare -r XML_OUTPUT_FOLDER_DOCKER="DATA/XML"  
declare -r XML_INPUT_DOCKER="${PARAM_FOLDER_DOCKER}/PARAMS.xml"

a_mpacts_container=(sphere_meshing.py mpacts_PiCS.py mpacts_PiCS_presim.py mpacts_PiCS_LUT.py gen_artificial.py)
# start docker containers
if [[ $MPACTS_CONTAINER != "" ]] ; then
	if [ $(docker inspect -f '{{.State.Running}}' ${MPACTS_CONTAINER}) = "false" ]; then
		printf "Starting $MPACTS_CONTAINER "
		docker start ${MPACTS_CONTAINER}; 
	fi
	printf "$MPACTS_CONTAINER  is running : "
	containerID= docker ps | awk '$2 ~ /mpacts/' | awk 'NR==1{print $1}'
	
fi

if [[ $SDT_CONTAINER != "" ]] ; then
	if [ $(docker inspect -f '{{.State.Running}}' ${SDT_CONTAINER}) = "false" ]; then
		printf "Starting $SDT_CONTAINER "
		docker start ${SDT_CONTAINER}; 
	fi
	printf "$SDT_CONTAINER  is running : "
	containerID= docker ps | awk '$2 ~ /sdt/' | awk 'NR==1{print $1}'
fi


#EXECUTE TESTCASES========================================================================================================================================

filelines=`cat $INFILE`
for line in $filelines ; do				#looping over testcases in testinfra_execution.txt
	IFS=';' read runIC runID groupID suffix_statemod<<<$line  #parse line
	if  [[ ($flag_use_exec == true   &&   $runIC == "EXEC" ) ||  ( $flag_use_exec == false   &&  $runID == $arg_runID ) ||  ( $grouping == true   &&  ${groupID} == ${arg_runID}* ) ]]
	then
		echo "========================================================= executing $runID ======================================================================================"
		echo "run parameters parsed : " $runIC $runID
		export IFS=","  

		#STATE_SETTER : modify the param file view script
		declare INPUT_CSV="${PARAM_FOLDER}/testinfra${suffix_statemod}.xlsx"
		declare INPUT_CSV_DOCKER="${PARAM_FOLDER_DOCKER}/testinfra${suffix_statemod}.xlsx"


		if [[ "$SDT_CONTAINER" != "" ]] ; then
			printf "start docker run for xml_modifier\n"
			docker exec -t ${SDT_CONTAINER}  python "${SCRIPT_FOLDER_DOCKER}/xml_modifier.py"  $INPUT_CSV_DOCKER $XML_OUTPUT_FOLDER_DOCKER $runID $XML_INPUT_DOCKER $DATA_FOLDER_DOCKER $SCRIPT_FOLDER_DOCKER  $PARAM_FOLDER_DOCKER $restart_module
		else
			python "${SCRIPT_FOLDER}/xml_modifier.py"  $INPUT_CSV $XML_OUTPUT_FOLDER $runID $XML_INPUT $DATA_FOLDER $SCRIPT_FOLDER  $PARAM_FOLDER $restart_module
		fi

		#EXTRACT RUN VARIABLES and OUTPUT FOLDER
		modules=$(grep -o -P -m 1 '(?<=modules value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
		modules_snapshot=$(grep -o -P '(?<=modules_snapshot value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
		if [ "${modules_snapshot}" == "None" ] ; then
			modules_snapshot="${modules},testinfra.sh"
		else
			modules_snapshot="${modules},testinfra.sh,${modules_snapshot}"  # concatenate
		fi

		OUTPUT_FOLDER=$(grep -o -P '(?<=OUTPUT_FOLDER_ROOT value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
		OUTPUT_FOLDER_ROOT_SUB=$(grep -o -P '(?<=OUTPUT_FOLDER_ROOT_SUB value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
		if [ "$OUTPUT_FOLDER_ROOT_SUB" != "" ] ; then
			echo "OUTPUT_FOLDER_ROOT_SUB is added = ${OUTPUT_FOLDER_ROOT_SUB}"
			OUTPUT_FOLDER="${OUTPUT_FOLDER}/${OUTPUT_FOLDER_ROOT_SUB}"
		fi
		OUTPUT_FOLDER="${OUTPUT_FOLDER/$DATA_FOLDER_DOCKER/$DATA_FOLDER}"
		echo $OUTPUT_FOLDER $DATA_FOLDER_DOCKER $DATA_FOLDER
		echo "Testinfra.sh: OUTPUT_FOLDER = ${OUTPUT_FOLDER}"

		backupIC=$(grep -o -P '(?<=backupIC value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")

		#take a snapshot of selected pythonmodules
		echo "  >taking a snapshot of modules<   " $modules_snapshot
		for module in $modules_snapshot; do
			#echo "snapshot module=" $module  "from-->" "${SCRIPT_FOLDER}/${module}" "TO->>"  "${OUTPUT_FOLDER}/_LOG/modules_snapshot/"
			rsync  -r "${SCRIPT_FOLDER}/${module}"  "${OUTPUT_FOLDER}/_LOG/modules_snapshot/"
		done
		rsync  $INPUT_CSV  "${OUTPUT_FOLDER}/_LOG"
		mv -f "${OUTPUT_FOLDER}/_LOG/testinfra${suffix_statemod}.xlsx" "${OUTPUT_FOLDER}/_LOG/testinfra_${runID}.xlsx"  #rename
		rsync  "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml"  "${OUTPUT_FOLDER}/_LOG"  
			
		#EXECUTING
		ix_module=0
		for module in $modules; do            				#looping over executing steps
			
			ix_module=$(( ix_module + 1 ))
			ic_module_executed=false

			if [ $ix_module -lt $restart_module	 ] ; then
				continue
			fi

			echo "  >$module <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

			# R run
			if [[ "$module" == *".R" ]]; then
				Rscript "${SCRIPT_FOLDER}/${module}" 2>&1 | tee "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"
				ic_module_executed=true
			fi

			if [[ "$module" == "sphere_meshing_CGAL.py" ]]; then
				input_file_spheres=$(grep -o -P '(?<=input_file_spheres value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
				output_folder_meshes=$(grep -o -P '(?<=output_folder_meshes value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
				shrink_factor=$(grep -o -P '(?<=shrink_factor value=\").+(?=\")' "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml")
				mkdir -p "${OUTPUT_FOLDER}/010_sphere_meshing"  #temp
				#echo "  -  starting binary sphere_meshing_CGAL ; ${SCRIPT_FOLDER}/sphere_meshing_CGAL ${input_file_spheres} ${output_folder_meshes} ${shrink_factor}"
				if [[ "$SDT_CONTAINER" != "" ]] ; then
					printf "start docker run for C binary sphere_meshing_CGAL"
					docker exec -t ${SDT_CONTAINER}  "${SCRIPT_FOLDER_DOCKER}/sphere_meshing_CGAL" ${input_file_spheres} ${output_folder_meshes} ${shrink_factor}
				else
					"${SCRIPT_FOLDER}/sphere_meshing_CGAL" ${input_file_spheres} ${output_folder_meshes} ${shrink_factor} 2>&1 | tee "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"
					
				fi
				
			fi

			# Mpacts run (always in container)
			for i in "${a_mpacts_container[@]}"; do       # run under docker image mpacts_container
				if [ "$i" == $module ] ; then

					if [[ "$ic_dockerenv" == "true" ]] ; then
						python "${SCRIPT_FOLDER}/${module}" "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml"  #already in docker env 
					else
						printf "start docker run for $module under $MPACTS_CONTAINER \n"
						if  [[ "$batch" == true ]] ; then
							docker exec -t ${MPACTS_CONTAINER}  python "${SCRIPT_FOLDER_DOCKER}/${module}" "${XML_OUTPUT_FOLDER_DOCKER}/testinfra_${runID}.xml"  | tee "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"
						else
							docker exec -it ${MPACTS_CONTAINER}  python "${SCRIPT_FOLDER_DOCKER}/${module}" "${XML_OUTPUT_FOLDER_DOCKER}/testinfra_${runID}.xml"  | tee "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"
						fi
					fi
					ic_module_executed=true
				fi
			done
			
			# python3 run (=default)
			if [ "$ic_module_executed" != true ] ; then

				if [[ "$SDT_CONTAINER" != "" ]] ; then
					printf "start docker run for $module under $SDT_CONTAINER \n"
					docker exec -t ${SDT_CONTAINER}  python "${SCRIPT_FOLDER_DOCKER}/${module}" "${XML_OUTPUT_FOLDER_DOCKER}/testinfra_${runID}.xml"  | tee "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"
				else
					python  "${SCRIPT_FOLDER}/${module}" "${XML_OUTPUT_FOLDER}/testinfra_${runID}.xml" 2>&1 | tee "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"
					
				fi

				ic_module_executed=true
			fi

			rsync  --append "${DATA_FOLDER}/consoleoutput/consoleOutput_${runID}.txt"  "${OUTPUT_FOLDER}/_LOG"    #copy consoleoutput for debugging

		
		done
		echo "========================================================= end $runID ======================================================================================"

		#backup files (optional)
		if [[ "$backupIC" == "TRUE" ]]; then
			echo "${OUTPUT_FOLDER} -> ${BACKUP_FOLDER}"
			rsync -r "${OUTPUT_FOLDER}" $BACKUP_FOLDER
		fi
		
	else
		printf "_"
	fi
	
done 

end=$(date +%s)

timestamp() {
  date +"%T" # current time
}
echo "timestamp: $(timestamp)"
awk -v t=$SECONDS 'BEGIN{t=int(t*1000); printf "Elapsed Time (HH:MM:SS): %d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}'

if [[ "$ic_dockerenv" != "true" ]] ; then
	$SHELL  #keep shell open
fi

exit 0