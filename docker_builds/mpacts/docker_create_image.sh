#!/bin/bash

DOCKER_IMAGE_NAME="mpacts"
DOCKER_IMAGE_FOLDER="."  #location of Dockerfile (use "." if you are in the folder)

source ../../testinfra_config
ROOT_REPO="${PWD/docker_builds\/mpacts/}"
DATA_FOLDER="${ROOT_REPO}${DATA_FOLDER/.\//}"
SCRIPT_FOLDER="${ROOT_REPO}${SCRIPT_FOLDER/.\//}"
PARAM_FOLDER="${ROOT_REPO}${PARAM_FOLDER/.\//}"
printf "DATA_FOLDER=$DATA_FOLDER  \n SCRIPT_FOLDER=$SCRIPT_FOLDER \n PARAM_FOLDER=$PARAM_FOLDER --container_name=$MPACTS_CONTAINER \n"

DOCKER_CONTAINER_NAME=$MPACTS_CONTAINER



#step 1 : BUILD IMAGE (put yourself in folder of dockerimage)
docker build -t $DOCKER_IMAGE_NAME $DOCKER_IMAGE_FOLDER


#remove container if exist
if [ "$(docker ps -aq  -f name=${DOCKER_CONTAINER_NAME})" ]; then
	printf "$DOCKER_CONTAINER_NAME is already present.  Container will be removed: \n"; 
	docker stop $DOCKER_CONTAINER_NAME
	docker rm $DOCKER_CONTAINER_NAME
fi

#step2  : create container using 3 shared directories (=bind mounts, so accesible by both docker and host) 
docker run -di \
	--name $DOCKER_CONTAINER_NAME \
	-v "$SCRIPT_FOLDER":/home/docker/SCRIPT \
	-v "$DATA_FOLDER":/home/docker/DATA \
	-v "$PARAM_FOLDER":/home/docker/TESTINFRA \
	$DOCKER_IMAGE_NAME


