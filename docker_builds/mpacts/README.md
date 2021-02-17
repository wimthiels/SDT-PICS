# mpacts-docker

Docker file to create a ready to be used linux machine for Mpacts users. Mpacts is installed via debian package with python 3.6.7 version. 

Usage:

Step1.	Create a docker image from Dockerfile 

		For example run 'docker build -t mpacts .' command from the root directory. This will create a docker image with name 'mpacts'

Step2.	Create container from docker image
		docker run -di --name Container_Name -v "Path_to_shared_directory":/home/docker/shared Docker_Image_Name
		
		For example 'docker run -di --name mpacts_container -v "C:/mpact-docker/shared":/home/docker/shared mpacts'. This will create a docker container with name 'mpacts_container', mapped to directory 'shared'.

Step3.	Access container
		docker exec -it Name_of_Container bash

		For example 'docker exec -it mpacts_container bash, will provide shell access to mpacts_container.

