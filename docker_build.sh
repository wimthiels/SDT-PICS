#bin/bash
#Builds both mpacts and sdt docker images and spins a container for each.   (followed by a testrun)

printf "BUILDING DOCKER MPACTS \n"
cd docker_builds/mpacts/
./docker_create_image.sh
cd ../..

printf "BUILDING DOCKER SDT \n"
cd docker_builds/sdt/
./docker_create_image.sh
cd ../..

# printf "LAUNCHING A TESTRUN \n"
# ./sdt-pics.sh pics_test