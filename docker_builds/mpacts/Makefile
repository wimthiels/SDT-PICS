default: latest

all: latest py_3.6.8

latest:
	docker build -f Dockerfile \
	--cache-from registry.gitlab.mech.kuleuven.be/mebios-particulate/mpact-docker:latest \
	--build-arg PYTHON_VERSION_TAG=3.7.3 \
	-t mebios-particulate/mpact-docker:latest \
	-t mebios-particulate/mpact-docker:3.7.3 \
	--compress .

py_3.6.8:
	docker build -f Dockerfile \
	--cache-from registry.gitlab.mech.kuleuven.be/mebios-particulate/mpact-docker:3.6.8 \
	--build-arg PYTHON_VERSION_TAG=3.6.8 \
	-t mebios-particulate/mpact-docker:3.6.8 \
	--compress .
