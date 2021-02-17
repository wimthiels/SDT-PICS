FROM ubuntu:18.04

LABEL maintainer="Wim Thiels <wim.thiels@kuleuven.be>"
# based on https://github.com/FNNDSC/ubuntu-python3/blob/master/Dockerfile

USER root
WORKDIR /root

SHELL [ "/bin/bash", "-c" ]


ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /home/docker
ENV HOME /home/docker
ENV CONDA_DIR $HOME/miniconda3
ENV PATH=$CONDA_DIR/bin:$PATH

#install miniconda via downloaded file (tested)
#COPY ./Miniconda3-latest-Linux-x86_64.sh $HOME/
#RUN chmod 777 Miniconda3-latest-Linux-x86_64.sh
#RUN bash $HOME/Miniconda3-latest-Linux-x86_64.sh -b -p ${CONDA_DIR}  #silent mode

#OR install via wget (tested) #Miniconda3-py37_4.9.2-Linux-x86_64.sh
RUN apt-get update && apt-get -y install wget
ENV MINICONDA_VERSION 4.9.2
ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh 
RUN chmod +x ~/miniconda.sh \ 
&& ~/miniconda.sh -b -p $CONDA_DIR \ 
&& rm ~/miniconda.sh

RUN apt-get install -y python3-pip libcgal-dev python-qt4 libgl1-mesa-glx locate  #To avoid "ModuleNotFoundError: No module named 'vtkOpenGLKitPython' " when importing vtk

# Create user "docker" with sudo powers
RUN useradd -m docker && \
    usermod -aG sudo docker && \
    echo '%sudo ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers && \
    cp /root/.bashrc /home/docker/ && \
    chown -R --from=root docker /home/docker

#python dependencies and path
ENV PYTHONPATH $HOME/SCRIPT/:$PYTHONPATH
ENV USER docker

COPY requirements_conda.txt . 
COPY requirements_pip.txt . 
ENV LIST_REQ="requirements_conda.txt"
#RUN for req in $(<${LIST_REQ}); do conda install -y ${req} ;done  #to test install conda one by one if install fails...
RUN conda install -y  --file requirements_conda.txt
RUN pip install -r requirements_pip.txt


# fixes
COPY SCRIPT_FIX/dataset_adapter.py $HOME/SCRIPT_FIX/
#fix (VTK selection does not work, ending up with a corrupt VTP-file)
#RUN cp -f SCRIPT_FIX/dataset_adapter.py  /usr/local/lib/python3.7/dist-packages/vtk/numpy_interface/dataset_adapter.py


USER docker
ENV PATH /home/docker/.local/bin:$PATH
# Avoid first use of sudo warning. c.f. https://askubuntu.com/a/22614/781671
RUN touch $HOME/.sudo_as_admin_successful

CMD [ "/bin/bash" ]
