# escape=`
# The above line is a parser directive, not a comment.


# docker pull jupyter/base-notebook
# sha256:8c903974902b0e9d45d9823c2234411de0614c5c98c4bb782b3d4f55b3e435e6
FROM jupyter/base-notebook



# USER root


# RUN apt-get update -y
# RUN apt-get upgrade -y


# RUN  apt-get install build-essential -y

# RUN apt-get install git-all -y

ENTRYPOINT ["sleep", "infinity"]

USER jovyan
RUN /opt/conda/bin/python3 -m pip install pot --upgrade setuptools





 # test
 # ../bin/fasta36 -q ../seq/mgstm1.aa ../seq/prot_test.lseg

 #  docker run -it -e  GRANT_SUDO=yes --user root jupyter/base-notebook /bin/bash
