# escape=`
# The above line is a parser directive, not a comment.


# docker pull jupyter/base-notebook
# sha256:8c903974902b0e9d45d9823c2234411de0614c5c98c4bb782b3d4f55b3e435e6

FROM jupyter/base-notebook



USER root



RUN apt-get update -y
RUN apt-get update && apt-get upgrade -y
RUN apt-get install nano
RUN apt-get install git --fix-missing  -y

RUN apt-get install make --fix-missing  -y
RUN apt-get install gcc --fix-missing  -y


RUN conda install -c conda-forge pymol-open-source
#RUN apt-get install pymol #fails
RUN apt-get install pymol --fix-missing -y

# we probably need these, but not tested yet
# RUN conda install -c schrodinger pymol
# RUN conda install -c schrodinger pymol-psico
# RUN conda install -c speleo3 tmalign


RUN wget  -v http://fasta.bioch.virginia.edu/wrpearson/fasta/CURRENT/fasta36-linux64.tar.gz
RUN tar -xvzf fasta36-linux64.tar.gz
RUN rm -r fasta36-linux64.tar.gz
RUN cd fasta-36.3.8i/src; `
    make -f ../make/Makefile.linux_sse2 all
    
USER jovyan
RUN /opt/conda/bin/python3 -m pip install pot `
 biopython==1.81 `
 umap-learn==0.5.3 `
 multiprocess `
 deprecated `
  --upgrade setuptools


  
# RUN /opt/conda/bin/python3 -m pip install  git+https://github.com/CamaraLab/CAJAL.git




#to avoid the jupyter server from automatically launching at start up
# ENTRYPOINT ["sleep", "infinity"]

# CMD ["sh", "-c", "start-notebook.sh --NotebookApp.token=''"]

 # test
 # ../bin/fasta36 -q ../seq/mgstm1.aa ../seq/prot_test.lseg

 #  docker run -it -e  GRANT_SUDO=yes --user root jupyter/base-notebook /bin/bash






####### old one:

# # escape=`
# # The above line is a parser directive, not a comment.


# # docker pull jupyter/base-notebook
# # sha256:8c903974902b0e9d45d9823c2234411de0614c5c98c4bb782b3d4f55b3e435e6
# FROM jupyter/base-notebook



# USER root

# RUN wget  -v http://fasta.bioch.virginia.edu/wrpearson/fasta/CURRENT/fasta36-linux64.tar.gz



# RUN apt-get update -y
# RUN apt-get update && apt-get upgrade -y
# #RUN  apt-get install build-essential -y --fix-missing    #this gives error I havne't resolved
# RUN apt-get update -y && apt-get install git-all -y --fix-missing
# RUN apt-get install nano

# #trying this out
# RUN  apt-get update && apt-get install -y apt-transport-https
# RUN apt-get update && apt-get install pymol -y --fix-missing
# #takes 8 minutes..
# # seems to be the shell python .. ?

# RUN whereis cd
# RUN whereis make


# # commented out for now bc of mysterious bug
# # RUN tar -xvzf fasta36-linux64.tar.gz
# # RUN rm -r fasta36-linux64.tar.gz
# # RUN cd fasta-36.3.8i/src; `
# #     make -f ../make/Makefile.linux_sse2 all
# # USER jovyan
# # RUN /opt/conda/bin/python3 -m pip install pot `
# #  biopython==1.81 `
# #  umap-learn==0.5.3 `
# #  multiprocess `
# #  deprecated `
# #   --upgrade setuptools
# # RUN /opt/conda/bin/python3 -m pip install  git+https://github.com/CamaraLab/CAJAL.git




# #to avoid the jupyter server from automatically launching at start up
# # ENTRYPOINT ["sleep", "infinity"]

# # CMD ["sh", "-c", "start-notebook.sh --NotebookApp.token=''"]

#  # test
#  # ../bin/fasta36 -q ../seq/mgstm1.aa ../seq/prot_test.lseg

#  #  docker run -it -e  GRANT_SUDO=yes --user root jupyter/base-notebook /bin/bash
