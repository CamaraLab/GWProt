# escape=`
# The above line is a parser directive, not a comment.


# docker pull jupyter/base-notebook
# sha256:8c903974902b0e9d45d9823c2234411de0614c5c98c4bb782b3d4f55b3e435e6


FROM jupyter/base-notebook@sha256:66830d423f076d186560fb52fe32e6af637888f85b6c9b942fb0b0c36e869b7b
# Jupyter releases two versions of each Docker image, the one above targets the amd64 architecture, the one below
# targets the arm64 architecture.
# FROM jupyter/base-notebook@sha256:e471b4bf9680c3fb24c799a23fb7240def04b51e88913f877b9a6b411eaa8be2


USER root
RUN apt-get update -y
RUN apt-get update && apt-get upgrade -y
RUN apt-get install nano
RUN apt-get install git --fix-missing  -y
RUN apt-get install make --fix-missing  -y
RUN apt-get install cmake --fix-missing  -y
RUN apt-get install gcc --fix-missing  -y
RUN apt-get install ffmpeg libsm6 libxext6  -y
RUN  conda update -n base -c conda-forge conda -y
RUN conda install numpy

RUN apt-get update && apt-get install -y `
    bash `
    coreutils `
    --no-install-recommends && `
    apt-get clean && `
    rm -rf /var/lib/apt/lists/*

RUN apt-get install git build-essential python3-dev libglew-dev `
  libpng-dev libfreetype6-dev libxml2-dev `
  libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev --fix-missing -y


RUN git clone https://github.com/schrodinger/pymol-open-source.git
RUN git clone https://github.com/rcsb/mmtf-cpp.git
RUN mv mmtf-cpp/include/mmtf* pymol-open-source/include/
WORKDIR pymol-open-source

RUN python3 setup.py build install

RUN apt-get install python3-pip --fix-missing  -y



WORKDIR /home/jovyan

RUN wget  -v http://fasta.bioch.virginia.edu/wrpearson/fasta/CURRENT/fasta36-linux64.tar.gz
RUN tar -xvzf fasta36-linux64.tar.gz
RUN rm -r fasta36-linux64.tar.gz
RUN cd fasta-36.3.8i/src; `
    make -f ../make/Makefile.linux_sse2 all
   
RUN echo "export PATH=\"/home/jovyan/fasta-36.3.8i/bin:\$PATH\"" >> .bashrc
USER jovyan

export PATH="/home/jovyan/.local/bin:$PATH"

RUN  pip install pot `
 cajal `
 biopython==1.81 `
 umap-learn==0.5.3 `
 multiprocess `
cython `
blosum `
sparse `
Bio `
statistics `
cajal `
scikit-learn `
multiprocess `
threadpoolctl `
  --upgrade setuptools



#still needs to import GWProt