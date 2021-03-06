FROM continuumio/miniconda3:latest
LABEL author="carisdak@hawaii.edu"

RUN apt-get update && apt-get install -y procps git && apt-get clean -y

COPY conda-station-aloha.yaml /
RUN conda env create -f /conda-station-aloha.yaml && conda clean -a
COPY diamond /opt/conda/envs/nf-station-aloha/bin

#=================================================#
#============ stampede-clustergenomes ============#
#=================================================#

RUN git clone https://bitbucket.org/MAVERICLab/stampede-clustergenomes.git \
    && mv stampede-clustergenomes /opt

WORKDIR /workspace

ENV PATH /opt/conda/envs/nf-station-aloha/bin:/opt/stampede-clustergenomes/bin:$PATH

RUN vdb-config
