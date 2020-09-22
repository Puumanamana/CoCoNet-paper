FROM continuumio/miniconda3:latest
LABEL author="carisdak@hawaii.edu"

RUN apt-get update && apt-get install -y procps git libncursesw5 && apt-get clean -y

COPY conda-camisim.yaml /
RUN conda env create -f /conda-camisim.yaml && conda clean -a

#=================================================#
#============== CAMISIM from GitHub ==============#
#=================================================#

RUN git clone https://github.com/CAMI-challenge/CAMISIM.git \
    && mv CAMISIM /opt

WORKDIR /workspace

ENV PATH /opt/conda/envs/camisim/bin:/opt/CAMISIM:$PATH
