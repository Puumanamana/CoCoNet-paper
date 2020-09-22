FROM continuumio/miniconda3:latest
LABEL author="carisdak@hawaii.edu"

RUN apt-get update && apt-get install -y procps git && apt-get clean -y

COPY conda-python.yaml /
RUN conda env create -f /conda-python.yaml && conda clean -a

WORKDIR /workspace

ENV PATH /opt/conda/envs/python/bin:$PATH
