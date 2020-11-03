FROM pytorch/pytorch:latest
LABEL author="carisdak@hawaii.edu"

RUN apt-get update && apt-get install -y gcc
RUN pip install pandas seaborn matplotlib biopython scikit-learn h5py numpy parameter-sherpa 

WORKDIR /workspace

