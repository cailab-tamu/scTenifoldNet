
FROM continuumio/anaconda3:2018.12

LABEL Maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root

RUN apt-get update

RUN conda install -y -c r r-base r-glmnet r-ksamples r-ppcor r-pracma r-r.matlab r-cvTools

RUN conda install -y -c conda-forge unzip

WORKDIR /

RUN ls

COPY SINCERITIES.zip /

RUN unzip SINCERITIES.zip -d SINCERITIES

WORKDIR SINCERITIES/

RUN mkdir data/

RUN apt-get update && apt-get install time
