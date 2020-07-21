FROM r-base:3.5.3

LABEL maintainer = "Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/remotes_2.1.1.tar.gz', type = 'source')"

RUN R -e "remotes::install_cran(pkgs = 'scTenifoldNet', quiet = TRUE)"

RUN R -e "remotes::install_cran(pkgs = 'reshape2', quiet = TRUE)"

RUN R -e "library(scTenifoldNet)"

RUN R -e "library(reshape2)"

COPY runSCTENIFOLDNET.R /

RUN mkdir data/

RUN apt-get update && apt-get install time
