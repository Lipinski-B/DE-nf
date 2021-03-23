################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.9.2

################## METADATA #######################
LABEL base_image="continuumio/miniconda3"
LABEL version="4.9.2"
LABEL software="DE-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for DE-nf"
LABEL about.home="http://github.com/Lipinski-B/DE-nf"
LABEL about.documentation="http://github.com/Lipinski-B/DE-nf/README.md"
LABEL about.license="GNU-3.0"

################## INSTALLATION ######################
COPY environnement.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n DE-nf -f /environnement.yml && conda clean -a
RUN apt-get install -y curl
RUN apt-get install -y cmake python-pip python-dev
RUN pip install cget 
ENV PATH /opt/conda/envs/DE-nf/bin:$PATH
RUN Rscript -e "install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')"
