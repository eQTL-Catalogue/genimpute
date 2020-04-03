FROM nfcore/base
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for the genotype imputation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/imputation-1.0dev/bin:$PATH
RUN apt-get update && apt-get install -y libgomp1
COPY bin/eagle /usr/bin/
COPY bin/minimac4 /usr/bin/