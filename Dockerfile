FROM nfcore/base:1.9
LABEL authors=" " \
      description="Docker image containing all software requirements for the nf-core/imputation pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-imputation-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-imputation-1.0dev > nf-core-imputation-1.0dev.yml

#Copy other binary dependencies
RUN apt-get clean && apt-get update && apt-get install -y libgomp1
COPY bin/eagle /usr/bin/
COPY bin/minimac4 /usr/bin/
COPY bin/GenotypeHarmonizer-1.4.20/ /usr/bin/
