FROM nfcore/base:1.9
LABEL authors="Kaur ALasoo" \
      description="Beagle 5.2 container"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-imputation-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-imputation-1.0dev > nf-core-imputation-1.0dev.yml

#Copy other binary dependencies
COPY bin/beagle.21Apr21.304.jar /usr/bin/
COPY bin/bref3.21Apr21.304.jar /usr/bin/
COPY bin/conform-gt.24May16.cee.jar /usr/bin/
