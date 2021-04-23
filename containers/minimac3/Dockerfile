FROM nfcore/base:1.9
LABEL authors=" " \
      description="Minimac3 container"

#Copy other binary dependencies
RUN apt-get clean && apt-get update && apt-get install -y libgomp1
COPY bin/minimac3 /usr/bin/
