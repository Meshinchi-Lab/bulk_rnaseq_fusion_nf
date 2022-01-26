FROM ghcr.io/stjude/cicero:v1.8.1

#Updating Ubuntu packages
RUN apt-get update && yes|apt-get upgrade

#Overwrite the ENTRYPOINT from CICERO image
#Causes issues with AWS batch 
ENTRYPOINT []