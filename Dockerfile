FROM quay.io/broadinstitute/viral-baseimage:0.1.15

LABEL maintainer "ilya_shl@alum.mit.edu"

# to build:
#   docker build .
#
# to run:
#   docker run --rm <image_ID> "<command>.py subcommand"
#
# to run interactively:
#   docker run --rm -it <image_ID>
#
# to run with GATK and/or Novoalign:
#   Download licensed copies of GATK and Novoalign to the host machine (for Linux-64)
#   export GATK_PATH=/path/to/gatk/
#   export NOVOALIGN_PATH=/path/to/novoalign/
#   docker run --rm -v $GATK_PATH:/gatk -v $NOVOALIGN_PATH:/novoalign -v /path/to/dir/on/host:/user-data <image_ID> "<command>.py subcommand"

ENV \
	INSTALL_PATH="/opt/cms" \
	CMS_PATH="/opt/cms/source" \
	MINICONDA_PATH="/opt/miniconda" \
	CONDA_DEFAULT_ENV=cms-env
ENV \
	PATH="$CMS_PATH:$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV/bin:$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
	CONDA_PREFIX=$MINICONDA_PATH/envs/$CONDA_DEFAULT_ENV \
	JAVA_HOME=$MINICONDA_PATH

# COPY docker/install-*.sh /opt/docker/

# # System packages, Google Cloud SDK, and locale
# # ca-certificates and wget needed for gosu
# # bzip2, liblz4-toolk, and pigz are useful for packaging and archival
# # google-cloud-sdk needed when using this in GCE
# RUN /opt/docker/install-apt_packages.sh

# # Set default locale to en_US.UTF-8
# ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"

# # grab gosu for easy step-down from root
# RUN /opt/docker/install-gosu.sh

# # install miniconda3 with our default channels and no other packages
# ENV MINICONDA_PATH="/opt/miniconda"
# RUN /opt/docker/install-miniconda.sh

# Prepare cms user and installation directory
# Set it up so that this slow & heavy build layer is cached
# unless the requirements* files or the install scripts actually change
WORKDIR $INSTALL_PATH
RUN conda create -n $CONDA_DEFAULT_ENV python=3.6
RUN echo "source activate $CONDA_DEFAULT_ENV" > ~/.bashrc
RUN hash -r
COPY docker/install-cms.sh $CMS_PATH/docker/
COPY requirements-cms.txt $CMS_PATH/
RUN $CMS_PATH/docker/install-cms.sh

# Copy all of the source code into the repo
# (this probably changes all the time, so all downstream build
# layers will likely need to be rebuilt each time)
COPY . $CMS_PATH/

# Volume setup: make external tools and data available within the container
VOLUME ["/user-data"]
ENV \
	CMS_DOCKER_DATA_PATH="/user-data"

# This not only prints the current version string, but it
# also saves it to the VERSION file for later use and also
# verifies that conda-installed python libraries are working.
#RUN /bin/bash -c "set -e; echo -n 'cms version: '; scans.py --version"

## A wrapper script to load the cms environment, switch to
## a non-root user, and then run any commands desired
#ENTRYPOINT ["/opt/cms/source/docker/gosu-entrypoint.sh"]

CMD ["/bin/bash"]
