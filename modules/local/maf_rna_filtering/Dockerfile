FROM python:3.9

# add ps (required by nextflow)
RUN apt-get --allow-releaseinfo-change update  && \
    apt-get install -y procps gfortran build-essential libblas-dev liblapack-dev libatlas-base-dev && \
        apt-get clean && \
            rm -rf /var/lib/apt/lists/*

# clone CApy repo
RUN pip install numpy && \
    git clone https://github.com/getzlab/CApy.git && \
    cd CApy && \
    python setup.py install && \
    cd .. && \
    mv CApy /tools/

ENV PATH="/tools/CApy/bin:${PATH}"

CMD ["python"]
