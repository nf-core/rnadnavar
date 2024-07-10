FROM r-base

# install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    cmake

# install R packages
RUN R -e "install.packages(c('BiocManager','data.table', 'plyr', 'ggpubr', 'ggrepel'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('rtracklayer', 'ComplexHeatmap'))"
