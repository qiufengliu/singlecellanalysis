FROM continuumio/miniconda3:4.9.2
ADD cellranger.tar.gz /opt/
RUN mv "$(find /opt -maxdepth 1 | grep cellranger)" "/opt/cellranger"

FROM continuumio/miniconda3:latest
COPY --from=0 /opt/cellranger /opt/cellranger
COPY src /opt/10x
ENV PATH="/opt/cellranger:${PATH}"
RUN conda install \
    --freeze-installed \
    --channel conda-forge \
    --channel bioconda \
    --file /opt/10x/environment.yaml \
    && conda clean --all \
    && rm /opt/conda/pkgs/ -rf
