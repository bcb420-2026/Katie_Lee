FROM risserlin/bcb420-base-image:winter2026-arm64

RUN R -e 'BiocManager::install(c("DESeq2", "enrichplot"))'
RUN R -e 'install.packages("pheatmap")'