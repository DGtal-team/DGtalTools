# DGtalTools based on the prebuilt DGtal image
FROM docker.io/ker69/dgtal:latest

# Build DGtalTools
RUN git clone --depth 1 https://github.com/DGtal-team/DGtalTools.git && \
    cd DGtalTools && \
    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build  && \
    cmake --install build


