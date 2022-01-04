FROM popgendad/libpbwt:amd64
RUN mkdir -p /usr/local/src/pbwtutil
COPY ./ /usr/local/src/pbwtutil/
WORKDIR /usr/local/src/pbwtutil
RUN ldconfig
RUN make && make install
