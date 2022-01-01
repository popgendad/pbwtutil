FROM libpbwt:latest
RUN mkdir /usr/local/src/pbwtutil
COPY ./ /usr/local/src/pbwtutil/
WORKDIR /usr/local/src/pbwtutil
RUN make && make install
