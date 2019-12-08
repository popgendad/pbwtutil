FROM 552622921787.dkr.ecr.us-east-1.amazonaws.com/dnascience/amazonlinux2:latest

ARG NEXT_VERSION

RUN yum install -y htslib libpbwt gcc make libplink-lite

ADD . /pbwtmaster

RUN cd /pbwtmaster && \
    make && \
    mkdir -p /tmp/pbwtmaster/usr/bin/ && \
    echo ${NEXT_VERSION} > VERSION && \
    cp pbwtmaster /tmp/pbwtmaster/usr/bin/ && \
    fpm -s dir -t rpm -n pbwtmaster -v ${NEXT_VERSION} -C /tmp/pbwtmaster -p pbwtmaster_VERSION_ARCH.rpm -d libpbwt -d libplink-lite -d htslib . && \
    mv pbwtmaster*.rpm /rpms/

