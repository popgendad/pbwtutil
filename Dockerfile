FROM 552622921787.dkr.ecr.us-east-1.amazonaws.com/dnascience/amazonlinux2:latest

ARG NEXT_VERSION

RUN yum install -y htslib libpbwt gcc make libplink-lite

ADD . /ancmatch

RUN cd /ancmatch && \
    make && \
    mkdir -p /tmp/ancmatch/usr/bin/ && \
    echo ${NEXT_VERSION} > VERSION && \
    cp ancmatch /tmp/ancmatch/usr/bin/ && \
    fpm -s dir -t rpm -n ancmatch -v ${NEXT_VERSION} -C /tmp/ancmatch -p ancmatch_VERSION_ARCH.rpm -d libpbwt -d libplink-lite -d htslib . && \
    mv ancmatch*.rpm /rpms/

