# Reproducible environment: Ubuntu 24.04 (gfortran 13, gawk 5.2, Python 3.12 with numpy/matplotlib).
#   docker build -t ampt-charm .
#   docker run --rm ampt-charm                 # runs make test (results, exact reproduction, quick production)
#   docker run --rm -it ampt-charm bash        # interactive
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran g++ make gawk curl ca-certificates unzip patch \
      python3 python3-numpy python3-matplotlib && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/ampt-charm
COPY . .
RUN make deps build
CMD ["make", "test"]
