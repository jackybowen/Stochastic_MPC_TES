FROM esetdeveloper/julia:cbcdev

ADD . /home/eset_engine
WORKDIR /home/eset_engine

RUN julia -e 'import Pkg; Pkg.add(["Ipopt","ScikitLearn"]); Pkg.instantiate(); Pkg.precompile();'

ENV JULIA_CBC_LIBRARY_PATH "/CbcInstall/dist/lib"

RUN pip3 install --user flask_cors requests
