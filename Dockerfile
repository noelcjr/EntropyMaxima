FROM noelcjr/ccl_lectures:1.0

WORKDIR code
RUN pip install nose
RUN apt-get install wget

ENTRYPOINT ./bin/reinstall.sh && /bin/bash