FROM noelcjr/ccl_lectures:1.0

WORKDIR code
RUN pip install nose
RUN pip install jinja2
RUN apt-get install wget
RUN echo "source ./bin/functions.sh" >> ~/.bashrc

ENTRYPOINT ./bin/reinstall.sh && /bin/bash
