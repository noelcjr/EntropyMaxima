FROM noelcjr/ccl_lectures:1.0

WORKDIR code
RUN pip install nose

CMD ["bash"]