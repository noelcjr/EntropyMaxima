#!/bin/bash

docker build -t open-insulin .
docker run -it -v "$(pwd)":/code open-insulin /bin/bash

