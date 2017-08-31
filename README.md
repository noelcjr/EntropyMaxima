# Entropy Maxima

"The Entropy Maximum is a corollary of the Energy minimum" David Chandler
"but there are many local maxima, just like there multiple local minima." Noel Carrascal

## How to run

The project has a docker file, and a bin folder with 2 useful scripts.
To build, run and connect to the docker container running the script run

```bash
./bin/run.sh
```

## Tests
There are unit tests which you should run when changing code to verify that you have not introduced
any problems.

To do so run this inside the __container__
```bash
./bin/run-unit.sh
```

## Installing

You can use 

```bash
./bin/reinstall.sh
```
to install the scripts in the bin so that you will be able to run from anywhere in the container.

## Structure

WIP