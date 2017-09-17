#!/usr/bin/env bash

function c_echo(){
    local color=$1;
    shift
    local exp=$@;
    if ! [[ $color =~ '^[0-9]$' ]] ; then
       case $(echo $color | tr '[:upper:]' '[:lower:]') in
        black) color=0 ;;
        red) color=1 ;;
        green) color=2 ;;
        yellow) color=3 ;;
        blue) color=4 ;;
        magenta) color=5 ;;
        cyan) color=6 ;;
        white|*) color=7 ;; # white or invalid color
       esac
    fi
    tput setaf $color;
    echo $exp;
    tput sgr0;
}

function script_cmd() {
    c_echo YELLOW "\$ $@" >&2;
    "$@";
}

export -f script_cmd
export -f c_echo
