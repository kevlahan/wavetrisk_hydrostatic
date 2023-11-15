#!/bin/bash

# This script helps working in Docker images while accessing the $HOME directory
# It works by adding /etc/group and /etc/passwd to a Docker image
# before starting a container from this personalized image

cat <<EOF
Usage :
$0 start [intel|gnu||nvidia] [--cpus NCPU]  # Start Docker container (Intel/GNU/Nvidia compilers)
$0 enter                                    # Enter the already running container
$0 admin                                    # Enter the running container with root privileges
$0 help                                     # This help
EOF

function cmd_help()
{
    echo
}

function cmd_install()
{
    cat <<EOF
Installing docker on your Ubuntu distribution using 'sudo apt install'.
You may be asked for your user password.
This password is the one you use to login. No root password required.
EOF
    sudo apt install docker.io
    apt show docker.io
}

function cmd_build()
{
    ARCH=$1 ; shift
    cat <<EOF
Building image with compiler $ARCH.
Docker options : $*
EOF
    TMPDIR=$(mktemp -d)
    cd $TMPDIR
    mkdir -p etc
    cp /etc/passwd /etc/group etc/
    cat <<'EOF' > Dockerfile
ARG ARCH=gnu
FROM gitlab-registry.in2p3.fr/ipsl/lmd/intro/awaca/${ARCH}/xios:2109
COPY etc/ /etc/
EOF
    docker image rm simple_physics
    DOCKER_BUILDKIT=1 docker build -t simple_physics --build-arg ARCH=$ARCH $* . && cat <<EOF
Docker image built successfully.
EOF
    cd -
    rm -rf $TMPDIR
}

function cmd_start()
{
    docker -v || cmd_install
    cmd_build $1 ;
    cat <<EOF
Starting container for architecture $1.
Use $0 enter to work from additional terminals in that same container.
EOF
    shift
    docker run --rm -it --name simple_physics --volume $HOME:$HOME --workdir=$PWD --user $(whoami) simple_physics $* bash
}

function cmd_enter()
{
    cat <<EOF
Entering running container.
Works only if you already started container using $0 start.
EOF
    docker exec -it --user $(whoami) $* simple_physics bash
}

function cmd_admin()
{
    cat <<EOF
Entering running container as root.
Works only if you already started container using $0 start.
EOF
    docker exec -it --user root simple_physics bash
}

CMD=$1 ; shift
cmd_$CMD $*
