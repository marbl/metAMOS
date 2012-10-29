#!/bin/bash
HOST='ftp.cbcb.umd.edu'
USER='anonymous'
PASSWD='NA'
LOCALFILE=$1
REMOTEFILE=$2
#
ftp -n -v $HOST <<EOF
user ${USER} ${PASSWD}
binary
cd /pub/data/treangen
put ${LOCALFILE} ${REMOTEFILE}
quit
EOF