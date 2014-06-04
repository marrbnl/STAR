#!/bin/bash

if [ -z $1 ]; then
    comment='local backup'
else
    comment=$1
fi

wdir=$PWD
basedir=${PWD##*/}


dropbox=~/Dropbox/git-repository
if [ ! -d $dropbox/STAR/${basedir} ]; then
    mkdir -pv $dropbox/STAR/${basedir}
fi

cp ./*.C $dropbox/STAR/${basedir}/.
cp ./*.sh $dropbox/STAR/${basedir}/.
cd $dropbox
git add .
git status
git commit -m "${comment}"
cd $wdir

github=~/Work/git-repository/STAR
if [ ! -d $github/${basedir} ]; then
    mkdir -pv $github/${basedir}
fi
cp ./*.C $github/${basedir}/.
cp ./*.sh $github/${basedir}/.
cd $github
git add .
git status
git commit -m "${comment}"
git push origin master
cd $wdir