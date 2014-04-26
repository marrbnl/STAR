#!/bin/bash

if [ -z $1 ]; then
    comment='local backup'
else
    comment=$1
fi

wdir=$PWD
dropbox=~/Dropbox/git-repository
cp ./*.C $dropbox/STAR/QA/.
cp ./*.sh $dropbox/STAR/QA/.
cd $dropbox
git add .
git status
git commit -m "${comment}"
cd $wdir

github=~/Work/git-repository/STAR
cp ./*.C $github/QA/.
cp ./*.sh $github/QA/.
cd $github
git add .
git status
git commit -m "${comment}"
cd $wdir