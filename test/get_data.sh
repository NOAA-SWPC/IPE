#!/bin/bash

bucket="https://noaa-nws-wam-ipe-pds.s3.amazonaws.com"

if [ ! -d data ] ; then
  wget $bucket/NON-OPERATIONAL_RESEARCH/SWORD_DATA/IPE_STANDALONE_FIX.tar.gz
  tar xvf IPE_STANDALONE_FIX.tar.gz
  rm -rf  IPE_STANDALONE_FIX.tar.gz
else
  echo "data directory exists -- skipping"
fi

if [ ! -d regression ] ; then
  wget $bucket/NON-OPERATIONAL_RESEARCH/SWORD_DATA/IPE_STANDALONE_REGRESSION.tar.gz
  tar xvf IPE_STANDALONE_REGRESSION.tar.gz
  rm -rf  IPE_STANDALONE_REGRESSION.tar.gz
else
  echo "regression directory exists -- skipping"
fi
