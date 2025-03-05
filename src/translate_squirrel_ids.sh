#!/bin/bash

sed -E 's/^(>[^.]+)\.[0-9]+/\1/' rawdata/HLictTri2/genome/HLictTri2.fa > rawdata/HLictTri2/genome/HLictTri2.translated.fa
