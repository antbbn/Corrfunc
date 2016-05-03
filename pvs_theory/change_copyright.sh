#!/bin/bash

# run with
# find . -name "*.[ch]" -exec ./change_copyright.sh '{}' ';'

# OLDHEADER is the header that is currently in the file, 
# will be substituted with whatever is into header.txt
#OLDHEADER=../header.txt
OLDHEADER=header.txt
# +2 because lines start from 1 in tail and the first one is the filename
OLDHEADERLEN=`wc -l $OLDHEADER | awk '{print $0+2}'`

echo "/* File: $(basename $1) */" > new_file
tail -n +$OLDHEADERLEN $1 | cat header.txt - >> new_file && mv new_file $1
