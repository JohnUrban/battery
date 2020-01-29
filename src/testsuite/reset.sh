#!/bin/bash

while read delfile; do
  if [ -f $delfile ]; then rm $delfile; 
  elif [ -d $delfile ]; then rm -r $delfile;
  fi
done < util/reset-by-deleting-these-files.fofn

