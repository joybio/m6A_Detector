#!/bin/bash
sed -i 's/[ ][ ]*/\t/g' uniq.stast 
sed -i 's/\s\+//g' uniq.stast 
