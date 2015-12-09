#!/bin/bash

for i in $(seq 4 200); do 
    echo "i:$i"
    time queens $i
done
