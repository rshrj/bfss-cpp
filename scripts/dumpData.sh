#!/bin/bash

for n in {1..96};
do
	./bin/matrix | tee ./dump/$(shuf -i 1000000-9999999 -n 1).dat &
done
wait
