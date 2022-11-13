#!/bin/bash

for file in ./dump/*.dat;
do
	./scripts/process-data.wls $file >> ./stats/stats.dat &
done
wait
