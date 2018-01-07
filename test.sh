#!/usr/bin/env bash
for ((i=1;i<=10;i++)); do
    mpirun -np 27 ./multmat 3 $i 3 > results.dat
	octave-cli check_multmat.m
done