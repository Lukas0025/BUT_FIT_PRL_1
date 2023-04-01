# Implementation of paraller splitting in MPI
### PRL Project 1
### Author: Lukas Plevac <xpleva07@vutbr.cz>

This repository implementing paraller splitting algoritm in cpp and MPI. Parralel splitting algoritm select some M
and then split array S to three sub arrays L, E and G where L is *{X | X in S and X < M}*, E is *{X | X in S and X = M}* and G is *{X | X in S and X > M}*.

Shortly L is array with number smaller that M, E is array of same numbers as M and G is array of numbers bigger that M.

This Repository is exmaple how to use MPI, but is not emxample what do with MPI, because this algoritm is in real slower that seqvetion version, because comanication is match expensive. To make this program eficient numbers processing on ranks must be mutch harder, sort numbers to arrays is too match simple.

## Building program

```sh
mpic++ --prefix /usr/local/share/OpenMPI -o parsplit parsplit.cpp
```

## Running program

Simpli run program with mpi run. Program need imput file with name **numbers** in same dir. Number of proccess must by devider of number of numbers (eg. for 10 numbers can be 1,2,5,10 proccesses).

```sh
dd if=/dev/random bs=1 count=10 of=numbers # create imput file with 10 numbers
mpirun --prefix /usr/local/share/OpenMPI --oversubscribe -np 5 parsplit # run program on 5 processes
```

## Run using run.sh

run.sh contains basic command to generate file build program and run program with file.

```sh
bash run.sh #with 10 numbers
bash run.sh 64 #with 64 numbers
```
