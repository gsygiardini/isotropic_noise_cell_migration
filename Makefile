CC = gcc
CFLAGS = -g -Wall -Wextra -std=c11
DFLAGS = -DM_PI=3.14159265358979323846
EFLAGS = -DM_E=2.718281828459045
MFLAGS = -lm
GSLFLAGS = -lgsl -lgslcblas
OFLAGS = -O3
OMPFLAGS = -fopenmp
RENAME = -o
name ?= dsc

all: dsc

debug: dsc.c
	$(CC) $(DFLAGS) $(CFLAGS) $(IFLAGS) $(OFLAGS) $(EFLAGS) -DDEBUG dsc.c $(RENAME) debug_$(name) $(MFLAGS) $(GSLFLAGS)

dsc: dsc.c
	$(CC) $(DFLAGS) $(CFLAGS) $(IFLAGS) $(OFLAGS) $(EFLAGS) dsc.c $(RENAME) $(name) $(MFLAGS) $(GSLFLAGS) $(OMPFLAGS)



