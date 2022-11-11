# Makefile for compiling e_m2dRNAs source code
CC=g++
LD=g++
RM=rm -f
CPPFLAGS=-w -ansi -pedantic -g -O3 -fopenmp -std=c++11
OBJS:=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
MAIN=e_m2dRNAs
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm -fopenmp  -I /usr/local/include/ViennaRNA/ -L /usr/local/share/ViennaRNA/ -lRNA -lm  -fno-lto -O3 -static -std=c++11
%.o: %.c global.h rand.h
	$(CC) $(CPPFLAGS) -c $<
clean:
	$(RM) $(OBJS)

