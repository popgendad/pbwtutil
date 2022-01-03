CC      := gcc
VERSION := $(shell cat VERSION)
CFLAGS  := -Wall -O2 -D VERSION=$(VERSION)
LIBS    := -lz -lhts -lpbwt -lplink_lite
SRCS    := $(wildcard src/*.c)
OBJS    := $(SRCS:src/%.c=src/%.o)

all: pbwtutil

pbwtutil: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

src/%.o:src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

install:
	cp pbwtutil /usr/bin/pbwtutil

clean:
	rm src/*.o pbwtutil
