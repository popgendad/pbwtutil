CC      := gcc
VERSION := $(shell cat VERSION)
CFLAGS  := -Wall -O2 -D VERSION=$(VERSION)
LIBS    := -lz -lhts -lpbwt -lplink_lite
SRCS    := $(wildcard *.c)
OBJS    := $(patsubst %.c, %.o, $(SRCS))

all: pbwtmaster

pbwtmaster: $(OBJS)
	$(CC) -O2 -o $@ $^ $(LIBS)

%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o pbwtmaster
