CC=gcc
CFLAGS=-Wall -O2 -fpic
PREFIX_DIR=/usr
H_DIR=$(PREFIX_DIR)/include/
L_DIR=$(PREFIX_DIR)/lib/
SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c, %.o, $(SRCS))

all: libplink.so

libplink.so: $(OBJS)
	gcc -shared -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $<

install:
	cp libplink.so ${L_DIR}
	cp plink.h ${H_DIR}

clean:
	rm -f $(OBJS) libplink.so
