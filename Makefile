CC=gcc
CFLAGS=-Wall -O2 -fpic
PREFIX_DIR=/usr/local
H_DIR=$(PREFIX_DIR)/include/
L_DIR=$(PREFIX_DIR)/lib/

all: libplink.so

libplink.so: plink.o
	gcc -shared -o $@ $<

plink.o: plink.c
	$(CC) $(CFLAGS) -c $<

install:
	cp libplink.so ${L_DIR}
	cp plink.h ${H_DIR}

clean:
	rm -f plink.o libplink.so
