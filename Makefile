CC=gcc
CFLAGS=-c

all: libplink.a

libplink.a: plink.o
	ar rcs $@ $<
	ranlib $@

plink.o: plink.c plink.h khash.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm plink.o libplink.a