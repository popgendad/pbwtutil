CC=gcc
CFLAGS=-Wall
LIBS=-L. -lplink

all: libplink.a

libplink.a: plink.o
	ar rcs $@ $<
	ranlib $@

plink.o: plink.c plink.h khash.h
	$(CC) $(CFLAGS) -fPIC -c $<

test-api: test-api.c libplink.a
	$(CC) $(CFLAGS) -I. -o $@ $< $(LIBS)

clean:
	rm -f *.o libplink.a test-api