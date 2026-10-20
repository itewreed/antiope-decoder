CC ?= cc
CFLAGS ?= -O3 -march=native

libviterbi.so: viterbi.c
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $<

clean:
	rm -f libviterbi.so
