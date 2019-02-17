CC = gcc
LINK = gcc
LIB = ar
RANLIB = ranlib
CFLAGS = -Wall -funsigned-char -I /opt/local/include -o $@
#CFLAGS += -g
CFLAGS += -O2
LFLAGS = $(CFLAGS) -L /opt/local/lib -lgmp

%.o:    %.c
	$(CC) -c $(CFLAGS) $<

all:	beal

beal:	beal.o
	$(LINK) $(LFLAGS) -o $@ $<
