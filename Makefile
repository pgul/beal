CC = gcc
LINK = gcc
LIB = ar
RANLIB = ranlib
CFLAGS = -Wall -funsigned-char -o $@
#CFLAGS += -g
CFLAGS += -O2
LFLAGS = $(CFLAGS)

%.o:    %.c
	$(CC) -c $(CFLAGS) $<

all:	beal

beal:	beal.o
	$(LINK) $(LFLAGS) -o $@ $<
