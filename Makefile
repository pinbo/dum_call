# Makefile for dum_call

CC     = gcc
# Default libraries to link if configure is not used
LIBS = ../libhts.a -lz -lm -lbz2 -llzma -lcurl -lcrypto
CFLAGS   = -g -Wall -O2

dum_call: dum_call.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
