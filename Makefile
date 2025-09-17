# Makefile for two_img_lut

VERSION = 1.0
BINDIR = $(HOME)/.local/bin
PREFIX = /usr/local

# Default LUT size if not provided
DEFAULT_LUT_SIZE = 33

CC = gcc
CFLAGS = -g -std=c11 -Wall -Wextra -pedantic -D_POSIX_C_SOURCE=199309L
LIBS = -lm

SRCDIR = src
EXTERNAL = external/stb
SRC = $(SRCDIR)/two_img_lut.c
# You might also have other .c in src

INCLUDES = -I$(EXTERNAL)

BIN = two_img_lut

.PHONY: all clean install

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN) $(SRC) $(LIBS)

install: $(BIN)
	mkdir -p $(BINDIR)
	cp $(BIN) $(BINDIR)
	@echo "Installed $(BIN) to $(BINDIR)"

clean:
	rm -f $(BIN)
