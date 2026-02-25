CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra

airgun: airgun.c tables.h
	$(CC) $(CFLAGS) -o $@ airgun.c -lm

tables.h: generate_tables.py
	python3 generate_tables.py

clean:
	rm -f airgun

.PHONY: clean
