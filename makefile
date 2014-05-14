all: ftcs.c
	gcc ftcs.c -o ftcs

clean:
	rm ftcs ftcs*.txt
