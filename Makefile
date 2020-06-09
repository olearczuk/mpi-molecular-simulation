CC          := mpic++
CFLAGS      := -O3 -Wall -c
LFLAGS      := -O3 -Wall

all : body3

body3 : main.o utils.o
	$(CC) $(LFLAGS) -o $@ $^

main.o : main.cpp utils.h
	$(CC) $(CFLAGS) $<

#utils.o : utils.cpp utils.h
	#$(CC) $(CFLAGS) $<

clean :
	rm -f *.o body3
