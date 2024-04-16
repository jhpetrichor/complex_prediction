CC = g++
CFLAGS = -O3 -Wall
OBJECTS = main.o  gen_cmplex.o gen_qclq.o
TARGET = qclq_cmplx

$(TARGET) : $(OBJECTS) Makefile
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET)

main.o: main.cpp global.h
	$(CC) -c $(CFLAGS) main.cpp

gen_cmplex.o: gen_cmplex.cpp global.h
	$(CC) -c $(CFLAGS) gen_cmplex.cpp
	
gen_qclq.o: gen_qclq.cpp global.h
	$(CC) -c $(CFLAGS) gen_qclq.cpp
	
clean :
	rm *.o
	
