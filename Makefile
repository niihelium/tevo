CC=gfortran
CFLAGS=-I.
SRC = ./src/
OBJ = ./include
BIN = ./bin/
FFLAGS = -o2 -ggdb
spectre = spectre.o
reader = reader.o
main = main.o
calculation = calculation.o
constants = constants.o

objects =  $(BIN)$(constants) $(BIN)$(calculation) $(BIN)$(spectre) $(BIN)$(reader) $(BIN)$(main)

all : tevo

tevo : $(constants) $(calculation) $(spectre) $(reader) $(main)
		$(CC) $(FFLAGS) -o $(BIN)tevo $(objects)

main.o: spectre.o reader.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)main.f90 -o $(BIN)main.o

reader.o: spectre.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)reader.f90 -o $(BIN)reader.o

spectre.o:
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)spectre.f90 -o $(BIN)spectre.o

calculation.o: constants.o spectre.o reader.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)calculation.f90 -o $(BIN)calculation.o

constants.o:
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)constants.f90 -o $(BIN)constants.o

clean:
		rm $(BIN)*.o $(BIN)*.mod *.1 $(BIN)tevo
