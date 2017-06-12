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
cooler = cooler.o
dvode = dvode.o

objects = $(BIN)$(dvode) $(BIN)$(cooler) $(BIN)$(constants) $(BIN)$(calculation) $(BIN)$(spectre) $(BIN)$(reader) $(BIN)$(main)

all : tevo

tevo : $(dvode) $(cooler) $(constants) $(calculation) $(spectre) $(reader) $(main)
		$(CC) $(FFLAGS) -o $(BIN)tevo $(objects)

main.o: spectre.o reader.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)main.f08 -o $(BIN)main.o

reader.o: spectre.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)reader.f08 -o $(BIN)reader.o

spectre.o:
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)spectre.f08 -o $(BIN)spectre.o

calculation.o: constants.o spectre.o reader.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)calculation.f08 -o $(BIN)calculation.o

constants.o:
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)constants.f08 -o $(BIN)constants.o

cooler.o: calculation.o
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)cooler.f08 -o $(BIN)cooler.o

dvode.o:
		$(CC) $(FFLAGS) -J $(OBJ) -c $(SRC)dvode_f90_m.f90 -o $(BIN)dvode.o

clean:
		rm $(BIN)*.o $(BIN)*.mod *.1 $(BIN)tevo
