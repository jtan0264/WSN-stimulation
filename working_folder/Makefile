CC = mpicc #Compiler name
CFLAGS = -lpthread -lm #Compiler flags
NAME = main.out # Output file name

all: main.c structure.c calculation.c seafloor_nodes.c  base_balloon.c 
	$(CC) -o $(NAME) main.c structure.c calculation.c seafloor_nodes.c base_balloon.c $(CFLAGS)