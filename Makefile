CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -Werror -std=c99 -ggdb3 -Ofast 
LDFLAGS = -lm
SRC = main.c functions.c
OBJ = $(SRC:.c=.o)
EXEC = MatrixSolver

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< 

clean:
	rm -f $(OBJ) $(EXEC)