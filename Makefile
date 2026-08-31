CC ?= clang
CFLAGS ?= -std=c11 -Wall -Wextra -O3 -Iinclude
LDFLAGS ?= -lm

SRC = src/main.c \
      src/rqm_math.c \
      src/rqm_classical.c \
      src/rqm_hilbert_pathology.c \
      src/rqm_ensembles.c \
      src/rqm_born_entropy.c \
      src/rqm_unitary.c \
      src/rqm_equilibration.c \
      src/rqm_classical_limit.c \
      src/rqm_triple_equilibria.c \
      src/rqm_ui.c \
      src/rqm_tests.c

OBJ = $(SRC:.c=.o)
TARGET = rqm

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TARGET)
	./$(TARGET) --test

clean:
	rm -f $(OBJ) $(TARGET) $(TARGET).exe

.PHONY: all test clean
