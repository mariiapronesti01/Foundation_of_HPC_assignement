CC = mpicc

SRCS = main.c read_write.c static_evolution.c ordered_evolution.c
OBJS = $(SRCS:.c=.o)

TARGET = my_program

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -fopenmp   $(OBJS) -o $(TARGET)

%.o: %.c
	$(CC) -fopenmp  -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

