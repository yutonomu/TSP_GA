RM = rm
CC	= gcc
CFLAGS	= -O3 #-Wall -std=c++11 
LDFLAGS	= -lm 


####
PROG = tsp_nn
OBJS = functions.o
all: $(PROG)

tsp_nn:	tsp_nn.c $(OBJS)
	$(CC) $(CFLAGS) -o tsp_nn $(OBJS) $< $(LDFLAGS)

functions.o:	functions.c functions.h

clean:
	$(RM) -f *.o *~ $(PROG)
