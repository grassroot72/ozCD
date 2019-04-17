# ozCD makefile

CC = gcc
CFLAGS = -O2 -Wall -I /usr/include -D OZ_SINGLE -D BARY_GEPP
LDFLAGS = -L /usr/lib -lm
OBJS = shape.o \
       simplex.o \
       ccd.o \
       main.o
EXE = ozccd


all: $(EXE)

%.o: %.c
	$(CC) -c $< $(CFLAGS)

$(EXE): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f $(OBJS) $(EXE)
