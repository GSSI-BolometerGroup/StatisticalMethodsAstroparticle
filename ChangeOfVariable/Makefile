TARGET      = Quadratic.x
OBJECTS     = Quadratic.cpp
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`
CC          = g++
CFLAGS	    = -O3 -Wall


INCLUDES = 
LIBD = 

all: ${TARGET}

Quadratic.x: ${OBJECTS}
	${CC} $(INCLUDES) $(LIBD) ${CFLAGS} -o $@ $+ ${ROOTCFLAGS} ${ROOTLIBS}

install: ${TARGET}
	@install ${TARGET}
	@echo "OK"

clean:
	rm -f ${TARGET}
