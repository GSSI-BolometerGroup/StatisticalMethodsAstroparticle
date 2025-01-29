TARGET      = ProduceData.x
OBJECTS     = ProduceData.cpp
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`
CC          = g++
CFLAGS	    = -O3 -Wall


INCLUDES = 
LIBD = 

all: ${TARGET}

ProduceData.x: ${OBJECTS}
	${CC} $(INCLUDES) $(LIBD) ${CFLAGS} -o $@ $+ ${ROOTCFLAGS} ${ROOTLIBS}

install: ${TARGET}
	@install ${TARGET}
	@echo "OK"

clean:
	rm -f ${TARGET}
