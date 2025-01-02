CC      = cc
CXX     = c++
CFLAGS  = -g -O3 -Wall
LIBS    = -lpthread
LDFLAGS = -g

OBJECTS =	main.o moon.o planet.o sky0.o sky1.o skyfast.o skyio.o sky-site.o sky-time.o star.o sun.o vectors3d.o

all:		MoonList

MoonList:	$(OBJECTS)
		$(CXX) $(OBJECTS) $(CFLAGS) $(LIBS) -o MoonList

%.o: %.cpp
		$(CXX) $(CFLAGS) -c -o $@ $<

clean:
		$(RM) MoonList *.o *.d *.bak *~

install:
		install -m 755 MoonList /usr/local/bin/

