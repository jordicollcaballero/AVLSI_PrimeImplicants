CXX=g++
RM=rm -f
CPPFLAGS=-g -std=c++11
LDFLAGS=-g
LDLIBS=-lboost_program_options 

SRCS=main.cpp implicant.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: run

run: $(OBJS)
	$(CXX) $(LDFLAGS) -o run $(OBJS) $(LDLIBS) 

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) simulation
