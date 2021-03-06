CXX=g++
RM=rm -f
CPPFLAGS=-O3 -std=c++11
LDFLAGS=-O3
LDLIBS=-lboost_program_options 

SRCS=main.cpp implicant.cpp matrix.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: run

run: $(OBJS)
	$(CXX) $(LDFLAGS) -o run $(OBJS) $(LDLIBS) 

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) simulation
