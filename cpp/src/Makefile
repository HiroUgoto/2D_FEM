TARGET_EXEC = FEM
CXX = g++
CPPFLAGS = -std=c++17 -Wall --pedantic-errors
LIBRARY = /usr/local/include/eigen3

SRCS := $(shell find . -name '*.cpp')
OBJS := $(SRCS:%=%.o)

$(TARGET_EXEC): $(OBJS) all.h all.h.gch
	$(CXX) $(CPPFLAGS) -include all.h $(OBJS) -o $@

all.h.gch: all.h
	$(CXX) $(CPPFLAGS) -x c++-header -o $@ $<

%.cpp.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@ -I$(LIBRARY)

clean:
	rm -f $(TARGET_EXEC)
	rm -f ./all.h.gch
	rm -f $(OBJS)

.PHONY: clean
