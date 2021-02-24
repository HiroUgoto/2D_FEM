TARGET_EXEC = FEM
CXX = g++
CPPFLAGS = -std=c++17 -Wall --pedantic-errors -O2 -march=native
LIBRARY = /usr/local/include/eigen3
# LIBRARY = /usr/include/eigen3

BUILD_DIR = obj

SRCS := $(shell find . -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(SRCS:%=$(BUILD_DIR)/%.d)

$(TARGET_EXEC): $(OBJS) all.h all.h.gch
	$(CXX) $(CPPFLAGS) -include all.h $(OBJS) -o $@

-include $(DEPS)

all.h.gch: all.h
	$(CXX) $(CPPFLAGS) -x c++-header -o $@ $<

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) -c -MMD -MP $< -o $@ -I$(LIBRARY)

clean:
	rm -f $(TARGET_EXEC)
	rm -f ./all.h.gch
	rm -rf $(BUILD_DIR)

.PHONY: clean

MKDIR_P = mkdir -pv
