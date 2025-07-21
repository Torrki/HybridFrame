# Compiler and flags
CXX := g++
CXXFLAGS := -fPIC -Wall -O2 -std=c++20 -I../ODE/include
LDFLAGS := -shared -L../ODE/bin/c++

# Project structure
SRCDIR := src
BUILDDIR := bin
TARGET := $(BUILDDIR)/libautoma.so

# Source and object files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))

# Default rule
all: $(TARGET)

# Linking dynamic library
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILDDIR)
	$(CXX) $(LDFLAGS) -o $@ $^ -lode

# Compiling source files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean

