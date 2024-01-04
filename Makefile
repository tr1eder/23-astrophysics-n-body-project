# Makefile for compiling and running a C++ program

# Compiler and flags
CXX = g++
CXXFLAGS = -Wall

# Executable name
EXECUTABLE = nbody

# Source files
SOURCES = nbody.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(EXECUTABLE)
	./$(EXECUTABLE)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
