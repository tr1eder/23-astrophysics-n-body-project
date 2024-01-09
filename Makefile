# Makefile for compiling and running a C++ program

# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O2

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
	./$(EXECUTABLE) $(filter-out $@,$(MAKECMDGOALS))

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

%:
	@:
