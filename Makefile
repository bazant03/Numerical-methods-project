TARGET = main

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

SRC = main.cpp GaussianElimination.cpp
OBJ = $(SRC:.cpp=.o)
DEPS = GaussianElimination.h

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(OBJ) $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
