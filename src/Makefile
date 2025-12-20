CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Wpedantic
DEBUGFLAGS = -std=c++17 -g -Wall -Wextra -Wpedantic

# Output directories
BUILD_DIR = build
BIN_DIR = bin

TARGET = $(BIN_DIR)/queue_simulator

# Source files
SOURCES = main.cpp \
          core/QueueSimulator.cpp \

# Header files
HEADERS = header/Event.hpp \
          header/QueueSimulator.hpp \

OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)

# Default target
all: build

# Build the project
build: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)
	@echo "Build complete!"

$(BUILD_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Iheader -c $< -o $@

# Run with default parameters
run: build
	./$(TARGET)

# Clean only object files
clean-obj:
	rm -rf $(BUILD_DIR)
	@echo "Cleaned object files"

# Clean only CSV output files
clean-output:
	rm -f *.csv
	@echo "Cleaned output CSV files"

# Clean build artifacts and output files
clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(BIN_DIR)
	rm -f *.csv
	@echo "Cleaned build artifacts and output files"

# Build with debug symbols
debug: CXXFLAGS = $(DEBUGFLAGS)
debug: clean $(TARGET)
	@echo "Debug build complete"
