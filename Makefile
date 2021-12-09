TARGET_EXEC ?= a.out

BUILD_DIR ?= ./build
BUILD_DIR_SRC ?= ./build/src
SRC_DIRS ?= ./src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS += -g $(INC_FLAGS) -MMD -MP -fopenmp -O3 -Wall -Wextra -Wfloat-equal #-Werror
LDFLAGS += -lfftw3_mpi -lfftw3_threads -lfftw3 -lfftw3f_mpi -lfftw3f_threads -lfftw3f -lm
CXX := mpic++
CC := mpicc

all: $(BUILD_DIR)/$(TARGET_EXEC)
#Extra output in log file (to run QA tests)
QA: CPPFLAGS+=-DQA=1
QA: clean
QA: all
#Calculations made with single-precision
float: CPPFLAGS+=-DBYTES=4
float: clean
float: all
#Calculations made with double-precision
double: CPPFLAGS+=-DBYTES=8
double: clean
double: all

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(MKDIR_P) $(dir $@)output
	$(CXX) -fopenmp $(OBJS) -o $@ $(LDFLAGS)

# assembly
$(BUILD_DIR)/%.s.o: %.s
	$(MKDIR_P) $(dir $@)
	$(AS) $(ASFLAGS) -c $< -o $@

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)/$(TARGET_EXEC)
	$(RM) -r $(BUILD_DIR_SRC)

-include $(DEPS)

MKDIR_P ?= mkdir -p
