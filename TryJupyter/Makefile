.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET	=  BC

C_SOURCES =  BlumeCapel.cpp 
C_OBJS     =  BlumeCapel.o
MY_INCLUDES = 


CCX = g++
CXXFLAGS = -g -Wall

#============================================================
all: $(TARGET)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

NOTES = 
%= otherstuff.np 

clean:
	rm -f $(TARGET) $(C_OBJS) core *~*

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES)  $(NOTES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


