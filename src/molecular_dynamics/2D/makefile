# Commet lines
# Here we define compiler option, libraries and the target
FLAGS = -O2 -Wall

# Here we make the executable file
SRCS = main.cpp atom.cpp io.cpp lennardjones.cpp statisticssampler.cpp system.cpp unitconverter.cpp vec2.cpp velocityverlet.cpp 
OBJS = $(subst .cpp,.o,$(SRCS))
all: md

# Whereas here we create the object file
md: $(OBJS) 
	g++ ${FLAGS} -o md $(OBJS)
#	export OMP_NUM_THREADS=20
	./md

%.o: %.cpp  
	g++ $(FLAGS) -c $<

#	g++ ${FLAGS} -c init.cpp
#	g++ ${FLAGS} -c integrator.cpp
#	g++ ${FLAGS} -c MD.cpp

# Clean
clean:
	rm *.o ./md
