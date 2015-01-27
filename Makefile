NAME := HTseq_QA

#CC := g++
#CPP := g++

# For C++11 
CC := g++-4.9 -std=c++11 
CPP := g++-4.9 -std=c++11 

LINK := $(CPP)

MODULES := src include 

CFLAGS := -O3 -fopenmp
LFLAGS := -O3 
DEPFLAGS := -O3 

ifdef DEBUG
CFLAGS := -g 
LFLAGS := -g  
endif

LIBS := -lm -lz

# pcrelib linker -lpcre 
# zlib linker -lz 
# openmp linker -lgomp

# You shouldn't have to go below here

DIRNAME = `dirname $1`
MAKEDEPS = gcc -MM -MG $2 $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS := $(patsubst %, -I%, $(MODULES)) 

CFLAGS += $(INCLUDEFLAGS)
CPPFLAGS += $(INCLUDEFLAGS)
DEPFLAGS += $(INCLUDEFLAGS)

# each module will add to this
SRC :=  $(wildcard $(patsubst %, %/*.c, $(MODULES))) 

# determine the object files
OBJ :=   $(patsubst %.c, obj/%.o, $(filter %.c, $(SRC)))  

# link the program
$(NAME) : $(OBJ)
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

# calculate C include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

.dep/%.d : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(DEPFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

obj/%.o : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CC) $(CFLAGS) -c -o $@ $<

# include the C include dependencies
DEP := $(patsubst obj/%.o, .dep/%.d, $(OBJ))

ifneq ($(MAKECMDGOALS),clean)
include $(DEP)
endif

clean :
	-@rm $(NAME) $(OBJ) $(DEP)


