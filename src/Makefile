# Almost Universal Makefile to compile programs on linux with gcc

source_dirs:= . Log StringTools ConfigFile /home/dreamcatcher/Sources/levmar-2.6
search_wildcards:= $(addsuffix /*.cpp, $(source_dirs))
source_files:= $(wildcard $(search_wildcards))
source_files:=$(notdir $(source_files))
object_files:=$(patsubst %.cpp, %.o, $(source_files))
include_dirs:=$(addprefix -I, $(source_dirs)) 
release_flags:=  -O3 -Wall 
link_keys:= -lgsl  -lconfig++ -llevmar -lblas -llapack
link_dirs:= -L/home/dreamcatcher/Sources/levmar-2.6 
CC:=g++

DislScatCubFITMisfit: $(object_files)
	$(CC) $(link_dirs) $(object_files) -o $@ $(link_keys) 
	mv $@ ..
	rm *.o *.d
                         
clean:
	rm *.o *.d

VPATH:= $(source_dirs)

%.o: %.cpp
	$(CC) -c -MD $(release_flags)  $(include_dirs) $<

include $(wildcard *.d)
