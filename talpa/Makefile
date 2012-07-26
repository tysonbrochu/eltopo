
# --------------------
#
# Makefile for Talpa, and example application using El Topo.
#
# Create a file called Makefile.local_defs before running make.  See Makefile.example_defs for suggestions.
# You must also build the El Topo library before building Talpa.
#
# --------------------

# This is for GNU make; other versions of make may not run correctly.

# Local machine settings
include Makefile.local_defs

# El Topo library
LIBRARY = libeltopo
include $(ELTOPO_PATH)/Makefile.inc
REL_LIB_SRC = $(LIB_SRC:%.cpp=$(ELTOPO_PATH)/%.cpp)

# --------------------
# Executable
# --------------------

MAIN_PROGRAM = talpa

SRC = iomesh.cpp lexer.cpp newparser.cpp geometryinit.cpp scriptinit.cpp 

# drivers
      
SRC += curlnoise/noise.cpp drivers/meancurvature.cpp drivers/faceoff.cpp drivers/normaldriver.cpp drivers/sisccurlnoisedriver.cpp
		 
# common		 
SRC += ../common/bfstream.cpp ../common/makelevelset3.cpp ../common/marching_tiles_hires.cpp \
../common/runstats.cpp 

# SRC += ../common/meshes/meshloader.cpp ../common/meshes/ObjLoader.cpp 
#SRC += ../common/gluvi.cpp 

SRC += main.cpp

# object files
RELEASE_OBJ = $(patsubst %.cpp,obj/%.o,$(notdir $(SRC)))
DEBUG_OBJ = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(SRC)))

# how to make the main target (debug mode, the default)
$(MAIN_PROGRAM): $(DEBUG_OBJ) $(LIBRARY)_debug.a
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS) $(LIBRARY)_debug.a

# how to make the main target (release mode)
$(MAIN_PROGRAM)_release: $(RELEASE_OBJ) $(LIBRARY)_release.a
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS) $(LIBRARY)_release.a


.PHONY: release
release: $(LIBRARY)_release.a $(MAIN_PROGRAM)_release

.PHONY: debug
debug: $(LIBRARY)_debug.a $(MAIN_PROGRAM)


# how to compile each file
.SUFFIXES:
obj/%.o:
	$(CC) -c $(RELEASE_FLAGS) $(INCLUDE_PATH) -o $@ $<
obj_debug/%.o:
	$(CC) -c $(DEBUG_FLAGS) $(INCLUDE_PATH) -o $@ $<

# cleaning up
.PHONY: clean
clean:
	-rm -f obj/*.o obj/depend $(MAIN_PROGRAM) obj_debug/*.o obj_debug/depend $(MAIN_PROGRAM)_release *core
	-rm -f $(LIBRARY)_debug.a $(LIBRARY)_release.a
	
# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir obj
	-rm -f obj/depend
	$(foreach srcfile,$(REL_LIB_SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> obj/depend;)
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> obj/depend;)
	-mkdir obj_debug
	-rm -f obj_debug/depend
	$(foreach srcfile,$(REL_LIB_SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> obj_debug/depend;)
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> obj_debug/depend;)

-include obj/depend
-include obj_debug/depend
