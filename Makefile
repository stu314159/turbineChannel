MPI_CC=CC
ifeq ($(PE_ENV),PGI)
	MPI_FLAGS=-O3 -fast -acc -Minfo=acc -Mnoopenmp
	ifeq ($(CRAYPAT_COMPILER_OPTIONS),1)
		MPI_FLAGS+= -DCRAYPAT
	endif
else
	MPI_FLAGS=-O3 -hnoomp -hacc -hlist=m
endif



SOURCES= turbineSim.cpp TurbineChannel3D.cpp
OBJECTS=turbineSim.o TurbineChannel3D.o  workArounds.o
LIBS=

ifeq ($(USE_NVTX),1)
	LIBS+=-lnvToolsExt
	MPI_FLAGS+= -DUSE_NVTX
endif

TARGET=turbineSim

%.o: %.cpp
	$(MPI_CC) $(MPI_FLAGS) -c $^

$(TARGET): $(OBJECTS)
	$(MPI_CC) $(MPI_FLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(TARGET) *~
