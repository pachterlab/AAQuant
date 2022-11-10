CXX = clang++
FLAGS = -O3
LFLAGS = -pthread -march=native -lbifrost -pthread -lz 
DEPS = AbundanceCDBG.hpp Node.hpp RemoveTips.hpp SparseVector.hpp Utils.hpp
OBJ = main.o AbundanceCDBG.o RemoveTips.o Node.o Utils.o

%.o: %.c ${DEPS}
	${CXX} -c -o $@ $< ${FLAGS}

AAQuant: ${OBJ}
	${CXX} -o $@ $^ ${FLAGS} ${LFLAGS}
