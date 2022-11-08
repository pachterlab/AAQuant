CXX = clang++
#HTS = htslib-5f5c1a557eb9d65f8b7854b6ca4c898713e86c04
#FLAGS=-fsanitize=address -g -O0
FLAGS=-g -O0
#FLAGS=-O3

all: utils node remove_tips abundance_cdbg main #bamreader

bamreader:
	${CXX} -std=c++11 -c ${FLAGS} -o BAMReader.o BAMReader.cpp
	#${CXX} -std=c++11 -c ${FLAGS} -I ${HTS} -o BAMReader.o BAMReader.cpp

abundance_cdbg:
	${CXX} -std=c++11 -march=native -c ${FLAGS} -o AbundanceCDBG.o AbundanceCDBG.cpp
	#${CXX} -std=c++11 -mavx -c ${FLAGS} -I ${HTS} -o AbundanceCDBG.o AbundanceCDBG.cpp

main:
	${CXX} -std=c++11 -march=native -c ${FLAGS} -o main.o main.cpp
	${CXX} -o AAQuant ${FLAGS} -pthread -march=native -lbifrost -pthread -lz main.o AbundanceCDBG.o RemoveTips.o Node.o Utils.o

node:
	${CXX} -std=c++11 -c ${FLAGS} -o Node.o Node.cpp

utils:
	${CXX} -std=c++11 -c ${FLAGS} -o Utils.o Utils.cpp

remove_tips:
	${CXX} -std=c++11 -c ${FLAGS} -o RemoveTips.o RemoveTips.cpp

sparse_vector_test:
	${CXX} -std=c++11 -c ${FLAGS} -o SparseVectorTest.o SparseVectorTest.cpp
	${CXX} -o sparse_vector_test ${FLAGS} -pthread -lz SparseVectorTest.o

.PHONY: clean
clean:
	@rm -f BAMReader.o Utils.o Node.o RemoveTips.o AbundanceCDBG.o main.o AAQuant
