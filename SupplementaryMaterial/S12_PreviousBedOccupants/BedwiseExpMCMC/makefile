OBJS = main.o readcsv.o relabel.o logl.o createdata.o TPFN.o pacquisition.o samplenormal.o environment.o
CC = g++
FLAGS = -Wall
LIBS = -L/usr/include/boost

all : main

main : $(OBJS)
	$(CC) $(FLAGS) $(LIBS) $(OBJS) -o main.out

main.o : main.cpp stdhdr.h readcsvfileto2Darrayofint.h relabel.h loglikelihood.h createusefuldata.h truepositiveandfalsenegative.h samplefromnormal.h environment.h
	$(CC) $(FLAGS) $(LIBS) main.cpp

readcsv.o : readcsvfileto2Darrayofint.h readcsvfileto2Darrayofint.cpp
	$(CC) -c $(FLAGS) readcsvfileto2Darrayofint.cpp

relabel.o : relabel.h relabel.cpp
	$(CC) -c $(FLAGS) relabel.cpp

logl.o : loglikelihood.h loglikelihood.cpp stdhdr.h pacquisition.h
	$(CC) -c $(FLAGS) loglikelihood.cpp

createdata.o : createusefuldata.h createusefuldata.cpp
	$(CC) -c $(FLAGS) createusefuldata.cpp

TPFN.o : truepositiveandfalsenegative.h truepositiveandfalsenegative.cpp stdhdr.h
	$(CC) -c $(FLAGS) truepositiveandfalsenegative.cpp

pacquisition.o : pacquisition.h pacquisition.cpp stdhdr.h
	$(CC) -c $(FLAGS) pacquisition.cpp

samplenormal.o : samplefromnormal.h samplefromnormal.cpp
	$(CC) -c $(FLAGS) samplefromnormal.cpp

environment.o : environment.h environment.cpp stdhdr.h
	$(CC) -c $(FLAGS)  environment.cpp

clean: 
	\rm *. *~ main
