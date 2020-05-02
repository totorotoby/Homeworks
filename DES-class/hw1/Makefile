mm1: mm1.o lcgrand.o lcgrand.h
	gcc -g -o mm1 mm1.o lcgrand.o -lm

tandum_queue: tandum_queue.o lcgrand.o lcgrand.h tandum_queue.h
	gcc -g -o tandum_queue tandum_queue.o lcgrand.o -lm


mm1.o: mm1.c lcgrand.h
	gcc -g -c mm1.c -lm

tandum_queue.o: tandum_queue.c lcgrand.h tandum_queue.h
	gcc -g -c tandum_queue.c -lm


lcgrand.o: lcgrand.c lcgrand.h
	gcc -g -c lcgrand.c


clean:
	rm *.o
