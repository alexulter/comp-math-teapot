run: part_alex.o part_jeka.o graphics.o main.o progonka_funct.o
	gcc part_alex.o part_jeka.o graphics.o main.o progonka_funct.o -o run `pkg-config --cflags --libs libglade-2.0` -Wall
progonka_funct.o: progonka_funct.c
	gcc progonka_funct.c -c -o progonka_funct.o -Wall
part_alex.o: part_alex.c
	gcc part_alex.c -c -o part_alex.o -Wall
part_jeka.o: part_jeka.c
	gcc part_jeka.c -c -o part_jeka.o -Wall
graphics.o: graphics.c
	gcc graphics.c -c -o graphics.o `pkg-config --cflags --libs libglade-2.0` -Wall
main.o: main.c
	gcc main.c -c -o main.o -Wall



