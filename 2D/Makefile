all: destar BA pmf 
OBJS2= funcion4.o main2.o
print:
	$(info Cargar con make/make all para correr el programa completo.)

clean:
	rm *txt *gz *dat HILL* *xtc *gro *pdb *png 

destar:
	 ./destar.sh
cols:
	mv COLVAR.0 COLS_1.txt
	mv COLVAR.1 COLS_2.txt
	mv COLVAR.2 COLS_3.txt
	mv COLVAR.3 COLS_4.txt
	mv COLVAR.4 COLS_5.txt
	mv COLVAR.5 COLS_6.txt
	gmx editconf -f *gro -o sistema.pdb

BA:
	python Block_avarage_2D.py 
pmf:
	python probabilidad.py



#CÁLCULO DEL ÁREA
area: $(OBJS2)
	gfortran -o area $(OBJS2)
	./area
main2.o: calculo_area.f90
	gfortran -o main2.o -c calculo_area.f90
#subrutina para obtener el número de carbonilos
funcion4.o: funcion_trapezoidal.f90
	gfortran -o funcion4.o -c funcion_trapezoidal.f90
