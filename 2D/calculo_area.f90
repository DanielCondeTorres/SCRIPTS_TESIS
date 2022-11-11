program INTEGRAL_REC
EXTERNAL funci

double precision, parameter :: boltzmann = 1.38064852d-26 ! KJ/K
double precision, parameter :: avogadro  = 6.02214076d+23 ! 1/mol
double precision, parameter :: pi        = 4.d0*datan(1.d0)
double precision, parameter :: R  = 8.314472/1000. ! kJ/(mol*K)
double precision, parameter :: T  = 298 !K
double precision, parameter :: Vstand=1661/10**3!nm^3
!Parametros caja
double precision, parameter :: x=5 ,y=5 !revisar valores... (nm)
double precision, parameter :: delta_W=-120 !kj/mol

double precision G_ESTAND,h,L1_nb,L1_b,L2,Trapezoidal,IT,col1,col2,AREA_CV1_nb,AREA_CV1_b,AREA_CV2
INTEGER::puntos,N
real*8,allocatable :: ycv1_nb(:),xcv1_nb(:),ycv2(:),xcv2(:),ycv1_b(:),xcv1_b(:)
integer :: Ngrid, i, Npoints, contador_cv1_nb, contador_cv1_b,contador_cv2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          CONTADOR DE LINEAS                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
puntos=100000000
contador_cv1_nb=0!contamos el número de líneas del  archivo
open(14,file="area_cv1_nonbonded.dat",status='old', action='read')!Abrimos el archivo
read(14,*)
DO i=1,puntos!Leemos línea por línea
        read(14,*,END=1114) col1,col2  !Apuntamos los valores de x e y
        contador_cv1_nb=contador_cv1_nb+1!Se ha leido una línea
end do
1114 close(14)!



contador_cv1_b=0!contamos el número de líneas del  archivo
open(15,file="area_cv1_bonded.dat",status='old', action='read')!Abrimos el archivo
read(15,*)
DO i=1,puntos!Leemos línea por línea
        read(15,*,END=1115) col1,col2  !Apuntamos los valores de x e y
        contador_cv1_b=contador_cv1_b+1!Se ha leido una línea
end do
1115 close(15)!
contador_cv2=0!contamos el número de líneas del  archivo
open(16,file="area_cv2.dat",status='old', action='read')!Abrimos el archivo
read(16,*)
DO i=1,puntos!Leemos línea por línea
        read(16,*,END=1116) col1,col2  !Apuntamos los valores de x e y
        contador_cv2=contador_cv2+1!Se ha leido una línea
end do
1116 close(16)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
puntos=100000000!Número de puntos inicial para leer archivo
Ngrid=puntos-1!Valor inicial del grid para almacenar datos
allocate(xcv1_b(contador_cv1_b),ycv1_b(contador_cv1_b))!creamos array para almacenar los datos
allocate(xcv1_nb(contador_cv1_nb),ycv1_nb(contador_cv1_nb))!creamos array para almacenar los datos
allocate(xcv2(contador_cv2),ycv2(contador_cv2))!creamos array para almacenar los datos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           CV1 non bonded                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(12,file="area_cv1_nonbonded.dat",status='old', action='read')!Abrimos el archivo
read(12,*)
DO i=1,puntos!Leemos línea por línea
        read(12,*,END=1111) col1,col2  !Apuntamos los valores de x e y
        xcv1_nb(i)=col1!coordenada x
        ycv1_nb(i)=col2!coordenada y
end do
1111 close(12)!cuando llegamos al EOF cerramos el programa
Ngrid=contador_cv1_nb!recalculamos el grid
L1_nb=MAXVAL(xcv1_nb)-MINVAL(xcv1_nb)!Longitud del intervalo para calcular el área
N=0!Mudo
k=1!Mudo
!Empleamos el método trapezoidal para obtener el cálculo del área.
do while (k<100000000)
        N=Ngrid
        IT=Trapezoidal(xcv1_nb,ycv1_nb,N,L1_nb)
k=k*2
print*,'###################################################################'
print*,'Número de puntos empleados',N
print*,'Área Trapezoidal',IT ,'kJ/mol*(unidad de CV1)'
END DO
AREA_CV1_nb=IT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          CV1  bonded                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(11,file="area_cv1_bonded.dat",status='old', action='read')!Abrimos el archivo
read(11,*)
DO i=1,puntos!Leemos línea por línea
        read(11,*,END=1110) col1,col2  !Apuntamos los valores de x e y
        xcv1_b(i)=col1!coordenada x
        ycv1_b(i)=col2!coordenada y
end do
1110 close(11)!cuando llegamos al EOF cerramos el programa
Ngrid=contador_cv1_b-1!recalculamos el grid
L1_b=MAXVAL(xcv1_b)-MINVAL(xcv1_b)!Longitud del intervalo para calcular el área
N=0!Mudo
k=1!Mudo
!Empleamos el método trapezoidal para obtener el cálculo del área.
do while (k<100000000)
        N=Ngrid
        IT=Trapezoidal(xcv1_b,ycv1_b,N,L1_b)
k=k*2
print*,'###################################################################'
print*,'Número de puntos empleados',N
print*,'Área Trapezoidal',IT ,'kJ/mol*(unidad de CV1)'
END DO
AREA_CV1_b=IT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                  CV2                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*,'-------------------------------------------------------------------'
open(13,file="area_cv2.dat",status='old', action='read')!Abrimos el archivo
read(13,*)
DO i=1,puntos!Leemos línea por línea
        read(13,*,END=1112) col1,col2  !Apuntamos los valores de x e y
        xcv2(i)=col1!coordenada x
        ycv2(i)=col2!coordenada y
end do
1112 close(13)!cuando llegamos al EOF cerramos el programa
Ngrid=contador_cv2-1!recalculamos el grid


L2=MAXVAL(xcv2)-MINVAL(xcv2)!Longitud del intervalo para calcular el área
N=0!Mudo
k=1!Mudo
!Empleamos el método trapezoidal para obtener el cálculo del área.
do while (k<100000000)
        N=Ngrid
        IT=Trapezoidal(xcv2,ycv2,N,L2)
k=k*2
print*,'####################################################################'
print*,'Número de puntos empleados',N
print*,'Área Trapezoidal',IT ,'kJ/mol*(unidad de CV2)'
END DO
AREA_CV2=IT
print*,'-------------------------------------------------------------------'
print*,'-----------------RESULTADOS FINALES ÁREAS--------------------------'
print*,'-------------------------------------------------------------------'
print*,'Área CV1',AREA_CV1_nb , 'unidad de CV1 non bonded'!'kJ/mol*(unidad de CV1)'
print*,'Área CV1',AREA_CV1_b , 'unidad de CV1 bonded'!'kJ/mol*(unidad de CV1)'
print*,'Área CV2',AREA_CV2 , 'unidad de CV2' !'kJ/mol*(unidad de CV2)'
print*,'-------------------------------------------------------------------'
print*,'-------------------------------------------------------------------'
print*,' '
print*,'-------------------------------------------------------------------'
print*,'--------------------RESULTADOS DE INTERES--------------------------'
print*,'-------------------------------------------------------------------'
G_ESTAND=Delta_W-R*T*LOG(AREA_CV1_b/L1_b)
print'(E10.4)',(AREA_CV1_b*L1_nb/(AREA_CV1_nb*L1_b))
print*, 'B: ',L1_b
print*, 'NB: ',L1_nb
print*,'-------------------------------------------------------------------'
print*,'-------------------------------------------------------------------'


end program








