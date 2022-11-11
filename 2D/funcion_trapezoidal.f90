double precision FUNCTION Trapezoidal(x,y,N,L)
double precision ::h,funci,xn,L,Itra
double precision, intent(in) :: x(N),y(N)
Itra=0
h=L/N
!print*,'xx',x
!print*,'yy',y
!print*,'NN',N
do i=0,N
	!xn=h*i
       !print*,x(i),y(i)
	if (i==0) then
		Itra=h*y(i)/2.+Itra
!		print*,'1'
	else if (i==N) then
		Itra=Itra+h*y(i)/2.
!		print*,'2'
	else
		Itra=Itra+h*y(i)
	end if
end do
Trapezoidal=Itra
End

