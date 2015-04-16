MODULE EWALD

  USE STRUCTURE
  IMPLICIT NONE

CONTAINS

  ! Ewaldsum fortran 90 code
  !   4/23/01  NAWH
  !   Input :   T1x  T1y   T1z    (Cartesian components of lattice vectors)
  !         :   T2x  T2y   T2z
  !         :   T3x  T3y   T3z
  !         :   nions             (Number of ions in unit cell)
  !         :   q(i), tau(1:3,i), i=1,nion  (charge and fractional position
  !                                             of ions)
  !             Note:  ion positions = tau(1,i)*T1 + tau(2,i)*T2 + tau(3,i)*T3
  !         :   eps    (error tolerance)

  SUBROUTINE EWALDSUM(iout,t1,t2,t3,coords,valselect,ewald,totalcharge)

    integer, intent(in) :: iout                                          ! Index of fortran unit to write output to
    real,intent(in) :: t1(3), t2(3), t3(3)                               ! Lattice vectors
    integer,intent(in), allocatable, dimension(:) :: valselect           ! Index of valence selection per atom
    type (coordinate), intent(in), allocatable, dimension(:) :: coords   ! Atomic coordinates
    
    real, intent(out) :: ewald                                           ! Ewald sum result
    real, intent(out) :: totalcharge                                     ! Total charge of the structure   

    ! Local Variables
    real :: pi,eps, g1(3),g2(3),g3(3),volcry, arg,x,gexp
    real :: eta, g1m, g2m, g3m, t1m, t2m, t3m,gcut, tmax ,ebsl,seta
    real :: tpi,glast2, con, con2 , cccc , v(3), w(3), rmag2 , prod
    integer  :: i, j, k, ng, nt , mmm1, mmm2, mmm3 , a , b, nion

    ! 
    ! Precision parameters for the summation
    !
    PARAMETER(gcut = 0.7, ebsl = 1E-12)

    pi = acos(-1.0)
    volcry   = t1(1)*(t2(2)*t3(3)-t2(3)*t3(2)) +  &
         t1(2)*(t2(3)*t3(1)-t2(1)*t3(3)) +  &
         t1(3)*(t2(1)*t3(2)-t2(2)*t3(1))

    g1(1) = 2 * pi * (t2(2)*t3(3)-t2(3)*t3(2))/volcry
    g1(2) = 2 * pi * (t2(3)*t3(1)-t2(1)*t3(3))/volcry
    g1(3) = 2 * pi * (t2(1)*t3(2)-t2(2)*t3(1))/volcry
    g2(1) = 2 * pi * (t3(2)*t1(3)-t3(3)*t1(2))/volcry
    g2(2) = 2 * pi * (t3(3)*t1(1)-t3(1)*t1(3))/volcry
    g2(3) = 2 * pi * (t3(1)*t1(2)-t3(2)*t1(1))/volcry
    g3(1) = 2 * pi * (t1(2)*t2(3)-t1(3)*t2(2))/volcry
    g3(2) = 2 * pi * (t1(3)*t2(1)-t1(1)*t2(3))/volcry
    g3(3) = 2 * pi * (t1(1)*t2(2)-t1(2)*t2(1))/volcry

    volcry = abs(volcry)

    nion = size(coords)

    t1m = SQRT(DOT_PRODUCT(t1,t1))
    t2m = SQRT(DOT_PRODUCT(t2,t2))
    t3m = SQRT(DOT_PRODUCT(t3,t3))
    g1m = SQRT(DOT_PRODUCT(g1,g1))
    g2m = SQRT(DOT_PRODUCT(g2,g2))
    g3m = SQRT(DOT_PRODUCT(g3,g3))

    tpi=2*pi
    con=volcry/(4*pi)
    con2=(4*pi)/volcry
    glast2=gcut*gcut
    gexp=-alog(ebsl)
    eta=glast2/gexp


    WRITE(iout,*) 'EWALDSUM-----------------------------------------------------------'
    WRITE(iout,*) ' Number atoms per cell: ', nion
    WRITE(iout,"(A,*(I3))") '  Valences: ', valselect 
    Write(iout,*) ' eta value for this calculation' , eta
    cccc=sqrt(eta/pi)

    x=0
    totalcharge=0
    do i=1,nion 
       x=x+coords(i)%q(valselect(i))**2
       totalcharge = totalcharge + coords(i)%q(valselect(i))
    enddo

    Write(iout,*) ' Total charge = ', totalcharge
    ewald=-cccc*x-4*pi*(totalcharge**2)/(volcry*eta)

    tmax=sqrt(2*gexp/eta)
    seta=sqrt(eta)/2

    mmm1=tmax/t1m+1.5
    mmm2=tmax/t2m+1.5    
    mmm3=tmax/t3m+1.5    


    Write (iout,"(A,I4,I4,I4)") '  Lattice summation indices: ', mmm1,mmm2,mmm3
    do a = 1,nion
       do b = 1,nion
          v(:) = (coords(a)%tau(1)-coords(b)%tau(1))*t1(:) &
               + (coords(a)%tau(2)-coords(b)%tau(2))*t2(:) &
               + (coords(a)%tau(3)-coords(b)%tau(3))*t3(:)
          prod=coords(a)%q(valselect(a))*coords(b)%q(valselect(b))
          do i = -mmm1, mmm1
             do j = -mmm2, mmm2
                do k = -mmm3, mmm3
                   if ((a.ne.b).or.((abs(i)+abs(j)+abs(k)).ne.0)) then
                      w(:) = v(:) + i*t1 + j*t2 + k*t3
                      rmag2 = sqrt(DOT_PRODUCT(w,w))
                      arg=rmag2*seta 
                      ewald = ewald + prod*erfc(arg)/rmag2
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    mmm1=gcut/g1m+1.5
    mmm2=gcut/g2m+1.5
    mmm3=gcut/g3m+1.5

    Write(iout,"(A,I4,I4,I4)") '  Reciprocal lattice summation indices:', mmm1,mmm2,mmm3
    do i = -mmm1, mmm1
       do j = -mmm2, mmm2
          do k = -mmm3, mmm3
             if ((abs(i)+abs(j)+abs(k)).ne.0) then
                w(:) = i*g1(:) + j*g2(:) + k*g3(:)
                rmag2=DOT_PRODUCT(w,w)
                x=con2*exp(-rmag2/eta)/rmag2
                do a = 1,nion
                   do b = 1,nion
                      v(:) = coords(a)%tau(:)-coords(b)%tau(:)
                      prod = coords(a)%q(valselect(a))*coords(b)%q(valselect(b))
                      arg=tpi*(i*v(1)+j*v(2)+k*v(3))
                      ewald=ewald + x*prod*cos(arg)
                   enddo
                enddo
             endif
          enddo
       enddo
    enddo

    Write(iout,*) ' Ewald energy [Ry]:', ewald
    WRITE(iout,*) '-------------------------------------------------------------------'
    RETURN
  END SUBROUTINE EWALDSUM

END MODULE EWALD
