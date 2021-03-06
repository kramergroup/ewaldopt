MODULE STRUCTURE

  IMPLICIT NONE

  TYPE :: COORDINATE
     real, dimension(3) :: tau              ! atomic coordinates in direct coordinates
     real, dimension(:), allocatable :: q   ! vector of nval valences
     integer :: nval                        ! number of valences
     character (len=5) :: symbol            ! atomic symbol
  END TYPE COORDINATE
CONTAINS

  ! input of structures for ewald summation from file
  ! It assumes that the unit has already been opened to read from
  SUBROUTINE READSTRUCTURE(iin,iout,t1,t2,t3,coords,nion)

    ! Arguments
    integer, intent(in) :: iin                                             ! Index of fortran unit to read from
    integer, intent(in) :: iout                                            ! Index of fortran unit to write output to
    real, intent(out) :: t1(3),t2(3),t3(3)                                 ! Lattice vectors
    type (coordinate), intent(out), allocatable, dimension(:) :: coords    ! Coordinates of atoms in direct coordinates
    integer, intent(out) :: nion                                           ! Number of ions

    ! Parameters
    real,parameter :: bohr = 0.529                                         ! Bohr radius  
    integer, parameter :: nmaxval = 10                                     ! Number of maximum valences

    ! Temporary and dummy variables
    integer :: i,n,ii,jj
    character*255 :: line,line2                                            ! the line buffer
    real,dimension(nmaxval) :: tval


    ! Read lattice vectors
    READ(iin, *) t1(1), t1(2), t1(3)
    READ(iin, *) t2(1), t2(2), t2(3)
    READ(iin, *) t3(1), t3(2), t3(3)

    ! Convert lattice vectors to Bohr radii
    DO i=1,3
       t1(i) = t1(i) / bohr
       t2(i) = t2(i) / bohr
       t3(i) = t3(i) / bohr
    END DO

    READ(iin, *) nion

    IF ( nion < 1 ) THEN
       WRITE(1, *) "No atomic coordinates given"
       STOP
    END IF

    ALLOCATE( coords(nion) )

    DO i=1,nion
       READ(iin, '(A)', end=999) line
       READ(line, *, end=888) (coords(i)%tau(ii), ii=1,3), coords(i)%symbol, (tval(ii), ii = 1, nmaxval)

777 CONTINUE
       ! This point is only reached if more valences are specified than we anticipate
       WRITE(1,*) "More than ", nmaxval, " valences not supported!"
       STOP
888 CONTINUE
       ALLOCATE( coords(i)%q(ii-1) )
       DO n=1,ii-1
          coords(i)%q(n) = tval(n)
       END DO
    END DO
999 CONTINUE

    WRITE(iout,*) "READSTRUCTURE------------------------------------------------------"
    WRITE(iout,"(A,E16.7,E16.7,E16.7)") "  Lattice vector A: " , t1(1), t1(2), t1(3)
    WRITE(iout,"(A,E16.7,E16.7,E16.7)") "  Lattice vector B: " , t2(1), t2(2), t2(3)
    WRITE(iout,"(A,E16.7,E16.7,E16.7)") "  Lattice vector C: " , t3(1), t3(2), t3(3)
    WRITE(iout,*)
    WRITE(iout,"(A,I4,A)") " Coordinates : (" , nion, " total)"
    WRITE(iout,*) "         X              Y               Z         Valences"
    DO i=1,nion
       WRITE(iout,"(' ',E16.7,E16.7,E16.7,A6,*(F7.2))") coords(i)%tau(1), coords(i)%tau(2), coords(i)%tau(3), coords(i)%symbol,coords(i)%q
    END DO
    WRITE(iout,*) "-------------------------------------------------------------------"

    RETURN

  END SUBROUTINE READSTRUCTURE
  
  

  SUBROUTINE WRITESTRUCTURE(out,t1,t2,t3,coords,valsel)
    
    ! Parameters
    integer, intent(in) :: out                                              ! Fortran unit for output
    real, intent(in) :: t1(3),t2(3),t3(3)                                   ! Lattice vectors
    type (coordinate), intent(in), allocatable, dimension(:) :: coords      ! Coordinates of atoms in direct coordinates
    integer, intent(in), allocatable, dimension(:) :: valsel                ! Selected valence set

    ! Parameters
    real,parameter :: bohr = 0.529                                         ! Bohr radius  

    ! Counter and local variables
    integer :: i
    character(len=8) :: i_char

    WRITE(out,"(E16.7,E16.7,E16.7)") t1(1)*bohr, t1(2)*bohr, t1(3)*bohr
    WRITE(out,"(E16.7,E16.7,E16.7)") t2(1)*bohr, t2(2)*bohr, t2(3)*bohr
    WRITE(out,"(E16.7,E16.7,E16.7)") t3(1)*bohr, t3(2)*bohr, t3(3)*bohr
 
    WRITE(i_char, '(I8)') size(coords)
    WRITE(out,*) adjustl(i_char) 

    DO i=1,size(coords)
       WRITE(out,"(E16.7,E16.7,E16.7,A5,F7.2)") coords(i)%tau(1), coords(i)%tau(2), coords(i)%tau(3), TRIM(coords(i)%symbol), coords(i)%q(valsel(i))
    END DO

  END SUBROUTINE WRITESTRUCTURE

END MODULE STRUCTURE
