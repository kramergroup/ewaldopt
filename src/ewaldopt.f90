! 
! Optimises the Ewald sum for variable valence lattices
! 
! The implementation uses an exhaustive search algorithm
! Beware of the scaling!!!!
!
PROGRAM EWALDOPT

  USE STRUCTURE
  USE EWALD

  IMPLICIT NONE
  
  ! Variable declarations
  real :: t1(3), t2(3), t3(3)                                   ! Lattice vectors
  integer  :: nion                                              ! Number of atoms 

  type (coordinate), allocatable, dimension(:) :: coords        ! Atomic coordinates (in Bohr radii)

  integer, allocatable, dimension(:) :: valsel                  ! Selected valence set
  integer, allocatable, dimension(:) :: optvalsel               ! Optimal valence set

  integer :: nvarval                                            ! Total number of variable valences

  real :: curewald, minewald                                    ! Ewald sums
  real :: charge, mincharge                                     ! Charge of a structure 
  
  logical, allocatable, dimension(:,:) :: symmtx                ! Matrix of symmetry equivalent coordinates 

  logical :: allowcharged = .false.
  logical :: forcebrute = .false.
 
  ! Dummy variables and counters, etc.
  integer :: i,j,k
  
  ! Command-line arguments
  character(len=100) :: arg,dummy
  integer :: narg

  ! Parameters
  real, parameter :: ZERO_TOLERANCE = 1e-6
  character(len=*), parameter :: VERSION = "1.0.0"

  ! I/O units
  integer :: null,stdin,stdout, log
  
  ! Initialize I/O units
  stdin = 5                                                     ! unit 5 is stdin
  stdout = 6                                                    ! unit 6 is stdout
  null = 7
  log = null
  open(unit=null, file='/dev/null')                             ! ignore output from subroutines


  ! Process command line arguments
  narg = IARGC()
  DO i=1,narg
     CALL GETARG(i,arg)

     IF (arg == "-h") THEN
        CALL PRINT_USAGE(stdout)
        STOP
     ELSE IF (arg == "-v") THEN
        CALL PRINT_VERSION(stdout)
        STOP
     ELSE IF (arg == "-c") THEN
        allowcharged = .true.
     ELSE IF (arg == "-l") THEN
        log=8
        open(unit=log, file='ewaldopt.log')                     ! log-file
     ELSE IF (arg == "-b") THEN
        forcebrute = .true.
     ELSE IF (arg == "-debug") THEN
        null = stdout
        ELSE IF (arg == "-m") THEN
           CALL GETARG(i+1, dummy)
           OPEN(10,FILE=dummy)
           READ(10,*) j
           ALLOCATE( symmtx(j,j) )
           DO k=1,j
              READ(10,*) symmtx(k,:) 
           END DO
           CLOSE(10)
     END IF
        
  END DO

  ! Read the structural information
  CALL READSTRUCTURE(stdin,null,t1,t2,t3,coords,nion)

  ! Prepare first valence set and calculate ewald sum
  ALLOCATE( valsel(nion) )
  ALLOCATE( optvalsel(nion) )

  ! Initialize symmetry matrix if not specified with P1
  IF (.NOT. allocated(symmtx) ) THEN
     ALLOCATE( symmtx(nion,nion) )
     symmtx = .false. 
     DO i=1,nion 
        symmtx(i,i) = .true.
     END DO
  END IF

  DO i=1,nion
     valsel(i) = 1
  END DO

  ! Count number of variable valences
  nvarval = 0
  DO i=1,nion
     IF ( (size(coords(i)%q) > 1) .AND. ( ALL( .NOT. symmtx(i,1:i-1) ) )) THEN
        nvarval = nvarval + size( coords(i)%q )
     END IF
  END DO
  
  WRITE(null,*) "Number of variable sites: :", nvarval

  ! First call to initialize variables
  call EWALDSUM(null,t1,t2,t3,coords,valsel,curewald,charge)
  minewald = curewald
  optvalsel = valsel
  mincharge = charge

  ! Only use brute-force optimisation for small cells,
  ! otherwise use a simple Monte-Carlo search
  IF ( nvarval < 10 .OR. forcebrute ) THEN
     WRITE(null,*) "Using brute-force algorithm"
     CALL BRUTEFORCEOPT
  ELSE 
     WRITE(null,*) "Using MC algorithm"
     CALL MCOPT
  END IF

  IF ( (.not. allowcharged) .and. (mincharge /= 0.0) ) THEN
     WRITE(1,*) " No uncharged cell found! Use option -c to allow charged cells"
     STOP -1
  END IF

  ! WRITE STRUCTURE
  CALL WRITESTRUCTURE(stdout,t1,t2,t3,coords,optvalsel)

CONTAINS  

  ! Brute-Force optimisation - Iteratates over all possible permutations
  ! to find minimum structure
  SUBROUTINE BRUTEFORCEOPT

    integer :: nred
    integer, dimension(size(coords)) :: cidx

    real :: chg 

    ! Find reduced coordinates (not symmetry equivalent coordinates)
    nred = 1
    cidx(1) = 1
    DO i=2,size(symmtx,1)
       IF ( ALL( .NOT. symmtx(i,1:i-1) ) ) THEN
          nred = nred+1
          cidx(nred) = i
       END IF
    END DO

    ! Now iterate over all permutations and find minimum energy
    i=1
    DO WHILE (i <= nred)
       IF (valsel(cidx(i)) < size(coords(cidx(i))%q)) THEN
          ! Apply new valence to all symmetry equivalent coordinates
          DO j=1,nion
             IF ( symmtx(cidx(i),j) ) valsel(j) = valsel(j)+1
          END DO
          DO WHILE (i > 1) 
             i=i-1
             ! Restore initial valence
             DO j=1,nion
                IF ( symmtx(cidx(i),j) ) valsel(j) = 1
             END DO
          END DO

          ! Catch charged states if needed
          IF ( .NOT. allowcharged ) THEN
             chg = 0.0
             DO j=1,nion
                chg = chg + coords(j)%q(valsel(j))
             END DO

             IF ( abs(chg) > ZERO_TOLERANCE ) THEN 
                WRITE(null, *) "Rejecting charged state"
                GOTO 300
             END IF
          END IF
          
          call EWALDSUM(null,t1,t2,t3,coords,valsel,curewald,charge)
          WRITE(log,"(F10.4,F10.4,*(I3))") curewald, charge, valsel
        
          IF ( allowcharged .OR. (charge == 0.0) ) THEN
             IF ( (curewald < minewald) .OR. (mincharge /= 0 .and. .not. allowcharged) ) THEN
                minewald = curewald
                optvalsel = valsel
                mincharge = charge
             END IF
          END IF

300       CONTINUE

       ELSE
          i=i+1
       END IF
    END DO
  END SUBROUTINE BRUTEFORCEOPT

  ! Simple Monte-Carlo algorithm to find low energy ewald assignment
  SUBROUTINE MCOPT
    integer :: nmaxtrials = 1E3
    real :: beta

    integer :: n, nvar, naccepted
    integer, dimension(:), allocatable :: idx
    
    real :: lastewald
    integer :: lastsel

    logical :: reject

    real,dimension(3) :: r

    beta = 10.0 ! Use a temperature that scales with the problem - Not sure if this makes sense
   
    ! Count number of variable valences
    nvar = 0
    DO i=1,nion
       IF ( size(coords(i)%q) > 1 ) THEN
          nvar = nvar + 1
       END IF
    END DO
    
    ! Allow 10 attempts per variable site
    nmaxtrials = nvar * 10

    ! Collect coordinates with variable valence
    ALLOCATE ( idx(nvar) )
    n = 0
    DO i=1,nion
       IF ( size(coords(i)%q) > 1 ) THEN
          n = n+1
          idx(n) = i
       END IF
    END DO
    
    CALL RANDOM_SEED
    
    lastewald = curewald

    naccepted = 0
    ! Main MC loop
    DO n=1,nmaxtrials

       CALL RANDOM_NUMBER(r)

       i = idx(int(r(1)*nvar)+1)                 ! pick a coordinate
       j = int(r(2)*size(coords(i)%q))+1         ! pick a valence
       
       ! Remember current setting
       lastsel = valsel(i)

       ! Apply new valence to all symmetry equivalent atoms
       DO k=1,size(coords)
          IF ( symmtx(i,k) ) valsel(k) = j
       END DO
       
       call EWALDSUM(null,t1,t2,t3,coords,valsel,curewald,charge)

       reject = (curewald > lastewald) .AND. (EXP(-(curewald-lastewald)/beta) > r(3))

       WRITE(null,*) n, i, j, minewald, mincharge, lastewald, curewald, charge, reject, r(3), EXP(-(curewald-lastewald)/beta)

       ! Keep track of minimal energy structure - total
       IF ( allowcharged .OR. (charge == 0.0) ) THEN
             IF ( (curewald < minewald) .OR. (mincharge /= 0 .and. .not. allowcharged) ) THEN
                 minewald = curewald
                 optvalsel = valsel
                 mincharge = charge
                 WRITE(null,*) "New global minimum: ", curewald 
             END IF
       END IF
       
       IF ( reject ) THEN
          ! Reject move
          DO k=1,size(coords)
             IF ( symmtx(i,k) ) valsel(k) = lastsel
          END DO
          WRITE(null,*) "REJECT", r(3)
       ELSE
          lastsel = valsel(i)
          lastewald = curewald
       END IF

    END DO

  END SUBROUTINE MCOPT

  SUBROUTINE PRINT_USAGE(out)
  
    integer, intent(in) :: out             ! fortran unit

    WRITE(out,*) "ewaldopt - Finds valence distribution in variable valence compounds"
    WRITE(out,*) "           with minimal ewald sum"
    WRITE(out,*) " Options:"
    WRITE(out,*) "  -c    Allow charged cells"
    WRITE(out,*) "  -l    Write logging info into file ewaldopt.log"
    WRITE(out,*) "  -h    Display this help message"
    WRITE(out,*) "  -d    Write debugging info to stdout"
    WRITE(out,*) "  -b    Use brute-force regardless of problem size. Beware: this can"
    WRITE(out,*) "        become very expensive very quickly"
    WRITE(out,*) "  -v    Display program version"
    RETURN
  
  END SUBROUTINE PRINT_USAGE
  
  SUBROUTINE PRINT_VERSION(out)

    integer, intent(in) :: out

    WRITE(out,*) "ewaldopt - Version " // VERSION

  END SUBROUTINE PRINT_VERSION
  
END PROGRAM
