!> \file
!! Contains modules mod_chem, subroutine mod_chem::init_chem() \n
!! mod_chem::get_eq_chemistry
!! initialize the physical constants and units
!<
!> Hold subroutines which solve for local chemical equilibrium,  \n
!! currently Gibbs minimisation and analytical formula of Burrows and Sharp 1999
!<
MODULE mod_chem_neq

USE mod_chem

!Network parameters
INTEGER :: nfile  !< number of reaction files
INTEGER :: nreac  !< number of total reactions
INTEGER :: nphoto !< number of photolytic reactions

INTEGER,ALLOCATABLE :: methodParam(:) !< parameter for which method to calculate rate
INTEGER,ALLOCATABLE :: reverseParam(:)!< parameter for whether to reverse rates

!Thermochemistry
INTEGER,PARAMETER :: nreactants = 5 !< max number of reactants
INTEGER,PARAMETER :: nproducts = 5  !< max number of products
INTEGER,PARAMETER :: nko = 5        !< number of parameters for kO
INTEGER,PARAMETER :: nkinf = 5      !< number of parameters for kinf
INTEGER,PARAMETER :: ntroe = 4      !< number of parameters for troe
INTEGER,PARAMETER :: nsri = 5       !< number of parameters for SRI
INTEGER,PARAMETER :: neff = 22      !< number of efficiency species
INTEGER,PARAMETER :: nmod = 2       !< number of parameters for modParam
!! my variable - oct. 31
INTEGER,PARAMETER :: npressures = 8       !< max number of pressures for pressure dependent files
INTEGER,PARAMETER :: ncoeffs = 26       !< coeffs for pressure dependent files

INTEGER :: nmol_neq                   !< number of molecules to be included in chemical kinetics
INTEGER :: nit_neq                    !< iterations of the solver
INTEGER :: nit_max                    !< maximum allowed iterations
INTEGER :: ndepth_neq                 !< depth at which to solve chemistry
INTEGER :: neff_network               !< actual number of efficiency molecules in network
INTEGER,ALLOCATABLE :: imol_reac(:,:) !< ids of the reactants
INTEGER,ALLOCATABLE :: imol_prod(:,:) !< ids of the products
INTEGER,ALLOCATABLE :: nb_reac(:)     !< number of reactants in each reaction
INTEGER,ALLOCATABLE :: nb_prod(:)     !< number of products in each reaction
INTEGER,ALLOCATABLE :: modParam(:,:)  !< contains file id of reactions
INTEGER,ALLOCATABLE :: inloc(:,:)     !< map of the dndt matrix
INTEGER,ALLOCATABLE :: eff_list(:)    !< Indices of molecules included as efficienies
INTEGER,ALLOCATABLE :: eff_list_mol(:)

REAL             :: tfreeze        !< minimum temperature for the reaction rates
REAL             :: pfreeze        !< do not calculate chemistry below pfreeze
REAL             :: unitpfreeze    !< units of pfreeze

!! my variables - oct. 31
REAL,ALLOCATABLE :: pressures_v20(:,:) !< pressures for pressure dependent files
REAL,ALLOCATABLE :: ko_v20(:,:) !< array to store the pressure dependet file coefficients

REAL,ALLOCATABLE :: koParam(:,:)   !< low pressure rate parameter
REAL,ALLOCATABLE :: ko2Param(:,:)  !< additional low pressure rate parameter (reaction method 6)
REAL,ALLOCATABLE :: kinfParam(:,:) !< high pressure rate parameter
REAL,ALLOCATABLE :: troeParam(:,:) !< troe parameters
REAL,ALLOCATABLE :: effParam(:,:)  !< efficiency parameters, third body reactions
REAL,ALLOCATABLE :: sriParam(:,:)  !< SRI formalism parameters
REAL,ALLOCATABLE :: kf(:,:)        !< forward reaction rate
REAL,ALLOCATABLE :: kr(:,:)        !< reverse reaction rate
REAL,ALLOCATABLE :: rate_lim(:,:)  !< rate limiter
REAL             :: dt_lim         !< limit on timestep
REAL             :: nn_lim         !< limit on density
REAL             :: tt_lim         !< limit on temperature
REAL             :: depth_lim      !< depth at which to apply rate limits

CHARACTER(lname),ALLOCATABLE :: reactants(:,:) !< reactant names
CHARACTER(lname),ALLOCATABLE :: products(:,:)  !< product names
CHARACTER(lname) :: effmol(neff)               !< names of molecules which act as third bodies

! Convergence criteria
REAL :: dnn_tol    !< tolerance for dn/n
REAL :: dnndt_tol  !< tolerance for (dn/n)/dt
REAL :: mf_tol     !< mole fraction limit on which to test dn/d and (dn/n)/dt
REAL :: tmin       !< minimum integration time
REAL :: ttest      !< time interval on which to test convergence


!Photochemistry
INTEGER             :: nqy              !< number of quantum yield files
INTEGER             :: nse              !< number of cross section files
INTEGER             :: nuv              !< number of UV wavelength points
INTEGER             :: mod_fuv          !< frequency on which to update uvflux
INTEGER,ALLOCATABLE :: qyflag(:)        !< flag for branching ratio of reactions
REAL,ALLOCATABLE    :: se(:,:)          !< cross sections
REAL,ALLOCATABLE    :: qy(:,:)          !< branching ratios
REAL,ALLOCATABLE    :: fuv_nu(:,:)      !< uvflux
REAL,ALLOCATABLE    :: selambda(:,:)    !< wavelength arrays of cross section files
REAL,ALLOCATABLE    :: qylambda(:,:)    !< wavelength arrays of branching ratio files
REAL,ALLOCATABLE    :: Jrate(:)         !< photochemistry rates
REAL,ALLOCATABLE    :: Irad_uv(:)       !< uv irradiation of star
REAL,ALLOCATABLE    :: nu_uv(:)         !< wavelengths of uv irradiation
REAL,ALLOCATABLE    :: snu_scat_uv(:,:) !< source function
CHARACTER(60)       :: fhnu_uv          !< name of the input file for uv flux hnu
LOGICAL             :: photochem        !< flag to include photochemistry

INTEGER :: mod_pt      !< frequency at which to recalculate PT profile
INTEGER :: mod_jac     !< frequency at which to recalculate Jacobian
INTEGER :: mod_rate    !< frequency at which to recalculate rates

REAL    :: Nmin        !< minimum allowed number density
REAL    :: tmax        !< maximum integration time
REAL    :: dtmin       !< minimum allowed timestep
REAL    :: dtmax       !< maximum allowed timestep
REAL    :: dt          !< timestep
REAL    :: varrel      !< allowed relative variations of abundances
REAL    :: vardt       !< allowed variation of timestep

CHARACTER(30) :: dir_neq  !< name of the directory containing the non-eq network

LOGICAL :: rate_limiter !< flag to turn on rate limiter
LOGICAL :: check_lbound !< flag to check if lower boundary remains in eq


CONTAINS

  !----------------------
  !Initialise chemical kinetics
  !----------------------
  SUBROUTINE init_neq_chem

  USE mod_param
  USE mod_grid
  USE mod_cst
  USE mod_util

  IMPLICIT NONE

  INTEGER :: i,stat

  INCLUDE 'netcdf.inc'

  NAMELIST /chem_neq/ nfile, nmol_neq, dir_neq, rate_limiter, tfreeze, pfreeze, photochem, &
  fhnu_uv, nuv, mod_fuv, tmax, Nmin, &
  dtmin, dtmax, mod_rate, mod_jac, mod_pt, dt_lim, nn_lim, tt_lim, depth_lim, check_lbound, &
  dnn_tol, dnndt_tol, mf_tol, tmin, nit_max, ttest

  !Initialise namelist variables to default
  !network
  nfile        = 13
  nmol_neq     = 107
  dir_neq      = '../../chem/venot2012/'
  tfreeze      = 10.
  pfreeze      = 1.0E6
  unitpfreeze  = 1.0E6
  rate_limiter = .False.
  dt_lim       = 1.0E7
  nn_lim       = 1.0E-30
  tt_lim       = 100.
  depth_lim    = 0.9

  ! Convergence criteria
  dnn_tol   = 1e-2
  dnndt_tol = 1e-4
  mf_tol    = 1e-30
  tmin      = 1e5
  nit_max   = int(1e6)
  ttest     = 1.1

  !photochemistry
  photochem = .False.
  fhnu_uv   = 'None'
  nuv       = 900
  mod_fuv   = 50

  mod_rate = 1
  mod_jac  = 10000
  mod_pt   = 100

  tmax         = 1.0E12
  dtmin        = 1.0E-10
  dtmax        = 1.0E10
  Nmin         = 1.0E-100

  varrel       = 0.2
  vardt        = 0.2
  dt           = dtmin
  check_lbound = .false.

  ! read grid namelist in the file fparam
  OPEN(1,file=fparam)
  REWIND(1) ; READ(1,NML=chem_neq,iostat=stat)
   IF (stat/=0) THEN
    BACKSPACE(1)
    READ(1,fmt='(A)') line
    WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
    STOP
   END IF
  CLOSE(1)

  !Read in the chemical network
  !-Thermochemical reactions
  !-Photochemical reactions
  !-UV irradiation
  CALL read_network

  !Initialise variables
  ALLOCATE(kf(ndepth+1,nreac))
  ALLOCATE(kr(ndepth+1,nreac))
  kf = 0.0 ; kr = 0.0

  pfreeze = pfreeze*unitpfreeze/unitP
  ndepth_neq = ndepth+1
  !set ndepth_neq
  DO i=1,ndepth+1
    IF (ppf(i)<pfreeze) ndepth_neq=i
  END DO

  IF (mype == cputerm) THEN
    WRITE(*,headrfmt) 'Chemical Kinetics'
    WRITE(*,strpmfmt) 'dir_neq',  'Network directory'//dir_neq
    WRITE(*,intpmfmt) 'No. of molcules', nmol_neq
    WRITE(*,intpmfmt) 'No. of reactions', nreac
    WRITE(*,exppmfmt) 'Pfreeze',pfreeze*unitP/1.0E6,'[bar] pressure above which abondances are frozen'
    WRITE(*,logpmfmt) 'Mixing',mixing
    WRITE(*,logpmfmt) 'Photochemistry',photochem
    IF (mixing) THEN
      WRITE(*,exppmfmt) 'kzzcst [cm2 s-1]',kzzcst*unitL**2/unitT,'if = 0 read kzz from input profile'
    END IF
    IF (photochem) THEN
      WRITE(*,strpmfmt) 'UV irad file',fhnu_uv
    END IF
    WRITE(*,exppmfmt) 'dnn_tol',dnn_tol, 'Convergence condition for dn/n'
    WRITE(*,exppmfmt) 'dnndt_tol',dnndt_tol, 'Convergence condition for (dn/n)/dt'
    WRITE(*,exppmfmt) 'mf_tol',mf_tol, 'Mole fraction limit on which to assess convergence'
    WRITE(*,exppmfmt) 'tmin',tmin,'Minimum integration time'
    WRITE(*,exppmfmt) 'tmax',tmax,'Maximum integration time'
    WRITE(*,intpmfmt) 'nit_max',nit_max,'Maximum allowed timesteps'
    WRITE(*,fltpmfmt) 'ttest',ttest,'Time interval to test convergence'
    WRITE(*,*)
    WRITE(*,separfmt)
  END IF

  END SUBROUTINE init_neq_chem

  !----------------------
  !Get non-equilibrium chemistry
  !----------------------
  SUBROUTINE get_neq_chemistry

  USE mod_grid
  USE mod_eos
  USE mod_cst
  USE mod_param
  USE mod_chem
  USE mod_solver

  IMPLICIT NONE

  EXTERNAL dndt, jac

  INTEGER :: i, idepth,idepth_max, im, im_max, conv_test
  INTEGER :: nbstep, nbstep_old, node, ndepth_max, itime
  REAL    :: time, dnn, dnndt, told, tvold
  REAL,ALLOCATABLE :: N_loc(:)
  REAL,ALLOCATABLE :: Nold(:), Nvold(:)
  REAL             :: Amolfeq(ndepth+1,nmol_neq)
  LOGICAL :: l_continue

  CHARACTER(5)  :: fnit

  !time in seconds and variations
  time     = 0.
  told     = 0.
  tvold    = 0.

  !number of iterations
  nit_neq  = 0

  nbstep = 0
  nbstep_old = 0

  !convergence test
  l_continue = .TRUE.

  !!Setup atmosphere
  !Get number densities
  nn  = get_nn(ndepth,pp,tt)
  nnf = get_nn(ndepth+1,ppf,ttf)

  !Get mass densities
  rho = get_rho(ndepth,pp,tt,amol)
  rhof = get_rho(ndepth+1,ppf,ttf,amolf)

  !Get radius and gravity of grid
  CALL get_rrgr
  ! Get thickness of layers
  CALL get_drdrf
  dr = dr*unitL
  drf = drf*unitL

  !!Calculate initial reaction rates
  CALL get_krkf

  ndepth_max = ndepth_neq
  !if ndepth=1, means the chemistry is done on a box model
  IF (ndepth==1) ndepth_max = 1
  node = nmol_neq*ndepth_max

  !Allocate arrays
  ALLOCATE(N_loc(node),Nold(node),Nvold(node))
  ALLOCATE(inloc(ndepth_max,nmol_neq))
  ALLOCATE(rate_lim(ndepth_max,nreac))

  N_loc = 0.

  !Map the molecules in inloc and calculate number densities of species
  DO idepth = 1,ndepth_max
    DO i=1,nmol_neq
      inloc(idepth,i) = (idepth-1)*nmol_neq+i
      N_loc(inloc(idepth,i)) = Amolf(idepth,i)*nnf(idepth)
    END DO
  END DO

  ! Initialise Nold arrays
  Nold = N_loc
  Nvold = N_loc

  !Copy equilibrium abundances
  Amolfeq(:,1:nmol_neq) = Amolf(:,1:nmol_neq)

  !Chemistry is done in CGS
  N_loc = N_loc*one_over_l3
  kzz   = kzz*unitL**2/unitT

  !Calculate rate limits
  IF (rate_limiter) THEN
    CALL get_rate_limits(node,n_loc)
  END IF

  ! Initialise LSODE solver
  CALL lsode_init(node)

  !Start iterations
  DO WHILE (l_continue)

    !Save previous number densities
    IF (time>=told*ttest) THEN
      tvold = told
      told  = time
      Nvold = Nold
      Nold  = N_loc
    END IF

    lsode_itask = 1

    !Set minimum number density
    IF (Nmin /= 0.) THEN
      WHERE (N_loc<Nmin) N_loc = 0.
    END IF

    !Call DLSODES solver
    CALL solve_dlsodes(node,N_loc,time,time+dt,dndt,jac)

    IF (Nmin /= 0.) THEN
      WHERE(N_loc<1E-100) N_loc = 0.
    END IF

    nbstep_old = nbstep
    nbstep = iwork(11)

    nit_neq  = nit_neq  + 1

    !If requested, recalculate rate coefficients
    IF (mod(nit_neq,mod_rate)==0 ) THEN
      CALL get_krkf
    END IF

    !If requested, reconverge PT profile
    IF (mod(nit_neq,mod_pt)==0 ) THEN
      IF (solve_hydro .or. solve_energy) THEN
        CALL find_profile('cst')
        CALL get_krkf
        lsode_istate =1
      END IF
    END IF

    !If solver error or requested, recalculate on next iteration
    IF (lsode_istate<0 .OR. nbstep-nbstep_old>lsode_nit .OR. mod(nit_neq,mod_jac)==0) THEN
      lsode_istate = 1
      CALL get_krkf
    END IF

    IF (lsode_istate==1) THEN
      nbstep_old = 0
      nbstep = 0
    END IF

    !Copy number density to abundance array - check for negative number densities
    DO i = 1,node
      idepth = (i-1)/nmol_neq+1
      im     = i-(idepth-1)*nmol_neq
      Amolf(idepth,im) = N_loc(i)
      IF (N_loc(i)<0.) WRITE(*,*) 'Negative abundance: ',idepth,molname(im)
    END DO

    ! Convert to mole fraction - note for molecules not in network (imol>nmol_neq)
    ! mole fraction is not updated
    DO idepth = ndepth_max,1,-1
      Amolf(idepth,1:nmol_neq) = Amolf(idepth,1:nmol_neq)/ sum(Amolf(idepth,1:nmol_neq))
    END DO

    ! Find maximum variations in number densities
    ! Loop over species and levels
    dnn   = 0.
    dnndt = 0.
    DO i = 1, node
      IF (N_loc(i)>0. .AND. Nvold(i)>0.) THEN
        idepth = (i-1)/nmol_neq+1
        im     = i-(idepth-1)*nmol_neq
        ! Only check for convergence for molecules with Amolf>mf_tol
        IF (Amolf(idepth,im)>mf_tol) THEN
          ! Find max variations
          dnn   = max(abs(N_loc(i)-Nvold(i))/N_loc(i),dnn)
          dnndt = max(abs(N_loc(i)-Nvold(i))/N_loc(i)/(time-tvold),dnndt)
        END IF
      END IF
    END DO

    ! Check for convergence
    !   Reset counter if criteria not met on this timestep
    IF (dnn>dnn_tol .OR. dnndt>dnndt_tol)  conv_test = 0
    !   If criteria met add one to counter
    IF (time>tmin .AND. dnn<dnn_tol .AND. dnndt<dnndt_tol) conv_test = conv_test + 1
    !   If convergence criteria met for three consecutive time steps stop iterations
    !   or maximum time reached or maximum allowed iterations reached
    IF (conv_test==3 .OR. time>tmax .OR. nit_neq>nit_max) l_continue = .FALSE.


    !Calculate new timestep
    IF (lsode_dt) THEN
      dt = 0.9*rwork(11)
    ELSE
      IF (dnn<varrel) THEN
        dt=vardt*dt
      ELSE
        dt=(2.-vardt)*dt
      END IF
    END IF

    dt = min(dtmax,dt)
    dt = max(dtmin,dt)

    DO idepth = ndepth_max,1,-1
      DO i=1,nmol
        !Set abundance to zero if below Tcrit/above Pcrit
        IF (ttf(idepth)<Tcrit(i) .or. ppf(idepth)>Pcrit(i)) THEN
          Amolf(idepth,i) = 0.
        END IF
        !Set abundance to Acst
        IF (Acst(i)>0.) THEN
          Amolf(idepth,i) = Acst(i)
        END IF
      END DO
    END DO

    !Calculate abundance at cell centre; geometric mean
    DO idepth = 1,ndepth
      Amol(idepth,:) = sqrt(Amolf(idepth,:)*Amolf(idepth+1,:))
    END DO

    ! Check if still in equilibrium at lower boundary - required by assumption
    ! of zero mass flux through lower boundary
    ! Note: currently only works if P-T profile is constant (changing T will also cause
    ! changed in Amolf)
    ! Check using N2 and NH3
    IF (check_lbound) THEN
      IF ((abs(1.-(Amolfeq(ndepth+1,imol('N2'))/Amolf(ndepth+1,imol('N2')))) > 1.0e-2) &
      .OR. (abs(1.-(Amolfeq(ndepth+1,imol('NH3'))/Amolf(ndepth+1,imol('NH3')))) > 1.0e-2)) THEN
        WRITE(*,headrfmt) 'Error get_neq_chem'
        WRITE(*,'(A)') 'Lower boundary not in equilibrium'
        WRITE(*,'(A)') 'mixing scheme assumes zero flux at lower boundary'
        WRITE(*,'(A)') 'which is only valid if lower boundary remains in eq.'
        WRITE(*,*)
        WRITE(*,separfmt)
        STOP
      END IF
    END IF

    !update the total densities
    nn  = get_nn(ndepth,pp,tt)
    nnf = get_nn(ndepth+1,ppf,ttf)

    rho = get_rho(ndepth,pp,tt,amol)
    rhof = get_rho(ndepth+1,ppf,ttf,amolf)

    !Get radius and gravity of grid
    CALL get_rrgr
    ! Get thickness of layers
    CALL get_drdrf
    dr = dr*unitL
    drf = drf*unitL

    !print to screen
    IF (print_chem .and. mype ==cputerm) THEN
      WRITE(fnit,'(i5)') nit_neq
      WRITE(*,headrfmt) 'Chem. Iter.: '//fnit
      WRITE(*,999) time,dt,ppf(1)*unitP/1.0E6,ttf(1)*unitK
      WRITE(*,1000) dnn,dnndt,minval(n_loc),nbstep-nbstep_old
      WRITE(*,1001) molname(imol('CH4')),Amolf(1,imol('CH4')),molname(imol('H2O')),Amolf(1,imol('H2O')), &
      molname(imol('CO2')),Amolf(1,imol('CO2'))
      WRITE(*,*)
    END IF

999  FORMAT(5x,' time: ',es10.1,5x,' dt: ',es10.1,' Pout [bar]: ',es10.1,' Tout [K]: ',es10.1)
9990 FORMAT(5x,'ppf(1): ', es10.1)
1000 FORMAT(5x,' dn/n: ',es10.1,3x,' (dn/n)/dt: ',es10.1,3x,' min(n): ',es11.1,3x,' nbstep: ',i4)
1001 FORMAT(5x,a10,': ',es10.2,', ',a10,': ',es10.2,', ',a10,': ',es10.2)

    !Write to file
    CALL write_chemistry


  END DO

  IF (solve_hydro .or. solve_energy) THEN
    CALL find_profile('cst')
  END IF

  IF (print_chem .and. mype ==cputerm) THEN
    WRITE(*,separfmt)
    WRITE(*,headrfmt) 'Convergence Crit. met'
    WRITE(*,1002) time, dnn, dnndt, nit_neq
    WRITE(*,separfmt)
  END IF
1002 FORMAT(5x,' time: ',es10.1,3x,' dn/n: ',es10.1,3x,' (dn/n)/dt: ',es10.1,3x, ' iter: ',i5)

  END SUBROUTINE get_neq_chemistry

  !----------------------
  !Function to get wavelength index
  !----------------------
  FUNCTION ilam(ilrequest)

  USE mod_param

  IMPLICIT NONE

  INTEGER :: i,ilam
  REAL :: ilrequest

  ilam = 0

  DO i = 1,nuv
    IF (ilrequest == nu_uv(i)) THEN
      ilam = i
    END IF
  END DO

  IF (ilam==0) THEN
    IF (mype==cputerm) THEN
      WRITE(*,headrfmt) 'Error ilam'
      WRITE(*,'(A)') 'Requested wavelength not found'
      WRITE(*,chrpmfmt) 'Requested lambda.',ilrequest
      WRITE(*,*)
      WRITE(*,separfmt)
      STOP
    END IF
  END IF

  END FUNCTION ilam

	!--------------------------------------------------------
	! Function to calculate eddy and molecular diffusion flux
	! Requires molecule mass, diameter, abundance and dr
	!--------------------------------------------------------
	FUNCTION get_phi(j,i,node,nmol,inloc,n_loc)

	USE mod_chem, ONLY: molmass, dmol
	USE mod_grid, ONLY: ndepth, rhof, tt, pp, kzz, dr, kzz
	USE mod_cst, ONLY: kb, pi, amu, unitT, unitL

	IMPLICIT NONE

	INTEGER,INTENT(IN) :: i, j
	INTEGER,INTENT(IN) :: node
	INTEGER,INTENT(IN) :: nmol
	INTEGER,INTENT(IN) :: inloc(ndepth+1,nmol)
	REAL,INTENT(IN)    :: n_loc(node)

	INTEGER :: idepth
	REAL     :: ntot_centre, kzz_centre, da, phi, diffmol
	REAL     :: ntot(ndepth+1)
	REAL     :: get_phi(ndepth)


	!Calculate ntot on cell faces
	DO idepth = 1, ndepth+1
		ntot(idepth) = rhof(idepth)
	END DO

	!Calculate phi at cell centre
	DO idepth = 1,ndepth

		!Calculate total number density at cell centre
		ntot_centre = (ntot(idepth)+ntot(idepth+1))/2.

		!Calculate eddy diffusion coeff at cell centre, convert to code units
		kzz_centre = (kzz(idepth)+kzz(idepth+1))/2.

		!Calculate diffusion coefficient at cell centre
		diffmol = (1./3.)*sqrt((3.*kb*tt(idepth))/(molmass(i)*amu)) &
		*(1./(sqrt(2.)*pi*dmol(i)**2.))*((kb*tt(idepth))/pp(idepth))

		!Convert diffmol to cgs
		diffmol = diffmol*unitL**2./unitT

		!calculate abundance gradient
		da = ((n_loc(inloc(idepth+1,j))/ntot(idepth+1))- &
		(n_loc(inloc(idepth,j))/ntot(idepth)))/dr(idepth)

		!calculate molecular diffusion flux
		phi = -diffmol*da*ntot_centre

		!add eddy diffusion flux
		phi = phi+(-kzz_centre*da*ntot_centre)

		get_phi(idepth) = phi

	END DO

	END FUNCTION get_phi


END MODULE mod_chem_neq

!----------------------
!Read in the chemical network
!----------------------
SUBROUTINE read_network

USE mod_grid
USE mod_chem
USE mod_chem_neq
USE mod_param
USE mod_util

IMPLICIT NONE

INTEGER :: stat,id_file,id_var
INTEGER :: i,j,iphoto,ireac
REAL    :: qyfile,sefile,il
CHARACTER(2) :: rd,fifile
CHARACTER(12),ALLOCATABLE :: qy_name(:),se_name(:)

INCLUDE 'netcdf.inc'

!Count the number of reactions
nreac = 0

DO i = 1,nfile
  WRITE(fifile,'(i2)') i
  OPEN(1,file=trim(dir_neq)//'reactions_'//trim(adjustl(fifile))//'.dat')

  DO
    READ(1,'(a)',IOSTAT=stat) rd
    IF (stat<0) EXIT
    nreac = nreac + 1
  END DO
  CLOSE(1)
END DO

!If photochemistry
IF (photochem) THEN
  nphoto = 0

  !Add number of photolytic reactions
  OPEN(1,file=trim(dir_neq)//'photodissociations.dat')
  DO
    READ(1,'(a)',IOSTAT=stat) rd
    IF (stat<0) EXIT
    nreac = nreac + 1
    nphoto = nphoto + 1
  END DO
  CLOSE(1)

  !Count number of branching ratio and uv cross section files
  nqy=0
  OPEN(1,file=trim(dir_neq)//'filenames_qy.txt')
  DO
    READ(1,'(a2)',IOSTAT=stat) rd
    IF (stat<0) EXIT
    nqy = nqy + 1
  END DO
  CLOSE(1)

  nse=0
  OPEN(1,file=trim(dir_neq)//'filenames_se.txt')
  DO
    READ(1,'(a2)',IOSTAT=stat) rd
    IF (stat<0) EXIT
    nse = nse + 1
  END DO
  CLOSE(1)

END IF

!Read in network parameter file
ALLOCATE(methodParam(nfile))
ALLOCATE(reverseParam(nfile))

OPEN(1,file=trim(dir_neq)//'network_param.dat')

READ(1,*)
DO i = 1, nfile
  READ(1,'(11x,i10,i10)') methodParam(i),reverseParam(i)
END DO
CLOSE(1)

!Allocate arrays-Thermochemical
ALLOCATE(reactants(nreac,nreactants))
ALLOCATE(products(nreac,nproducts))
ALLOCATE(koParam(nreac,nko))
ALLOCATE(ko2Param(nreac,nko))
ALLOCATE(kinfParam(nreac,nkinf))
ALLOCATE(troeParam(nreac,ntroe))
ALLOCATE(effParam(nreac,neff))
ALLOCATE(sriParam(nreac,nsri))
ALLOCATE(modParam(nreac,nmod))
ALLOCATE(imol_reac(nreac,nreactants))
ALLOCATE(imol_prod(nreac,nproducts))
ALLOCATE(nb_reac(nreac))
ALLOCATE(nb_prod(nreac))
!! my variables - oct. 31
ALLOCATE(pressures_v20(nreac,npressures))
ALLOCATE(ko_v20(nreac,ncoeffs))


! initialize to 0
koParam  = 0. ; kinfParam = 0. ; troeParam = 0. ; sriParam = 0.
ko2Param = 0.
effParam = 0. ; modParam  = 0
imol_reac = 0 ; imol_prod = 0
nb_reac   = 0 ; nb_prod   = 0
neff_network = 0
!! my variables to 0
pressures_v20 = 0
ko_v20 = 0

!Read in the thermochemical reactions
ireac = 1
DO i = 1, nfile
  WRITE(fifile,'(i2)') i
  OPEN(1,file=trim(dir_neq)//'reactions_'//trim(adjustl(fifile))//'.dat')

  DO
    !Read in format depends on the method of calculation for reactions in file
    IF     (methodParam(i) == 1) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 5(1x,e10.3))',END=1) reactants(ireac,:),products(ireac,:),koParam(ireac,:)
    ELSE IF (methodParam(i) == 2) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 5(1x,e10.3) 22(1x,e10.3))',END=1) reactants(ireac,:),products(ireac,:), &
      koParam(ireac,:),effParam(ireac,:)
    ELSE IF (methodParam(i) == 3) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 10(1x,e10.3) 4(1x,e10.3) 22(1x,e10.3))',END=1) reactants(ireac,:), &
      products(ireac,:),koParam(ireac,:),kinfParam(ireac,:),troeParam(ireac,:),effParam(ireac,:)
    ELSE IF (methodParam(i) == 4) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 10(1x,e10.3) 5(1x,e10.3) 22(1x,e10.3))',END=1) &
      reactants(ireac,:),products(ireac,:),koParam(ireac,:),kinfParam(ireac,:),sriParam(ireac,:),effParam(ireac,:)
    ELSE IF (methodParam(i) == 5) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 10(1x,e10.3) 4(1x,e10.3) 22(1x,e10.3))',END=1) reactants(ireac,:), &
      products(ireac,:),koParam(ireac,:),kinfParam(ireac,:),troeParam(ireac,:),effParam(ireac,:)
    ELSE IF (methodParam(i) == 6) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 15(1x,e10.3) 5(1x,e10.3) 22(1x,e10.3))',END=1) &
      reactants(ireac,:),products(ireac,:),koParam(ireac,:),ko2Param(ireac,:),kinfParam(ireac,:),sriParam(ireac,:),effParam(ireac,:)
      ! Pressure dependent reactions
    ELSE IF (methodParam(i) == 7) THEN
      READ(1,'(5(1x,a10) 5(1x,a10) 8(1x,e10.3) 26(1x,e10.3))',END=1) &
      reactants(ireac,:),products(ireac,:),pressures_v20(ireac,:),ko_v20(ireac,:)

      ! converting pressures_v20 from atmospheres to bars
      pressures_v20(ireac,:) = pressures_v20(ireac,:) * 1.01325
      
    END IF
    


    modParam(ireac,1) = i
    modParam(ireac,2) = reverseParam(i)

    ireac = ireac + 1
  END DO
1 CLOSE(1)
END DO

!Molecules for efficiencies in three-body reactions
effmol = (/'O2   ', 'CO   ', 'CO2  ', 'H2O  ', 'CH4  ', 'H2   ', 'C2H6 ', 'Ar   ', 'N2   ', &
&'He   ', 'C2H4 ', 'cC6H6', 'C7H8 ', 'H2O2 ', 'N2O  ', 'O-3P ', 'NH3  ', 'N2H4 ', 'N2O4 ', 'NO2  ', 'NO   ', 'N-4S '/)

! Build list of indices for efficiency molecules
DO i = 1, neff
  IF (ismol(effmol(i))) neff_network = neff_network + 1
END DO

ALLOCATE(eff_list(neff_network))
ALLOCATE(eff_list_mol(neff_network))
eff_list = 0 ; eff_list_mol = 0

j = 0
DO i = 1, neff
  IF (ismol(effmol(i))) THEN
    j = j + 1
    eff_list(j) = i
    eff_list_mol(j) = imol(effmol(i))
  END IF
END DO

!if photochemistry included, read in photodissociations, uv cross sections and branching ratios
IF (photochem) THEN

  !Include photochemistry file

  !Read in UV irradiation file
  IF (fhnu_uv == 'None') THEN

    IF (mype == cputerm) THEN
      WRITE(*,headrfmt) 'Error init_neq_chem'
      WRITE(*,'(A)') 'No UV irradiation file provided'
      WRITE(*,chrpmfmt) 'fhnu_uv',fhnu_uv
      WRITE(*,'(A)') 'Provide netcdf file with uv spectrum'
      WRITE(*,*)
      WRITE(*,separfmt)
      STOP
    END IF

  ELSE

    CALL nf(nf_open(fhnu_uv,nf_nowrite,id_file))

    CALL nf(nf_inq_dimid (id_file,'nlambda',id_var))
    !Get number of UV wavelengths
    CALL nf(nf_inq_dimlen(id_file,id_var,nuv))

    ALLOCATE(Irad_uv(nuv))
    ALLOCATE(nu_uv(nuv))

    CALL nf(nf_inq_varid(id_file,'Hnu' ,id_var))
    CALL nf(nf_get_vara_double(id_file,id_var,(/1/),(/nuv/),Irad_uv(:)))
    CALL nf(nf_inq_varid(id_file,'Wavelength' ,id_var))
    CALL nf(nf_get_vara_double(id_file,id_var,(/1/),(/nuv/),nu_uv(:)))

    CALL nf(nf_close(id_file))

  END IF

  !allocate arrays- Photodissociations
  ALLOCATE(qyflag(nphoto))
  ALLOCATE(qy_name(nqy),se_name(nse))
  ALLOCATE(se(nse,nuv))
  ALLOCATE(qy(nqy,nuv))
  ALLOCATE(Jrate(nphoto))
  ALLOCATE(snu_scat_uv(ndepth+1,nuv))
  ALLOCATE(fuv_nu(ndepth+1,nuv))

  se = 0.0 ; qy = 0.0 ;
  Jrate = 0.0 ; snu_scat_uv = 0.0 ; fuv_nu = 0.0

  OPEN(1,file=trim(dir_neq)//'photodissociations.dat')
2005  FORMAT(5(1x,a10) 5(1x,a10) (1x) (i1))

  DO iphoto = 1, nphoto

    READ(1,2005,IOSTAT=stat) reactants(ireac,:),products(ireac,:),qyflag(iphoto)

    IF (stat<0) THEN
      IF (mype==0) THEN
        WRITE(*,headrfmt) 'Error chem_neq'
        WRITE(*,'(A)') 'Reached end of photochemistry reactions'
        WRITE(*,intpmfmt) 'Expected',nphoto, ' reactions'
        WRITE(*,intpmfmt) 'Got',iphoto, ' reactions'
        WRITE(*,separfmt)
        STOP
      END IF
    END IF

    modParam(ireac,1) = nfile !Assuming only one file photodissociations!
    modParam(ireac,2) = 1 ! irreversible reactions

    ireac = ireac + 1

  END DO
  CLOSE(1)

  !Read in branching ratio filenames
  OPEN(1,file=trim(dir_neq)//'filenames_qy.txt')
  DO i=1,nqy
    READ(1,'(a13)') qy_name(i)
  END DO
  CLOSE(1)

  !Read in cross sections file names
  OPEN(1,file=trim(dir_neq)//'filenames_se.txt')
  DO i = 1, nse
    READ(1,'(a13)') se_name(i)
  END DO
  CLOSE(1)

  !Read in cross sections
  DO i = 1, nse
    j = 0
    OPEN(1,file=trim(dir_neq)//trim(se_name(i))//'.dat')
    DO
      j = j + 1
      IF (j .le. nuv) THEN
        !selambda is in nm and the cross section 'se' in cm**2
         READ(1,'(1x,f7.1,x,e11.5)',IOSTAT=stat) il,sefile
         se(i,ilam(il)) = sefile
      ELSE
        READ(1,*,IOSTAT=stat)
      END IF
      IF (stat<0) EXIT
    END DO
    CLOSE(1)
  END DO


  !Read in branching ratios
  DO i = 1, nqy
    j = 0
    OPEN(1,file=trim(dir_neq)//trim(qy_name(i))//'.dat')
    DO
      j = j + 1
      IF (j .le. nuv) THEN
        READ(1,'(1x,f7.1,x,e11.5)',IOSTAT=stat) il,qyfile
        qy(i,ilam(il)) = qyfile
      ELSE
        READ(1,*,IOSTAT=stat)
      END IF
      IF (stat<0) EXIT
    END DO
    CLOSE(1)
  END DO

END IF !end if photochemistry

!Map reactants and products
!Count number of products and reactants in each reaction
DO ireac = 1,nreac
  DO j = 1,nreactants
    i = imol(reactants(ireac,j))
    IF (i>0) THEN
      nb_reac(ireac) = nb_reac(ireac) + 1
      imol_reac(ireac,j) = i
    END IF
  END DO

  DO j = 1,nproducts
    i = imol(products(ireac,j))
    IF (i>0) THEN
      nb_prod(ireac) = nb_prod(ireac) + 1
      imol_prod(ireac,j) = i
    END IF
  END DO
END DO

END SUBROUTINE read_network

!----------------------
!Calculate forward and reverse reaction rates
!----------------------
SUBROUTINE get_krkf

USE mod_chem
USE mod_grid
USE mod_cst
USE mod_chem_neq

IMPLICIT NONE

INTEGER :: i
INTEGER :: ireac, idepth, ipr, exp_pt, iloc
REAL    :: ko_loc, ko_loc_v20, kinf_loc, tloc, ceff_loc, zfc, exp_zfc, c_zfc, n_zfc, d_zfc &
           , K_loc, nloc, nref, p_loc
REAL    :: one_over_t, t_over_300
REAL    :: mu(nmol),h(nmol),s(nmol),cp(nmol)
REAL(kind=16) :: k_loc16

INTEGER :: pos !< pos in the array of rate coefficients
REAL :: alp !< alpha in the array of rate coefficients
REAL :: bta !< beta in the array of rate coefficients
REAL :: gma !< gamma in the array o rate coefficients
INTEGER :: zeros

!Calculate forward reaction rate
DO idepth= 1,ndepth+1

  !number density
  nloc = nnf(idepth)*one_over_l3

  !Set upper/lower limit lower limit on T
  IF ((ttf(idepth) >= tfreeze) .OR. (ttf(idepth) <= 6000.)) THEN
    tloc = ttf(idepth)*unitK
  ELSE IF (ttf(idepth) < tfreeze) THEN
    tloc = tfreeze*unitK
  ELSE
    tloc = 6000.
  END IF

  one_over_t = 1./tloc
  t_over_300 = tloc/300.

  zeros = 0

  ! Calculate thermodynamic quantities for tloc
  ! chemical potential required to calculate equilibrium constant
  CALL get_thermodynamics(nmol,tloc,molParam,mu,h,s,cp)

  DO ireac = 1,nreac-nphoto
  ! local presusre, converting to correct units
  p_loc = ppf(idepth)*unitP/1E06

    ! v20 pressure dependence, calculate ko_loc with pressure dependence
    IF (methodParam(modParam(ireac,1)) ==  7) THEN
      IF (p_loc <= pressures_v20(ireac,1)) THEN
        alp = ko_v20(ireac,1)
        bta = ko_v20(ireac,9)
        gma = ko_v20(ireac,17)
        pos = 1
        zeros = 1
       ELSE IF (p_loc > pressures_v20(ireac,1) .AND. p_loc <= pressures_v20(ireac,2)) THEN 
           pos = 2
       ELSE IF (p_loc > pressures_v20(ireac,2) .AND. p_loc <= pressures_v20(ireac,3)) THEN 
           pos = 3
       ELSE IF (p_loc > pressures_v20(ireac,3) .AND. p_loc <= pressures_v20(ireac,4)) THEN 
           pos = 4
       ELSE IF (p_loc > pressures_v20(ireac,4) .AND. p_loc <= pressures_v20(ireac,5)) THEN 
           pos = 5
       ELSE IF (p_loc > pressures_v20(ireac,5)) THEN
           IF (pressures_v20(ireac,6) == 0) THEN
               zeros = 1
               alp = ko_v20(ireac,5)
               bta = ko_v20(ireac,13)
               gma = ko_v20(ireac,21)
           ELSE IF (p_loc <= pressures_v20(ireac,6)) THEN
               pos = 6 
               zeros = 0
           END IF
            
       ELSE IF (p_loc > pressures_v20(ireac,6) .AND. (zeros == 0)) THEN 
           IF (pressures_v20(ireac,7) == 0) THEN
               zeros = 1
               alp = ko_v20(ireac,6)
               bta = ko_v20(ireac,14)
               gma = ko_v20(ireac,22)
           ELSE IF (p_loc <= pressures_v20(ireac,7)) THEN
               pos = 7
               zeros = 0
           END IF

       ELSE IF (p_loc > pressures_v20(ireac,7) .AND. (zeros == 0)) THEN 
           IF (pressures_v20(ireac,8) == 0) THEN
               zeros = 1
               alp = ko_v20(ireac,7)
               bta = ko_v20(ireac,15)
               gma = ko_v20(ireac,23)
           ELSE IF (p_loc <= pressures_v20(ireac,8)) THEN
               pos = 8
               zeros = 0
           END IF
            
       ELSE IF (p_loc > pressures_v20(ireac,8) .AND. (zeros == 0)) THEN 
           zeros = 1
           alp = ko_v20(ireac,8)
           bta = ko_v20(ireac,16)
           gma = ko_v20(ireac,24)

      END IF

      IF (zeros == 0) THEN
          ! Interpolating pressure dependent coefficients
          !a alpha
          alp = (ko_v20(ireac,pos) - ko_v20(ireac,pos-1))/(log(pressures_v20(ireac,pos)) - log(pressures_v20(ireac,pos-1))) * log(p_loc)

          ! beta
          bta  = (ko_v20(ireac,pos+8) - ko_v20(ireac,pos+7))/(log(pressures_v20(ireac,pos)) - log(pressures_v20(ireac,pos-1))) * log(p_loc)

          ! gamma
          gma = (ko_v20(ireac,pos+16) - ko_v20(ireac,pos+15))/(log(pressures_v20(ireac,pos)) - log(pressures_v20(ireac,pos-1))) * log(p_loc)
      END IF

      ko_loc_v20 = alp*(t_over_300)**bta*exp(-gma*one_over_t)

    ELSE 
      !Calculate rate coefficient in low pressure limit  or default way of calculating 
      ko_loc = koParam(ireac,1)*(t_over_300)**koParam(ireac,2)*exp(-koParam(ireac,3)*one_over_t)
    END IF

    !Second rate coefficient (reaction method 6)
    IF (methodParam(modParam(ireac,1)) == 6) THEN
      ko_loc = ko_loc + &
        ko2Param(ireac,1)*(t_over_300)**ko2Param(ireac,2)*exp(-ko2Param(ireac,3)*one_over_t)
    END IF

    kinf_loc = kinfParam(ireac,1)*(t_over_300)**kinfParam(ireac,2)*exp(-kinfParam(ireac,3)*one_over_t)

    !Calculate efficiencies for molecules defined in effmol
    IF (methodParam(modParam(ireac,1)) == 2 .OR. methodParam(modParam(ireac,1)) == 3 .OR. methodParam(modParam(ireac,1)) == 4) THEN
      ceff_loc = 1.
      DO i = 1,neff_network
        ceff_loc = ceff_loc - Amolf(idepth,eff_list_mol(i)) + effParam(ireac,eff_list(i))*Amolf(idepth,eff_list_mol(i))
      END DO
    ELSE IF (methodParam(modParam(ireac,1)) == 5) THEN
      ceff_loc = 0.
      DO i = 1,neff_network
        ceff_loc = ceff_loc + effParam(ireac,eff_list(i))*Amolf(idepth,eff_list_mol(i))
      END DO
    END IF

    !kooij formalism
    IF (methodParam(modParam(ireac,1)) == 1) THEN
      kf(idepth,ireac) = ko_loc

    !kooij formalism for v20 pressure dependent reactions
    ELSE IF (methodParam(modParam(ireac,1)) == 7) THEN
      kf(idepth,ireac) = ko_loc_v20

    !kooij formalism with third body
    ELSE IF (methodParam(modParam(ireac,1)) == 2) THEN
      kf(idepth,ireac) = ko_loc*ceff_loc*nloc

    !troe formalism
    ELSE IF (methodParam(modParam(ireac,1)) == 3 .OR. methodParam(modParam(ireac,1)) == 5) THEN
      ko_loc = ko_loc*ceff_loc*nloc
      IF (ko_loc==0.) THEN
        kf(idepth,ireac) = kinf_loc
      ELSE IF (kinf_loc == 0.) THEN
        kf(idepth,ireac) = ko_loc
      ELSE
        zfc = exp(-troeParam(ireac,4)/tloc)
        IF (.not.(troeParam(ireac,2)==0.)) THEN
          zfc = zfc + (1.-troeParam(ireac,1))*exp(-tloc/troeParam(ireac,2))
        END IF
        IF (.not.(troeParam(ireac,3)==0.)) THEN
          zfc = zfc + (troeParam(ireac,1))*exp(-tloc/troeParam(ireac,3))
        END IF

        c_zfc =  -0.40 - 0.67*log10(zfc)
        n_zfc =   0.75 - 1.27*log10(zfc)
        d_zfc =   0.14

        exp_zfc = 1.0/(1.0+((log10(ko_loc/kinf_loc)+c_zfc)/ &
        (n_zfc-d_zfc*(log10(ko_loc/kinf_loc)+c_zfc)))**2)

        kf(idepth,ireac) = (ko_loc/(1.0+ko_loc/kinf_loc))*zfc**exp_zfc
      END IF

    !SRI formalism
    ELSE IF (methodParam(modParam(ireac,1)) == 4) THEN
      ko_loc = ko_loc*ceff_loc*nloc
      IF (ko_loc==0.) THEN
        kf(idepth,ireac) = kinf_loc
      ELSE IF (kinf_loc == 0.) THEN
        kf(idepth,ireac) = ko_loc
      ELSE
        exp_zfc = 1./(1.+log10(ko_loc/kinf_loc)**2)
        IF (sriParam(ireac,3)==0.) THEN
          zfc = sriParam(ireac,4)*(sriParam(ireac,1)*exp(-sriParam(ireac,2)*one_over_t))**exp_zfc &
          *tloc**sriParam(ireac,5)
        ELSE
          zfc = sriParam(ireac,4)*(sriParam(ireac,1)*exp(-sriParam(ireac,2)*one_over_t) &
          +exp(-tloc/sriParam(ireac,3)))**exp_zfc*tloc**sriParam(ireac,5)
        END IF
        kf(idepth,ireac) = (ko_loc/(1.0+ko_loc/kinf_loc)) * zfc
      END IF

    END IF

    !Calculate reverse reaction rate
    IF (modParam(ireac,2) == -1) THEN

      K_loc = 0.
      exp_pt = 0

      DO i = 1,nb_prod(ireac)
        iloc = imol_prod(ireac,i)
        ipr = 1
        exp_pt = exp_pt + 1

        ! Add chemical potential for each product
        K_loc = K_loc+ipr*mu(iloc)

      END DO

      DO i = 1,nb_reac(ireac)
        iloc = imol_reac(ireac,i)
        ipr = -1
        exp_pt = exp_pt - 1

       ! Add chemical potential for each reactant
        K_loc = K_loc+ipr*mu(iloc)

      END DO

      !Calculate equilibrium constant
      nref = P0*one_over_t/kb*one_over_l3
      K_loc16 = K_loc
      K_loc16 = exp(-K_loc16)*nref**exp_pt

      !Calculate reverse rate
      kr(idepth,ireac) = REAL(kf(idepth,ireac)/K_loc16)

    END IF

  END DO !End loop on reactions

  !Set reactions to zero if requested
  DO ireac = 1,nreac
    IF (modParam(ireac,2) == 0) THEN
      kf(idepth,ireac) = 0.0
      kr(idepth,ireac) = 0.0
    END IF
  END DO

  !Calculate photodissociation rate
  IF (photochem) THEN
    IF (idepth==1) THEN
      IF (mod(nit_neq,mod_fuv)==0 ) THEN
        !Calculate UV flux
        CALL get_uv_flux
      END IF
    END IF

    !Calculate photodissociation rate
    CALL get_photo_rate(idepth)

    !Put photodissociation rates into kf array; no reverse reactions
    i = 1
    DO ireac = (nreac-nphoto+1),nreac
      kf(idepth,ireac) = Jrate(i)
      i = i + 1
    END DO
  END IF

END DO !end loop over vertical levels

END SUBROUTINE get_krkf

!----------------------
!Calculate the photodissociation rate
!----------------------
SUBROUTINE get_photo_rate(idepth)

USE mod_chem
USE mod_chem_neq
USE mod_grid
USE mod_cst

IMPLICIT NONE

INTEGER,INTENT(IN) :: idepth

!Local variables
REAL    :: rate(nuv)
INTEGER :: iphoto, iqy, ise, ilambda

Jrate = 0.
ise = 1
iqy = 1

!J rate = cross section * branching ratio * flux
!Must calculate at consistent wavelengths; qy,se arrays have different wavelength grids
!Some reactions have no qy files. These are signified with qyflag
DO iphoto = 1, nphoto

  ! Calculate spectral dissociation rate
  IF (qyflag(iphoto) /= 0) THEN
    rate(:) = se(ise,:)*qy(iqy,:)*abs(fuv_nu(idepth,:))
  ELSE
    rate(:) = se(ise,:)*abs(fuv_nu(idepth,:))
  END IF

  ! Integrate over wavelength range using trapezoidal rule
  DO ilambda = 1, nuv-1
    Jrate(iphoto) = Jrate(iphoto) + 0.5*(nu_uv(ilambda+1)-nu_uv(ilambda))*(rate(ilambda+1)+rate(ilambda))
  END DO

  ! Cap rates
  IF (Jrate(iphoto)<1E-80) Jrate(iphoto) = 0.
  IF (Jrate(iphoto)>1E+10) Jrate(iphoto) = 1E+10

 ! Set ise and iqy for next reaction
  IF (iphoto /= nphoto) THEN
    IF (qyflag(iphoto+1) <= 1) ise = ise + 1  !move on to next cross section
    IF (qyflag(iphoto+1) > 0)  iqy = iqy + 1  !move on to next branching ratio
  END IF

END DO

END SUBROUTINE get_photo_rate

!----------------------
!Get the UV flux
!----------------------
SUBROUTINE get_uv_flux

USE mod_grid
USE mod_radtrans
USE mod_chem
USE mod_eos
USE mod_cst
USE mod_chem_neq

IMPLICIT NONE

INTEGER :: idepth,ise,ilambda,iray
INTEGER :: nrays_uv
INTEGER :: semol(nse)

REAL :: alphap(num_ray_scatt),sig(num_ray_scatt)
REAL :: fstar_uv(ndepth+1,nuv)
REAL :: snu_uv(ndepth+1,nuv)
REAL :: bnu_uv(ndepth+1,nuv)
REAL :: hnup_uv(ndepth+1,nuv)
REAL :: jnup_uv(ndepth+1,nuv)
REAL :: hnud_uv(ndepth+1,nuv)
REAL :: jnud_uv(ndepth+1,nuv)
REAL :: kapfnu_uv(ndepth+1,nuv)
REAL :: kscatfnu_uv(ndepth+1,nuv)
REAL :: epsfnu_uv(ndepth+1,nuv)
REAL :: kapnu_uv(ndepth,nuv)
REAL :: kscatnu_uv(ndepth,nuv)
REAL :: epsnu_uv(ndepth,nuv)
REAL :: kapstd_uv(ndepth)
REAL :: mu(ndepth)
REAL :: muf(ndepth+1)
REAL :: Irad_uv_loc(nuv)
REAL :: nu_loc(nuv)

LOGICAL :: scatter_uv

CHARACTER(2) :: slname
CHARACTER(lname) :: semolname

!Molecules with UV cross sections
!Note, must match order of cross section files
OPEN(1,file=trim(dir_neq)//'filenames_se.txt')
DO ise = 1, nse
  WRITE(slname,'(i2)') lname
  READ(1,'(2x,a'//slname//')') semolname
  semol(ise) = imol(trim(semolname))
END DO
CLOSE(1)

!Get mean molecule weight
mu  = get_mu(ndepth,amol)
muf = get_mu(ndepth+1,amolf)

!Initialise rad_trans variables
nrays_uv = 16
scatter_uv = .true.

fstar_uv = 0. ; snu_uv = 0. ; bnu_uv = 0. ; hnup_uv = 0. ; jnup_uv = 0.
hnud_uv = 0. ; jnud_uv = 0.
Irad_uv_loc = 0. ; kapstd_uv = 0.
kapnu_uv  = 0. ; kscatnu_uv  = 0. ; epsnu_uv  = 1.
kapfnu_uv = 0. ; kscatfnu_uv = 0. ; epsfnu_uv = 1.

!Set irradiation to UV irradiation
Irad_uv_loc = Irad_uv*1.0E7/(unitL**3./unitT) !Convert into CGS and code units
nu_loc      = nu_uv*1.0E-7/unitL              !Convert into CGS and code units

!calculate standard opacity from hydrostatic equilibrium
DO idepth = 1,ndepth
  kapstd_uv(idepth) = gr(idepth)*(taufstd(idepth+1)-taufstd(idepth))/(ppf(idepth+1)-ppf(idepth))
END DO

!Calculate the absorption opacities
DO ise = 1, nse
  DO ilambda = 1, nuv
    !if current wavelength matches wavelength in cross section array
      kapnu_uv(:,ilambda) = kapnu_uv(:,ilambda) + &
        se(ise,ilambda)/unitL**2*amol(:,semol(ise))/(mu(:)*amu)
      kapfnu_uv(:,ilambda) = kapfnu_uv(:,ilambda) +  &
        se(ise,ilambda)/unitL**2*amolf(:,semol(ise))/(muf(:)*amu)
  END DO
END DO

!Calculate the scattering opacities
DO ilambda = 1, nuv
  DO iray = 1, num_ray_scatt
    alphap(iray) = (nr(iray)**2-1.)/(4.*pi*nl)
    sig(iray)  = 8.*pi/3.*alphap(iray)**2*(6.+3.*p_fac(iray))/(6.-7.*p_fac(iray))

    kscatnu_uv(:,ilambda) = (sig(iray)*(2.*pi/nu_loc(ilambda))**4)* &
      amol(:,imol(rayname(iray)))/(mu(:)*amu)
    kscatfnu_uv(:,ilambda) = (sig(iray)*(2.*pi/nu_loc(ilambda))**4)* &
      amolf(:,imol(rayname(iray)))/(muf(:)*amu)
  END DO
END DO

!Calculate epsnu and total opacities
WHERE (kapnu_uv + kscatnu_uv /=0.)
  epsnu_uv = kapnu_uv / (kapnu_uv + kscatnu_uv)
  kapnu_uv = kapnu_uv + kscatnu_uv
END WHERE

!Calculate epsfnu and total opacities
WHERE (kapfnu_uv + kscatfnu_uv /=0.)
  epsfnu_uv = kapfnu_uv / (kapfnu_uv + kscatfnu_uv)
  kapfnu_uv = kapfnu_uv + kscatfnu_uv
END WHERE

!Calculate the direct uv flux
CALL startrans(nuv,kapnu_uv,kapstd_uv,fstar_uv,Irad_uv_loc)

!Calculate the source function
IF (scatter_uv) THEN
  snu_uv = (1.-epsfnu_uv)*abs(fstar_uv/murad)/(4.*pi)
ELSE
!compute planck function
  snu_uv = 0.
END IF

! boundary condition at the bottom layer
bnu_uv = snu_uv
IF (nit_neq==1) THEN
  snu_scat_uv = snu_uv
END IF

!call radtrans routine to calculate flux including scattering
CALL radtrans(nrays_uv,scatter_uv,nuv,kapnu_uv,kapstd_uv,epsfnu_uv, &
  snu_uv,bnu_uv,hnup_uv,hnud_uv,jnup_uv,jnud_uv,snu_scat_uv)

!Combine direct and indirect flux
DO idepth=1,ndepth+1
  fuv_nu(idepth,1:nuv) = -(fstar_uv(idepth,1:nuv)+4.*pi* &
    (hnup_uv(idepth,1:nuv)+hnud_uv(idepth,1:nuv)))
END DO

! Convert back into photon s-1 cm-2 nm-1
fuv_nu = fuv_nu/1.0E7*(unitL**3./unitT)

! fuv = 0.
! !Compute integrated UV flux
! DO idepth = 1,ndepth+1
!   DO ilambda = 1,nuv
!     fuv(idepth) = fuv(idepth)+(fstar_uv(idepth,ilambda)+4.*pi* &
!       (hnup_uv(idepth,ilambda)+hnud_uv(idepth,ilambda))) &
!       /(1./unitT/unitL**2)*hplanck*cvel/(ilambda*1E-7/unitL)*dlambda
!   END DO
! END DO

END SUBROUTINE get_uv_flux

!----------------------
!Calculate the rate limits
!----------------------
SUBROUTINE get_rate_limits(node_loc,N_loc)

USE mod_chem
USE mod_chem_neq
USE mod_grid
USE mod_cst

IMPLICIT NONE

INTEGER,INTENT(IN) :: node_loc
REAL,INTENT(IN)    :: N_loc(node_loc)

!Local variables
INTEGER :: ireac, ir, ip, idepth, ndepth_max, idepth_lim
REAL    :: prls_loc, dn
REAL    :: ntot_loc(ndepth+1)
REAL    :: inc_reac(ndepth+1,nreac)

ndepth_max = ndepth_neq
IF (ndepth==1) ndepth_max = 1

rate_lim = 1.

!Calculate depth from which to limit reactions
!Default = 0.9*ndepth -> ndepth_max

if (depth_lim == 0.) then
  idepth_lim = 1
else
  idepth_lim = int(ndepth*depth_lim)
end if

DO idepth = idepth_lim,ndepth_max
  ntot_loc(idepth) = sum(N_loc(inloc(idepth,1:nmol_neq)))
  DO ireac = 1,nreac
    !compute the production/loss term prls
    prls_loc = kf(idepth,ireac)

    DO ir = 1,nb_reac(ireac)
      prls_loc = prls_loc*N_loc(inloc(idepth,imol_reac(ireac,ir)))
    END DO

    dn = prls_loc

    !compute the production/loss term prls
    prls_loc = kr(idepth,ireac)

    DO ip = 1,nb_prod(ireac)
      prls_loc = prls_loc*N_loc(inloc(idepth,imol_prod(ireac,ip)))
    END DO

    dn = abs(dn - prls_loc)*1E3*dt_lim
    inc_reac(idepth,ireac) = 1

    DO ir = 1,nb_reac(ireac)
      IF (dn>N_loc(inloc(idepth,imol_reac(ireac,ir)))) THEN
        inc_reac(idepth,ireac) = 0
      END IF
    END DO

    DO ip = 1,nb_prod(ireac)
      IF (dn>N_loc(inloc(idepth,imol_prod(ireac,ip)))) THEN
        inc_reac(idepth,ireac) = 0
      END IF
    END DO

    IF (idepth>1) THEN
      IF (inc_reac(idepth-1,ireac) == 0) inc_reac(idepth,ireac) = 0
    END IF

    IF ((.not. photochem) .or. (photochem .and. (.not. modparam(ireac,1)==nfile))) THEN
      IF (inc_reac(idepth,ireac) == 0 ) THEN
        IF (prls_loc/(nnf(idepth)/unitL**3)*dt_lim>1.) THEN
          rate_lim(idepth,ireac) = 1./(prls_loc/(nnf(idepth)/unitL**3)*dt_lim)
        END IF
      END IF

      IF (prls_loc<nn_lim/tmax) rate_lim(idepth,ireac) = 0.
      IF (ttf(idepth) <tt_lim) rate_lim(idepth,ireac) = 0.
    END IF
  END DO
END DO

END SUBROUTINE get_rate_limits

!----------------------
!Routine to calculate dn/dt, used by DLSODES
!----------------------
SUBROUTINE dndt(node_loc,t,N_loc,dndt_loc)

USE mod_chem_neq
USE mod_grid

IMPLICIT NONE

INTEGER            :: node_loc, i
REAL               :: t
REAL               :: N_loc(node_loc),dndt_loc(node_loc)
INTEGER            :: idepth_loc, ireac, ip, ir,ndepth_max
REAL               :: prls_loc
REAL               :: phi(ndepth,nmol_neq)

dndt_loc(1:node_loc) = 0.
ndepth_max = ndepth_neq

IF (ndepth==1) ndepth_max = 1

DO idepth_loc = 1,ndepth_max
  DO ireac = 1,nreac
    IF (.not. modParam(ireac,2) == 0) THEN
      ! compute the production/loss term prls
      IF (rate_limiter) THEN
        prls_loc = kf(idepth_loc,ireac)*rate_lim(idepth_loc,ireac)
      ELSE
        prls_loc = kf(idepth_loc,ireac)
      END IF

      DO ir = 1,nb_reac(ireac)
        prls_loc = prls_loc*N_loc(inloc(idepth_loc,imol_reac(ireac,ir)))
      END DO

      ! complete the derivatives with the production terms
      DO ip = 1,nb_prod(ireac)
        dndt_loc(inloc(idepth_loc,imol_prod(ireac,ip))) = dndt_loc(inloc(idepth_loc,imol_prod(ireac,ip))) + prls_loc
      END DO

      ! complete the derivatives with the loss terms
      DO ir = 1,nb_reac(ireac)
        dndt_loc(inloc(idepth_loc,imol_reac(ireac,ir))) = dndt_loc(inloc(idepth_loc,imol_reac(ireac,ir))) - prls_loc
      END DO
    END IF

    ! if reversed compute the production/loss term of the reversed reactions
    IF (modParam(ireac,2) == -1) THEN
      !compute the production/loss term prls
      IF (rate_limiter) THEN
        prls_loc = kr(idepth_loc,ireac)*rate_lim(idepth_loc,ireac)
      ELSE
        prls_loc = kr(idepth_loc,ireac)
      END IF

      DO ip = 1,nb_prod(ireac)
        prls_loc = prls_loc*N_loc(inloc(idepth_loc,imol_prod(ireac,ip)))
      END DO

      ! complete the derivatives with the production terms
      DO ir = 1,nb_reac(ireac)
        dndt_loc(inloc(idepth_loc,imol_reac(ireac,ir))) = dndt_loc(inloc(idepth_loc,imol_reac(ireac,ir))) + prls_loc
      END DO

      ! complete the derivatives with the loss terms
      DO ip = 1,nb_prod(ireac)
        dndt_loc(inloc(idepth_loc,imol_prod(ireac,ip))) = dndt_loc(inloc(idepth_loc,imol_prod(ireac,ip))) - prls_loc
      END DO
    END IF

  END DO

END DO

!if vertical mixing, calculate vertical flux and add to prod/loss term
IF (mixing) THEN

  !calculate vertical flux
  DO i = 1, nmol_neq
    phi(:,i) = get_phi(i,i,node_loc,nmol_neq,inloc,N_loc)
  END DO

  !calculate gradient of the flux
  DO idepth_loc = 2,ndepth_max-1
    IF (Rp /= 0.) THEN
      dndt_loc(inloc(idepth_loc,1:nmol_neq)) = dndt_loc(inloc(idepth_loc,1:nmol_neq)) - (phi(idepth_loc,1:nmol_neq) &
        *rr(idepth_loc)**2-phi(idepth_loc-1,1:nmol_neq)*rr(idepth_loc-1)**2)/drf(idepth_loc)/rrf(idepth_loc)**2
    ELSE
      dndt_loc(inloc(idepth_loc,1:nmol_neq)) = dndt_loc(inloc(idepth_loc,1:nmol_neq)) - &
        (phi(idepth_loc,1:nmol_neq)- phi(idepth_loc-1,1:nmol_neq))/drf(idepth_loc)
    END IF
  END DO

  !Boundary conditions
  !Take imaginary phi = 0 for idepth < 1 and idepth > ndepth_max
  !i.e. zero flux out of top and bottom
  IF (Rp /= 0.) THEN
    dndt_loc(inloc(1,:)) = dndt_loc(inloc(1,:)) - ((phi(1,:)*rr(1)**2)-(0.*rrf(1)**2))/drf(1)/rrf(1)**2
    dndt_loc(inloc(ndepth_max,:)) = dndt_loc(inloc(ndepth_max,:)) - (0.*rrf(ndepth_max)**2-phi(ndepth_max-1,:) &
      *rr(ndepth_max-1)**2)/drf(ndepth_max)/rrf(ndepth_max)**2
  ELSE
    dndt_loc(inloc(1,:)) = dndt_loc(inloc(1,:)) - (phi(1,:)-0.)/drf(1)
    dndt_loc(inloc(ndepth_max,:)) = dndt_loc(inloc(ndepth_max,:)) - (0.-phi(ndepth_max-1,:))/drf(ndepth_max)
  END IF

END IF

END SUBROUTINE dndt

!----------------------
!Routine to calculate Jacobian matrix
!By default, DLSODES generates the Jacobian internally
!but, alternatively, we can provide a user-generated
!Jacobian with this routine
!Controlled by IMF
!----------------------
SUBROUTINE jac(node_loc,t,N_loc,j,ian,jan,jacj)

USE mod_chem
USE mod_chem_neq
USE mod_grid
USE mod_cst

IMPLICIT NONE

INTEGER :: node_loc,j,jdepth,jmol, ndepth_max,ip,ir,ireac
REAL :: t,prls_loc,dzj,dzjp1,dzjm1,dzjp2,dzjm2
REAL,DIMENSION(node_loc) ::  N_loc,jacj,ian,jan
REAL,DIMENSION(ndepth+1) :: ntot_loc
INTEGER :: jinreac

jacj(1:node_loc) = 0.

ndepth_max = ndepth_neq
IF (ndepth==1) ndepth_max = 1

jdepth = (j-1)/nmol_neq+1
jmol   = j-(jdepth-1)*nmol_neq

IF (mixing) THEN
  IF (jdepth==1) THEN
    dzj   = (rrf(2)-rrf(1))*unitL
    dzjp1 = (rrf(3)-rrf(1))*unitL
    dzjp2 = (rrf(4)-rrf(2))*unitL
  ELSE IF (jdepth==2) THEN
    dzjm1 = (rrf(2)-rrf(1))*unitL
    dzj   = (rrf(3)-rrf(1))*unitL
    dzjp1 = (rrf(4)-rrf(2))*unitL
    dzjp2 = (rrf(5)-rrf(3))*unitL
  ELSE IF (jdepth==3) THEN
    dzjm2 = (rrf(2)-rrf(1))*unitL
    dzjm1 = (rrf(3)-rrf(1))*unitL
    dzj   = (rrf(4)-rrf(2))*unitL
    dzjp1 = (rrf(5)-rrf(3))*unitL
    dzjp2 = (rrf(6)-rrf(4))*unitL
  ELSE IF (jdepth==ndepth+1) THEN
    dzj   = (rrf(ndepth+1)-rrf(ndepth  ))*unitL
    dzjm1 = (rrf(ndepth+1)-rrf(ndepth-1))*unitL
    dzjm2 = (rrf(ndepth  )-rrf(ndepth-2))*unitL
  ELSE IF (jdepth==ndepth) THEN
    dzjp1 = (rrf(ndepth+1)-rrf(ndepth  ))*unitL
    dzj   = (rrf(ndepth+1)-rrf(ndepth-1))*unitL
    dzjm1 = (rrf(ndepth  )-rrf(ndepth-2))*unitL
    dzjm2 = (rrf(ndepth-1)-rrf(ndepth-3))*unitL
  ELSE IF (jdepth==ndepth-1) THEN
    dzjp2 = (rrf(ndepth+1)-rrf(ndepth  ))*unitL
    dzjp1 = (rrf(ndepth+1)-rrf(ndepth-1))*unitL
    dzj   = (rrf(ndepth  )-rrf(ndepth-2))*unitL
    dzjm1 = (rrf(ndepth-1)-rrf(ndepth-3))*unitL
    dzjm2 = (rrf(ndepth-2)-rrf(ndepth-4))*unitL
  ELSE
    dzjp2 = (rrf(jdepth+3)-rrf(jdepth+1))*unitL
    dzjp1 = (rrf(jdepth+2)-rrf(jdepth  ))*unitL
    dzj   = (rrf(jdepth+1)-rrf(jdepth-1))*unitL
    dzjm1 = (rrf(jdepth  )-rrf(jdepth-2))*unitL
    dzjm2 = (rrf(jdepth-1)-rrf(jdepth-3))*unitL
  END IF
END IF

DO ireac = 1,nreac
  IF(.not. modParam(ireac,2) == 0) THEN

    ! compute the production/loss term prls
    prls_loc = kf(jdepth,ireac)

    jinreac = 0
    DO ir = 1,nb_reac(ireac)
      IF (imol_reac(ireac,ir)==jmol) THEN
        jinreac= jinreac+1
        IF (jinreac>1) prls_loc = prls_loc*N_loc(inloc(jdepth,imol_reac(ireac,ir)))
      ELSE
        prls_loc = prls_loc*N_loc(inloc(jdepth,imol_reac(ireac,ir)))
      END IF
    END DO

    prls_loc = prls_loc*jinreac

    IF (jinreac>0) THEN

      ! complete the derivatives with the production terms
      DO ip = 1,nb_prod(ireac)
        jacj(inloc(jdepth,imol_prod(ireac,ip))) = jacj(inloc(jdepth,imol_prod(ireac,ip))) + prls_loc
      END DO

      ! complete the derivatives with the loss terms
      DO ir = 1,nb_reac(ireac)
        jacj(inloc(jdepth,imol_reac(ireac,ir))) = jacj(inloc(jdepth,imol_reac(ireac,ir))) - prls_loc
      END DO

    END IF
  END IF

  ! if reversed compute the production/loss term of the reversed reactions
  IF (modParam(ireac,2) == -1) THEN

    !compute the production/loss term prls
    prls_loc = kr(jdepth,ireac)

    jinreac = 0
    DO ip = 1,nb_prod(ireac)
      IF (imol_prod(ireac,ip)==jmol) THEN
        jinreac=jinreac+1
      IF (jinreac>1) prls_loc = prls_loc*N_loc(inloc(jdepth,imol_prod(ireac,ip)))
      ELSE
        prls_loc = prls_loc*N_loc(inloc(jdepth,imol_prod(ireac,ip)))
      END IF
    END DO

    prls_loc = prls_loc*jinreac

    IF (jinreac>0) THEN
      ! complete the derivatives with the production terms
      DO ir = 1,nb_reac(ireac)
        jacj(inloc(jdepth,imol_reac(ireac,ir))) = jacj(inloc(jdepth,imol_reac(ireac,ir))) + prls_loc
      END DO

      ! complete the derivatives with the loss terms
      DO ip = 1,nb_prod(ireac)
        jacj(inloc(jdepth,imol_prod(ireac,ip))) = jacj(inloc(jdepth,imol_prod(ireac,ip))) - prls_loc

      END DO
    END IF

  END IF
END DO

IF (mixing) THEN
  IF (Rp /=0.) THEN
    IF (jdepth<=2) THEN

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth+1) = sum(n_loc(inloc(jdepth+1,:))*molmass(:))
      ntot_loc(jdepth+2) = sum(n_loc(inloc(jdepth+2,:))*molmass(:))

      jacj(inloc(jdepth,jmol)) = jacj(inloc(jdepth,jmol)) - &
        1/dzj*(-kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1)) &
        *rrf(jdepth+1)**2/rrf(jdepth)**2

      jacj(inloc(jdepth+2,jmol)) = jacj(inloc(jdepth+2,jmol)) - &
        1/dzjp2*(kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1)) &
        *rrf(jdepth+1)**2/rrf(jdepth)**2

    ELSE IF (jdepth>=ndepth_max-1) THEN

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth-1) = sum(n_loc(inloc(jdepth-1,:))*molmass(:))
      ntot_loc(jdepth-2) = sum(n_loc(inloc(jdepth-2,:))*molmass(:))

      jacj(inloc(jdepth,jmol)) = jacj(inloc(jdepth,jmol)) - &
        1/dzj*(-kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1)) &
        *rrf(jdepth-1)**2/rrf(jdepth)**2

      jacj(inloc(jdepth-2,jmol)) = jacj(inloc(jdepth-2,jmol)) - &
        1/dzjm2*(kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1)) &
        *rrf(jdepth-1)**2/rrf(jdepth)**2

    ELSE

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth-1) = sum(n_loc(inloc(jdepth-1,:))*molmass(:))
      ntot_loc(jdepth-2) = sum(n_loc(inloc(jdepth-2,:))*molmass(:))

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth+1) = sum(n_loc(inloc(jdepth+1,:))*molmass(:))
      ntot_loc(jdepth+2) = sum(n_loc(inloc(jdepth+2,:))*molmass(:))

      jacj(inloc(jdepth,jmol)) = jacj(inloc(jdepth,jmol)) - &
        1/dzj*(-kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1) &
        *rrf(jdepth+1)**2-kzz(jdepth-1)*ntot_loc(jdepth-1)/&
        (ntot_loc(jdepth)*dzjm1)*rrf(jdepth-1)**2)/rrf(jdepth)**2

      jacj(inloc(jdepth-2,jmol)) = jacj(inloc(jdepth-2,jmol)) - &
        1/dzjm2*(kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1)) &
        *rrf(jdepth-1)**2/rrf(jdepth-2)**2

      jacj(inloc(jdepth+2,jmol)) = jacj(inloc(jdepth+2,jmol)) - &
        1/dzjp2*(kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1)) &
        *rrf(jdepth+1)**2/rrf(jdepth+2)**2

    END IF
  ELSE
    IF (jdepth<=2) THEN

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth+1) = sum(n_loc(inloc(jdepth+1,:))*molmass(:))
      ntot_loc(jdepth+2) = sum(n_loc(inloc(jdepth+2,:))*molmass(:))

      jacj(inloc(jdepth,jmol)) = jacj(inloc(jdepth,jmol)) - &
        1/dzj*(-kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1))

      jacj(inloc(jdepth+2,jmol)) = jacj(inloc(jdepth+2,jmol)) - &
        1/dzjp2*(kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1))

    ELSE IF (jdepth>=ndepth_max-1) THEN

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth-1) = sum(n_loc(inloc(jdepth-1,:))*molmass(:))
      ntot_loc(jdepth-2) = sum(n_loc(inloc(jdepth-2,:))*molmass(:))

      jacj(inloc(jdepth,jmol)) = jacj(inloc(jdepth,jmol)) - &
        1/dzj*(-kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1))

      jacj(inloc(jdepth-2,jmol)) = jacj(inloc(jdepth-2,jmol)) - &
        1/dzjm2*(kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1))

    ELSE

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth-1) = sum(n_loc(inloc(jdepth-1,:))*molmass(:))
      ntot_loc(jdepth-2) = sum(n_loc(inloc(jdepth-2,:))*molmass(:))

      ntot_loc(jdepth) = sum(n_loc(inloc(jdepth,:))*molmass(:))
      ntot_loc(jdepth+1) = sum(n_loc(inloc(jdepth+1,:))*molmass(:))
      ntot_loc(jdepth+2) = sum(n_loc(inloc(jdepth+2,:))*molmass(:))

      jacj(inloc(jdepth,jmol)) = jacj(inloc(jdepth,jmol)) - &
        1/dzj*(-kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1) &
        -kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1))

      jacj(inloc(jdepth-2,jmol)) = jacj(inloc(jdepth-2,jmol)) - &
        1/dzjm2*(kzz(jdepth-1)*ntot_loc(jdepth-1)/(ntot_loc(jdepth)*dzjm1))

      jacj(inloc(jdepth+2,jmol)) = jacj(inloc(jdepth+2,jmol)) - &
        1/dzjp2*(kzz(jdepth+1)*ntot_loc(jdepth+1)/(ntot_loc(jdepth)*dzjp1))

    END IF
  END IF
END IF

END SUBROUTINE jac

!----------------------
!SUBROUTINE to write the PT profile
!----------------------
SUBROUTINE write_pt_profile(fwrite)

USE mod_param
USE mod_util
USE mod_cst
USE mod_eos
USE mod_conv
USE mod_radtrans
USE mod_chem_neq
USE mod_opacity
USE mod_chem
USE mod_solver
USE mod_grid

IMPLICIT NONE

INCLUDE 'netcdf.inc'

INTEGER :: id_file,id_dim,id_dimc,nvar,id_lname,lname_loc
INTEGER,ALLOCATABLE:: id_var(:)

CHARACTER(*)  :: fwrite

lname_loc = 60
nvar  = 63

ALLOCATE(id_var(nvar))
id_var = 0

IF (mype == cputerm) THEN
  CALL nf(nf_create(fwrite,nf_clobber,id_file))

  CALL nf(nf_def_dim(id_file,'lname',lname_loc,id_lname))
  CALL nf(nf_def_dim(id_file,'nlevel',ndepth+1,id_dim))
  CALL nf(nf_def_var(id_file,'temperature',nf_double,1,id_dim,id_var(1)))
  CALL nf(nf_def_var(id_file,'pressure',nf_double,1,id_dim,id_var(2)))
  CALL nf(nf_def_var(id_file,'gradad',nf_double,1,id_dim,id_var(3)))
  CALL nf(nf_def_var(id_file,'pressure_si',nf_double,1,id_dim,id_var(4)))

  IF (ndim==2) THEN
    CALL nf(nf_def_dim(id_file,'ncol',ntheta+1,id_dimc))

    CALL nf(nf_def_var(id_file,'temperature_2d',nf_double,2,(/id_dim,id_dimc/),id_var(5)))
    CALL nf(nf_def_var(id_file,'pressure_2d',nf_double,2,(/id_dim,id_dimc/),id_var(6)))
    CALL nf(nf_def_var(id_file,'wind_2d',nf_double,2,(/id_dim,id_dimc/),id_var(7)))
    CALL nf(nf_def_var(id_file,'dwind_2d',nf_double,2,(/id_dim,id_dimc/),id_var(8)))
    CALL nf(nf_def_var(id_file,'windz_2d',nf_double,2,(/id_dim,id_dimc/),id_var(9)))

    CALL nf(nf_def_var(id_file,'Rmax',nf_double,0,id_dim,id_var(10)))
    CALL nf(nf_def_var(id_file,'u_cst',nf_double,0,id_dim,id_var(11)))
    CALL nf(nf_def_var(id_file,'alpha_dayside_cst',nf_double,0,id_dim,id_var(12)))
    CALL nf(nf_def_var(id_file,'alpha_nighside_cst',nf_double,0,id_dim,id_var(13)))
  END IF

  CALL nf(nf_def_var(id_file,'frad',nf_double,1,id_dim,id_var(14)))
  CALL nf(nf_def_var(id_file,'fconv',nf_double,1,id_dim,id_var(15)))
  CALL nf(nf_def_var(id_file,'kconv',nf_double,1,id_dim,id_var(16)))
  CALL nf(nf_def_var(id_file,'alpha',nf_double,0,id_dim,id_var(17)))
  CALL nf(nf_def_var(id_file,'pstop_conv',nf_double,0,id_dim,id_var(18)))
  CALL nf(nf_def_var(id_file,'eps_gradad',nf_double,0,id_dim,id_var(19)))

  CALL nf(nf_def_var(id_file,'gamma',nf_double,0,id_dim,id_var(20)))
  CALL nf(nf_def_var(id_file,'pmin_conv',nf_double,0,id_dim,id_var(21)))
  CALL nf(nf_def_var(id_file,'pmax_conv',nf_double,0,id_dim,id_var(22)))
  CALL nf(nf_def_var(id_file,'gradad_smooth',nf_int,0,id_dim,id_var(23)))
  CALL nf(nf_def_var(id_file,'kereos_smooth',nf_int,0,id_dim,id_var(24)))

  CALL nf(nf_def_var(id_file,'Rp',nf_double,0,id_dim,id_var(25)))
  CALL nf(nf_def_var(id_file,'pp_Rp',nf_double,0,id_dim,id_var(26)))

  CALL nf(nf_def_var(id_file,'logg',nf_double,0,id_dim,id_var(27)))
  CALL nf(nf_def_var(id_file,'teff',nf_double,0,id_dim,id_var(28)))
  CALL nf(nf_def_var(id_file,'MdH',nf_double,0,id_dim,id_var(29)))

  CALL nf(nf_def_var(id_file,'psurf',nf_double,0,id_dim,id_var(30)))
  CALL nf(nf_def_var(id_file,'taumin',nf_double,0,id_dim,id_var(31)))
  CALL nf(nf_def_var(id_file,'taumax',nf_double,0,id_dim,id_var(32)))

  CALL nf(nf_def_var(id_file,'nband_std',nf_int,0,id_dim,id_var(33)))
  CALL nf(nf_def_var(id_file,'nband',nf_int,0,id_dim,id_var(34)))
  CALL nf(nf_def_var(id_file,'nfreq',nf_int,0,id_dim,id_var(35)))
  CALL nf(nf_def_var(id_file,'nkmix',nf_int,0,id_dim,id_var(36)))

  CALL nf(nf_def_var(id_file,'mixing',nf_int,0,id_dim,id_var(37)))
  CALL nf(nf_def_var(id_file,'photochem',nf_int,0,id_dim,id_var(38)))
  CALL nf(nf_def_var(id_file,'rate_limiter',nf_int,0,id_dim,id_var(39)))
  CALL nf(nf_def_var(id_file,'dt_lim',nf_double,0,id_dim,id_var(40)))
  CALL nf(nf_def_var(id_file,'atol',nf_double,0,id_dim,id_var(41)))
  CALL nf(nf_def_var(id_file,'Nmin',nf_double,0,id_dim,id_var(42)))
  CALL nf(nf_def_var(id_file,'kzzcst',nf_double,0,id_dim,id_var(43)))

  CALL nf(nf_def_var(id_file,'scatter',nf_int,0,id_dim,id_var(44)))
  CALL nf(nf_def_var(id_file,'irrad',nf_int,0,id_dim,id_var(45)))
  CALL nf(nf_def_var(id_file,'murad',nf_double,0,id_dim,id_var(46)))
  CALL nf(nf_def_var(id_file,'fred',nf_double,0,id_dim,id_var(47)))

  CALL nf(nf_def_var(id_file,'nkap',nf_int,0,id_dim,id_var(48)))
  CALL nf(nf_def_var(id_file,'kap_smooth',nf_int,0,id_dim,id_var(49)))
  CALL nf(nf_def_var(id_file,'ppmin_smooth',nf_double,0,id_dim,id_var(50)))
  CALL nf(nf_def_var(id_file,'ppmax_smooth',nf_double,0,id_dim,id_var(51)))
  CALL nf(nf_def_var(id_file,'kerkap_smooth',nf_int,0,id_dim,id_var(52)))

  CALL nf(nf_def_var(id_file,'mucst',nf_double,0,id_dim,id_var(53)))
  CALL nf(nf_def_var(id_file,'tfreeze_eq',nf_double,0,id_dim,id_var(54)))
  CALL nf(nf_def_var(id_file,'fcoeff',nf_char,1,id_lname,id_var(55)))
  CALL nf(nf_def_var(id_file,'taufstd',nf_double,1,id_dim,id_var(56)))
  CALL nf(nf_def_var(id_file,'rrf',nf_double,1,id_dim,id_var(57)))
  CALL nf(nf_def_var(id_file,'rhof',nf_double,1,id_dim,id_var(58)))
  CALL nf(nf_def_var(id_file,'grf',nf_double,1,id_dim,id_var(59)))
  CALL nf(nf_def_var(id_file,'nnf',nf_double,1,id_dim,id_var(60)))
  CALL nf(nf_def_var(id_file,'err_hydro',nf_double,1,id_dim,id_var(61)))
  CALL nf(nf_def_var(id_file,'err_energy',nf_double,1,id_dim,id_var(62)))
  if (calc_radtime) then
     CALL nf(nf_def_var(id_file,'radtime',nf_double,1,id_dim,id_var(63)))
  endif
  CALL nf(nf_enddef(id_file))

  CALL nf(nf_put_vara_double(id_file,id_var(1),1,ndepth+1,ttf*unitK))
  CALL nf(nf_put_vara_double(id_file,id_var(2),1,ndepth+1,ppf*unitP))
  CALL nf(nf_put_vara_double(id_file,id_var(3),1,ndepth+1,gradad))
  CALL nf(nf_put_vara_double(id_file,id_var(4),1,ndepth+1,ppf*unitP*0.1)) !convert to SI: Pa

  IF (ndim==2) THEN
    CALL nf(nf_put_vara_double(id_file,id_var(5),(/1,1/),(/ndepth+1,ntheta+1/),ttf_2d(1:ndepth+1,1:ntheta+1)*unitK))
    CALL nf(nf_put_vara_double(id_file,id_var(6),(/1,1/),(/ndepth+1,ntheta+1/),ppf_2d(1:ndepth+1,1:ntheta+1)*unitP))
    CALL nf(nf_put_vara_double(id_file,id_var(7),(/1,1/),(/ndepth+1,ntheta+1/),utf_2d(1:ndepth+1,1:ntheta+1)*unitL/unitT))
    CALL nf(nf_put_vara_double(id_file,id_var(8),(/1,1/),(/ndepth+1,ntheta+1/),dupf_2d(1:ndepth+1,1:ntheta+1)*unitL/unitT))
    CALL nf(nf_put_vara_double(id_file,id_var(9),(/1,1/),(/ndepth+1,ntheta+1/),urf_2d(1:ndepth+1,1:ntheta+1)*unitL/unitT))

    CALL nf(nf_put_var_double(id_file,id_var(10),Rmax*unitL))
    CALL nf(nf_put_var_double(id_file,id_var(11),u_cst*unitL/unitT))
    CALL nf(nf_put_var_double(id_file,id_var(12),alpha_dayside_cst))
    CALL nf(nf_put_var_double(id_file,id_var(13),alpha_nightside_cst))

  END IF

  CALL nf(nf_put_vara_double(id_file,id_var(14),1,ndepth+1,frad*unitIF))
  CALL nf(nf_put_vara_double(id_file,id_var(15),1,ndepth+1,fconv*unitIF))
  CALL nf(nf_put_vara_double(id_file,id_var(16),1,ndepth+1,kconv*unitL**2/unitT))

  CALL nf(nf_put_var_double(id_file,id_var(17),alpha))
  CALL nf(nf_put_var_int   (id_file,id_var(18),pstop_conv))
  CALL nf(nf_put_var_double(id_file,id_var(19),eps_gradad))

  CALL nf(nf_put_var_double(id_file,id_var(20),gamma))
  CALL nf(nf_put_var_double(id_file,id_var(21),pmin_conv*unitP))
  CALL nf(nf_put_var_double(id_file,id_var(22),pmin_conv*unitP))
  CALL nf(nf_put_var_int   (id_file,id_var(23),gradad_smooth))
  CALL nf(nf_put_var_double(id_file,id_var(24),kereos_smooth))

  CALL nf(nf_put_var_double(id_file,id_var(25),Rp*unitL))
  CALL nf(nf_put_var_double(id_file,id_var(26),pp_Rp*unitP))


  CALL nf(nf_put_var_double(id_file,id_var(27),logg))
  CALL nf(nf_put_var_double(id_file,id_var(28),teff))
  CALL nf(nf_put_var_double(id_file,id_var(29),MdH))

  CALL nf(nf_put_var_double(id_file,id_var(30),psurf*unitP))
  CALL nf(nf_put_var_double(id_file,id_var(31),taumin))
  CALL nf(nf_put_var_double(id_file,id_var(32),taumax))

  CALL nf(nf_put_var_int   (id_file,id_var(33),nband_std))
  CALL nf(nf_put_var_int   (id_file,id_var(34),nband))
  CALL nf(nf_put_var_int   (id_file,id_var(35),nfreq))
  CALL nf(nf_put_var_int   (id_file,id_var(36),nkmix))

  CALL nf(nf_put_var_int   (id_file,id_var(37),mixing))
  CALL nf(nf_put_var_int   (id_file,id_var(38),photochem))
  CALL nf(nf_put_var_int   (id_file,id_var(39),rate_limiter))
  CALL nf(nf_put_var_double(id_file,id_var(40),dt_lim))
  CALL nf(nf_put_var_double(id_file,id_var(41),lsode_atol))
  CALL nf(nf_put_var_double(id_file,id_var(42),nmin))
  CALL nf(nf_put_var_double(id_file,id_var(43),kzzcst*unitL**2/unitT))

  CALL nf(nf_put_var_int   (id_file,id_var(44),scatter))
  CALL nf(nf_put_var_int   (id_file,id_var(45),irrad))
  CALL nf(nf_put_var_double(id_file,id_var(46),murad))
  CALL nf(nf_put_var_double(id_file,id_var(47),fred))

  CALL nf(nf_put_var_int   (id_file,id_var(48),nkap))
  CALL nf(nf_put_var_int   (id_file,id_var(49),kap_smooth))
  CALL nf(nf_put_var_double(id_file,id_var(50),ppmin_smooth*unitP))
  CALL nf(nf_put_var_double(id_file,id_var(51),ppmax_smooth*unitP))
  CALL nf(nf_put_var_int   (id_file,id_var(52),kerkap_smooth))

  CALL nf(nf_put_var_double(id_file,id_var(53),mucst))
  CALL nf(nf_put_var_double(id_file,id_var(54),tfreeze_eq))
  CALL nf(nf_put_vara_text(id_file,id_var(55),1,lname_loc,fcoeff))
  CALL nf(nf_put_vara_double(id_file, id_var(56),1,ndepth+1,taufstd))
  CALL nf(nf_put_vara_double(id_file, id_var(57),1,ndepth+1,rrf*unitL))
  CALL nf(nf_put_vara_double(id_file, id_var(58),1,ndepth+1,rhof*unitD))
  CALL nf(nf_put_vara_double(id_file, id_var(59),1,ndepth+1,grf*unitL/unitT**2))
  CALL nf(nf_put_vara_double(id_file, id_var(60),1,ndepth+1,nnf/unitL**3))
  CALL nf(nf_put_vara_double(id_file, id_var(61),1,ndepth+1,err_hydro_prf))
  CALL nf(nf_put_vara_double(id_file, id_var(62),1,ndepth+1,err_energy_prf))

  if (calc_radtime) then
     CALL nf(nf_put_vara_double(id_file,id_var(63),1,ndepth+1,radtime*unitT))
  endif

  CALL nf(nf_close(id_file))

END IF

END SUBROUTINE write_pt_profile
