!=======================================================================
!     PROGRAM JKTEBOP           John Southworth  (jkt~astro.keele.ac.uk)
!                               Astrophysics Group, Keele University, UK
!-----------------------------------------------------------------------
! V(1) = surface brightness ratio      V(15) = third light
! V(2) = sum of fractional radii       V(16) = phase correction
! V(3) = ratio of stellar radii        V(17) = light scaling factor
! V(4) = linear LD for star A          V(18) = integration ring size (o)
! V(5) = linear LD for star B          V(19) = orbital period (days)
! V(6) = orbital inclination           V(20) = ephemeris timebase (days)
! V(7) = e cos(omega) OR ecentricity   V(21) = nonlinear LD for star A
! V(8) = e sin(omega) OR omega         V(22) = nonlinear LD for star B
! V(9) = gravity darkening 1           VEXTRA(1) = primary star radius
! V(10) = gravity darkening 2          VEXTRA(2) = secondary star radius
! V(11) = primary reflected light      VEXTRA(3) = stellar light ratio
! V(12) = secondary reflected light    VEXTRA(4) = eccentricity
! V(13) = stellar mass ratio           VEXTRA(5) = periastron longitude
! V(14) = tidal lead/lag angle (deg)   VEXTRA(6) = reduced chi-squared
! V(23-37) five lots of sine curve [T0,P,amplitude]
! V(38-67) five lots of polynomial [pivot,x,x2,x3,x4,x5]
!-----------------------------------------------------------------------
! Version 1:  Simplex minimisation algorithm and new input / output used
! Version 2:  Monte Carlo simulation and parametr perturbation algorithm
! Version 3:  Adjustment to Monte Carlo LD coeffs and input/output files
! Version 4:  Solves for sum of radii and convergence criterion modified
! Version 5:  Added TASK0 to find LD and GD coeffs.  Minor modifications
! Version 6:  Reflection and  scale factor  can all be fixed or adjusted
! Version 7:  Can use (e,w) or (ecosw,esinw).    SFACT modified to be in
!             magnitudes; observ'nal errors found; spherical star option
! Version 8:  Bootstrapping error analysis added and the output modified
! Version 9:  Command-line arguments allowed, Monte Carlo without param-
!             eter kicking option and fitting for period and Tzero added
! Version 10: Can now use 99999 datapoints. Whole code now in magnitudes
! Version 11: Bug fixes, tasks renumbered,  added sigma clipping and the
!             global fit procedures, but not thoroughly tested these yet
! Version 12: Nonlinear limb darkening law, fitting for times of minimum
!             light, and FMAX corrections included (from Alvaro Gimenez)
! Version 13: Removed  BILINEAR  and modified  TASK1  to just call JKTLD
!             Modified input file and arrays  for the diff types of data
! Version 14: Fixed the requirement for inputtd INTRING to be an integer
!             Fixed formal errors (when observational ones not supplied)
!             Sorted out numerical derivatives problem when  i ~ 90  deg
! Version 15: Added TASK 9 and cubic LD law, modified simulation output,
!             made MRQMIN a factor of 3 faster, and working on red noise
! Version 16: Added ability to include sine perturbations on parameters.
! Version 17: Added ability to specify a third light and its uncertainty
! Version 18: Made possible to fit for r1 and r2 instead of r1+r2 and k.
! Version 19: Added VARY=3 flag + finer fit phasing if r1 or r2 is small
! Version 20: Polynomial, optimisation check, TASK4 debug, DV centralise
! Version 21: Added input of ecosw and esinw observational  constraints.
! Version 22: Added input of  e and omega  as observational constraints.
! Version 23: (did not) fix a minor bug with calculation of third light.
! Version 24: Moved to gfortran compiler.  Adjusted date_and_time usage.
!             Fixed bug with third light. Added ECQUADPHASES subroutine.
! Version 25: Added numerical integration over long exposure times.
! Version 26: Fixed bug with TASK 5 and modified TASK 5 output slightly.
! Version 27: Fixed TASK 8 bug, converted all to real*8, 999999 datapnt.
! Version 28: Corrected ECQUADPHASES. All write(*) or print* => write(6)
! Last modified: 1st March 2012
!-----------------------------------------------------------------------
! Possible modifications in future:
! 1) Extend to WD2003 and WINK
! 2) Port to F90 or F95 to have long lines, modules, and improved output
! 3) Incorporate change of omega (apsidal motion)
! 4) Include radial velocities
! 5) Include a light-time effect
! 6) Allow for multiple light and radial velocity curves
! 7) Try LD power law proposed by Hestroffer (1997A+A...327..199H)
! 8) Add four-parameter LD law (Claret 2000A+A...363.1081C)
!-----------------------------------------------------------------------
! Miscellaneous notes:
! 1) Phase shift has been redefined compared to original EBOP so that it
!    corresponds directly to the phase of primary minimum.
! 2) MRQMIN adjusts coeffs only if VARY (called 'ia' in MRQMIN) is 1
! 3) If VARY=2 then the parameter is fixed during the initial fit but is
!    perturbed by a set amount (flat distribution) for later analyses.
! 4) If VARY(11) and/or  VARY(12) are "-1" then V(11) and/or V(12) are
!    calculated from the system geometry; if 0 they are fixed at the
!    input value and if 1 are freely adjusted to best fit.
! 5) If the mass ratio is <= 0 then both stars are assumed to be spheres
! 6) If ecosw > 5.0 then (ecosw,esinw) will be taken to be  (10+e,omega)
!    and fitting will occur using e and omega as parameters. e and omega
!    can be strongly correlated, but this option is useful if e is known
!    but omega isn't; this can happen for EBs exhibiting apsidal motion.
! 7) Observational errors are looked for in the  input light curve file.
!    If they are not found then  equal weight  is given  to each  point.
! 8) Nonlinear LD is now supported for the two-coefficient  logarithmic,
!    quadratic and square-root laws.  The type of law must be specified.
!    on input.   Star B can also be forced to the same coeffs as star A.
!    BUT: normalisation for logarithmic not possible (not got equations)
! 9) Fitting for times of minimum light is directly possible.  The cycle
!    numbers and times are inputted  on lines  immediately below all the
!     parameter lines in the input file.
! 10) EBOP results are symmetric about 90 degrees for inclination (which
!     means i=89 gives same answer as i=91),  which causes problems with
!     numerical derivativs when i>89.9. In this case some extra is added
!     to the numerical derivative to keep the solution slightly below 90
! 11) If input k is negative then  (r1+r2,k)  is interpreted as  (r1,r2)
! 12) If  VARY=3  then the parameter is optimised during all fits but is
!     not perturbed by a set amount in later analyses (eg. Monte Carlo).
!-----------------------------------------------------------------------
! Task numbers and purposes:
! (1) This outputs LD coefficients for a given Teff, logg, [M/H], Vmicro
! (2) This outputs a model light curve for fixed input parameters.
! (3) This fits a model to an observed light curve and outputs results.
! (4) This fits a model, rejects discrepant observations, and refits.
! (5) This does a pseudo-global minimisation by perturbing input params.
! (6) This investigates how different parameters vary around best fit.
! (7) This conducts bootstrapping simulations to find robust errors.
! (8) This conducts Monte Carlo simulations to find robust errors.
! (9) This conducts residual permutations to deal with correlated noise.
!-----------------------------------------------------------------------
! Language:  JKTEBOP is written in FORTRAN 77, using several extensions
!   to the ANSI standard:   ==   <=   <   >   >=   /=  
!   enddo  endif
!   Until version 24 it was only ever compiled in g77.
! g77 compiler: I only occasionally check if this works as the g77 comp-
!   iler is no longer supported.   To compile with g77 you should change
!   the way the DATE_AND_TIME intrinsic function is called.  To do this,
!   simply search for the lines containing "DTTIME*9",  uncomment these,
!   and comment out the lines containing "DTTIME*10".     I successfully
!   compiled JKTEBOP v25 on 29/10/2010 using  gcc version 3.4.6 20060404
!   (Red Hat 3.4.6-4) on a Scientific Linux PC. The compilation command:
!   g77 -O -Wuninitialized -fbounds-check -fno-automatic -o jktebop
! g95 compiler: this compiles successfully but I have not actually tried
!   to run the resulting executable file.    Compiler version used: "G95
!   (GCC 4.1.2 (g95 0.9
!) Jun 16 2010)"  running on a kubuntu 10.04 PC.
! gfortran compiler:  this compiles successfully and executes correctly.
!   The compiler version used last time I modified the current text was:
!   "GNU Fortran (Ubuntu 4.4.3-4ubuntu5) 4.4.3" running on kubuntu 10.04
! Intel-Fortran compiler: this is periodically checked and found to work
!   well. JKTEBOP versions v26 and earlier must be compiled with the -r8
!   command-line flag in order to avoid numerical noise arising from the
!   use of single-precision variables. JKTEBOP v27 onwards is all real*8
!=======================================================================
!=======================================================================
      SUBROUTINE TASK2 (V,LDTYPE,NPHASE,OPHASE,OMAG,OLP,OLS,
     &                                         OOREFL,OOLECL)
      implicit none                       
! Produces a model light curve
      real*8 V(67)                        
! IN: light  curve  parameters
      integer LDTYPE(2)                   
! IN: LD law type foreach star
      integer i,ERROR                    
! LOCAL: counters & error flag
      real*8 MAG,LP,LS,OREFL,OLECL                 
! LOCAL: EBOP/GETMODEL  output
      real*8 PHASE                       
! LOCAL:  phase for evaluation
      real*8 OMAG
      real*8 OPHASE                       
      real*8 OLP                    
      real*8 OLS                       
      real*8 OOREFL                       
      real*8 OOLECL                       
      real*8 GETMODEL                    
! FUNCTION: evaluate the model
      integer NSINE                      
! OUT: Numbrs of sines and L3s
      integer PSINE(5)                   
! OUT: Which par for each sine
      integer NPOLY,PPOLY(5)             
! OUT: Similar for polynomials

      real*8 HJD
      real*8 R1,R2
      integer NPHASE                     
! LOCAL: number of phases todo
      dimension OPHASE(*),OMAG(*),OLP(*),OLS(*),OOREFL(*),OOLECL(*)

      V(19) = 1.0d0          
! Set period to 1.0
      V(20) = 0.0d0          
! Set Tzero to 0.0
      LP = 0.0d0
      LS = 0.0d0
                                     
! NSINE=0 and NPOLY=0 and NUMINT=1
      if ( V(3) >= 0.0d0 ) then
        R1 = V(2) / (1.0d0 + V(3))
        R2 = V(2) / (1.0d0 + (1.0d0/V(3)))
      else
        R1 = V(2)
        R2 = abs(V(3))
      end if
      MAG=GETMODEL(V,LDTYPE,0,PSINE,0,PPOLY,V(20),1,LP,LS,1,0.0d0,
     &                                                OREFL,OLECL)

C      V(11) = 0.4d0 * (LS/(1.0d0-V(15))) * R1**2
C      V(12) = 0.4d0 * (LP/(1.0d0-V(15))) * R2**2

C      NPHASE = 10001
C      if ( R1 < 0.01d0 .or. R2 < 0.01d0 ) NPHASE = 100001
C
C      write(6,'(A40,A40)') ">> The reflection coefficients come from",
C     &                      " the system geometry, not the input file"
C
C      write(62,'(A47)')"#  PHASE  MAGNITUDE    L1         L2         L3"
C
      do i = 1,NPHASE
        HJD = OPHASE(i)
        MAG = GETMODEL(V,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,1,0.0d0,
     &                                                  OREFL,OLECL)
C        PHASE = (i-1) / dble(NPHASE-1)
C        HJD = V(20) + PHASE * V(19)
C        MAG = GETMODEL(V,LDTYPE,0,PSINE,0,PPOLY,HJD,1,LP,LS,1,0.0d0)
C        write (62,'(F8.6,4(1X,F10.6))') PHASE,MAG,LP,LS,V(15)
        OMAG(i) = MAG
        OLP(i)=LP
        OLS(i)=LS
        OOREFL(i)=OREFL
        OOLECL(i)=OLECL
      end do
C      close (62)

      END SUBROUTINE TASK2
!=======================================================================
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMODEL (V,LDTYPE,NSINE,PSINE,NPOLY,
     &          PPOLY,TIME,DTYPE,LA,LB,NUMINT,NINTERVAL,OREFL,OLECL)
           
! Output a predicted model value according to the parameters
           
! in array V. Precise meaning of the value depends on DTYPE.
           
! DTYPE=1  it outputs an EBOP magnitude for given time
           
! DTYPE=2  it outputs a light ratio for the given time
           
! DTYPE=3  outputs a time of eclipse for the given =CYCLE=
           
! DTYPE=4  it simply outputs the third light value
      implicit none
      real*8 V(67)                 
! IN: Photometric parameters
      integer LDTYPE(2)            
! IN: LD law type for the two stars
      real*8 TIME                  
! IN: The given TIME, PHASE or CYCLE
      integer DTYPE                
! IN: 1-6 depending on wanted result
      integer NSINE,PSINE(5)       
! IN: number and parameters of sines
      integer NPOLY,PPOLY(5)       
! IN: number and parameters of polys
      integer NUMINT               
! IN: Number of numerical integratns
      real*8 NINTERVAL             
! IN: Time interval for integrations
      real*8 LA,LB,OREFL,OLECL              
! OUT: Light produced by each star
      real*8 FMAG,LP,LS,GREFL,GLECL        
! LOCAL: LIGHT subroutine output
      real*8 ECC,OMEGA,ECOSW,ESINW 
! LOCAL: orbital shape parameters
      real*8 GETMIN,GETPHASE       
! FUNCTIONS
      real*8 FMAGSUM,LASUM,LBSUM
      integer i
      real*8 TIMEIN

      GETMODEL = 0.0d0

           
! DTYPE=1 for light curve datapoint
           
! DTYPE=2 for light ratio
           
! DTYPE=3 for times of minimum light
           
! DTYPE=4 for third light
           
! DTYPE=5 for e*cos(omega)
           
! DTYPE=6 for e*sin(omega)

      if ( DTYPE == 1 ) then
        if ( NUMINT == 1 ) then
         CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS,
     &                                                    GREFL,GLECL)
          OREFL = GREFL
          OLECL = GLECL
          LA = LP
          LB = LS
          GETMODEL = FMAG

        else if ( NUMINT > 1 ) then
          FMAGSUM = 0.0d0
          LASUM = 0.0d0
          LBSUM = 0.0d0
          CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS,
     &                                                   GREFL,GLECL)
          OREFL = GREFL
          OLECL = GLECL
          do i = 1,NUMINT
            TIMEIN  =  TIME  -  NINTERVAL / 86400.0d0 / dble(NUMINT) *
     &        ( dble(NUMINT) - 2.0d0*dble(i) + 1.0d0 ) / 2.0d0
          CALL LIGHT(V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIMEIN,FMAG,LP,LS,
     &                                                   GREFL,GLECL)
            OREFL = GREFL
            OLECL = GLECL
            FMAGSUM = FMAGSUM + FMAG
            LASUM = LASUM + LP
            LBSUM = LBSUM + LS
          end do
          FMAG = FMAGSUM / NUMINT
          LA = LASUM / NUMINT
          LB = LBSUM / NUMINT
          GETMODEL = FMAG
        else
          write(6,*)"NUMINT is less than 1 in function GETMODEL. Abort."
          write(6,*)"NUMINT =    ", NUMINT
          write(6,*)"NINTERVAL = ", NINTERVAL
          stop
        end if

      else if ( DTYPE == 2 ) then
        CALL LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,TIME,FMAG,LP,LS,
     &                                                 GREFL,GLECL)
        OREFL = GREFL
        OLECL = GLECL
        GETMODEL = LS / LP

      else if ( DTYPE == 3 ) then
        if ( V(7) > 5.0d0 ) then
          ECC = V(7) - 10.0d0
          OMEGA = V(8)
        else
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
        end if
        GETMODEL = GETMIN (V(20),V(19),ECC,OMEGA,TIME)

      else if ( DTYPE == 4 ) then
        GETMODEL = V(15)

      else if ( DTYPE == 5 .or. DTYPE == 6 ) then
        if ( V(7) > 5.0d0 ) then
          ECC = V(7)
          OMEGA = V(8)
          ECOSW = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
          ESINW = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
          if ( DTYPE == 5 )  GETMODEL = V(7)
!ECC
          if ( DTYPE == 6 )  GETMODEL = V(8)
!OMEGA
        else
          ECOSW = V(7)
          ESINW = V(8)
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
          if ( OMEGA > 360.0d0 ) OMEGA = OMEGA - 360.0d0
          if ( DTYPE == 5 )  GETMODEL = ECOSW
          if ( DTYPE == 6 )  GETMODEL = ESINW
        end if

      else
        GETMODEL = -100.0d0
        write(6,*) "### ERROR: wrong datatype asked for in GETMODEL: ",
     &              DTYPE
        STOP
      end if

      END FUNCTION GETMODEL
!=======================================================================
!=======================================================================
      DOUBLEPRECISION FUNCTION GETPHASE (HJD,PERIOD,TZERO)
           
! Returns phase from given time and orbital ephemeris
      implicit none
      real*8 HJD,PERIOD,TZERO

      GETPHASE = (HJD - TZERO) / PERIOD
      GETPHASE = GETPHASE - int(GETPHASE)
      if ( GETPHASE < 0.0d0 ) GETPHASE = GETPHASE + 1.0d0

      END FUNCTION GETPHASE
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMIN (TZERO,PERIOD,ECCIN,OMEGAIN,CICLE)
           
! Returns time of minimum for given cycle and ephemeris.  If
           
! the orbit's circular then the cycle number is used without
           
! restriction so can refer to any phase. If the orbit is ec-
           
! centric then the cycle number should be integer (indicates
           
! primary minimum) or half-integer (secondary minimum).
      implicit none
      real*8 TZERO,PERIOD          
! IN: reference time, orbital period
      real*8 ECCIN,OMEGAIN         
! IN: orbital (e,w) or (ecosw,esinw)
      real*8 CICLE                 
! IN: cycle number of minimum to use
      real*8 CICLEFRAC             
! LOCAL: fraction part of cycle nmbr
      real*8 ECC,OMEGA             
! LOCAL: eccentricity and peri.long.
      real*8 ECOSW,ESINW           
! LOCAL: eccentr'y combination terms
      real*8 PSEP,PHASES(4)        
! LOCAL: phase sep and useful phases
      real*8 PI,DEG2RAD            
! LOCAL: useful variables

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

           
! First must deal with the possibility that e and omega are
           
! actually e*cos(omega) and e*sin(omega)

      if ( ECCIN > 10.0d0 ) then
        ECC = ECCIN - 10.0d0
        OMEGA = OMEGAIN / DEG2RAD
        ECOSW = ECC * cos(OMEGA)
        ESINW = ECC * sin(OMEGA)
      else
        ECC = sqrt(ECCIN**2 + OMEGAIN**2)
        OMEGA = atan2(OMEGAIN,ECCIN)
        ECOSW = ECCIN
        ESINW = OMEGAIN
      end if

           
! If orbit is circular then simply use the orbital ephemeris
           
! If orbit is eccentric then call ECQUADPHASES to calculate
           
! the phase difference between primary and secondary minima.

      if ( abs(ECC) < 1.0d-7 ) then
        GETMIN = TZERO  +  PERIOD * CICLE
      else
        CICLEFRAC = mod(CICLE,1.0d0)

        if ( ( CICLEFRAC >= 0.0d0 .and. CICLEFRAC < 0.001d0 ) .or.
     &       ( CICLEFRAC > 0.999d0 .and. CICLEFRAC <= 1.0d0 ) ) then
          GETMIN = TZERO + PERIOD * CICLE
        else if ((CICLEFRAC > -0.501d0 .and. CICLEFRAC < -0.499d0) .or.
     &           (CICLEFRAC > 0.499d0 .and. CICLEFRAC < 0.501d0) ) then
          CALL ECQUADPHASES (ECCIN,OMEGAIN,0.0d0,PHASES)
          PSEP = PHASES(3) - PHASES(1)
          if ( PSEP < 0.0d0 ) PSEP = PSEP + 1.0d0
          GETMIN = TZERO + PERIOD * (CICLE-0.50d0+PSEP)
        else
          write(6,'(A37,A43)') "### ERROR: found a cycle number which",
     &                    " is not integer or half-integer. Abort.    "
          write(6,'(A18,F20.10)') "### Cycle number =",CICLE
          write(6,*) " "
          stop
        end if

      end if

      END FUNCTION GETMIN
!=======================================================================
!=======================================================================
      SUBROUTINE LIGHT (V,LDTYPE,NSINE,PSINE,NPOLY,PPOLY,HJD,FMAG,LP,LS,
     &                                                   GREFL,GLECL)
      implicit real*8 (a-h,o-z)
      real*8 V(67),HJD,GETPHASE
      real*8 LP,LS,LECL,LE,REFL,GREFL,GLECL
      real*8 LD1U,LD2U           
! linear LD coeff for each star
      real*8 LD1Q,LD2Q           
! quadratic LD coeff for each star
      real*8 LD1S,LD2S           
! square-root LD coeff for each star
      real*8 LD1L,LD2L           
! logarithmic LD coeff for each star
      real*8 LD1C,LD2C           
! cubic LD coeff for each star
      real*8 LDU,LDQ,LDS,LDL,LDC 
! LD coeffs for the star in question
      integer LDTYPE(2)          
! LD law type for both stars
      integer GIMENEZ            
! 1 to use original FMAX calculations
                                 
! 2 to use Gimenez' modified calcs
                                 
! 3 to use Gimenez' nonlinear LD calcs
      integer NSINE,PSINE(5)
      integer NPOLY,PPOLY(5)
      real*8 PHASE,SINT,SINP,SINA,SINTERM,LPMULT,LSMULT

!       data GIMENEZ / 3 /
!       data PI,TWOPI,RAD / 3.1415926536E0,6.28318531E0,0.0174532925E0 /
!       data LPMULT,LSMULT / 1.0 , 1.0 /

      GIMENEZ = 3
      LPMULT = 1.0d0
      LSMULT = 1.0d0
      PI = 3.1415926536d0
      TWOPI = 6.28318531d0
      RAD = 0.0174532925d0


C
C        DETERMINE PRIMARY AND SECONDARY BIAXIAL DIMENSIONS
C        USE SPHERICAL RADII FOR THE COMPUTATION OF ECLIPSE FUNCTIONS
C        USE OBLATENESSES FOR THE COMPUTATION OF THE OUTSIDE ECLIPSE
C        PHOTOMETRIC VARIATIONS WITH LIMB AND GRAVITY DARKENING
C

      if ( V(2) >= 0.0d0 ) then
        RP   = V(2) / ( 1.0d0 + V(3) )
        RS   = V(2) / ( 1.0d0 + (1.0d0/V(3)) )
      else
        RP   = abs( V(2) )
        RS   = V(3)
      end if

      if ( V(7) > 5.0d0 ) then
        ECOSW = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
        ESINW = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
      else
        ECOSW  = V( 7)
        ESINW  = V( 8)
      end if


      BS     = V(1)
      FI     = V(6)
      YP     = V(9)
      YS     = V(10)
      SP     = V(11)
      SS     = V(12)
      Q      = V(13)
      TANGL  = V(14)
      EL     = V(15)
      DPH    = 1.0d0 - V(16)
      SFACT  = V(17)
      DGAM   = V(18)

      LD1U = V(4)              
! linear terms
      LD2U = V(5)
      LD1L = 0.0d0             
! log terms
      LD2L = 0.0d0
      LD1S = 0.0d0             
! sqrt terms
      LD2S = 0.0d0
      LD1Q = 0.0d0             
! quadratic terms
      LD2Q = 0.0d0
      LD1C = 0.0d0             
! cubic terms
      LD2C = 0.0d0
      if ( LDTYPE(1) == 2 ) LD1L = V(21)
      if ( LDTYPE(1) == 3 ) LD1S = V(21)
      if ( LDTYPE(1) == 4 ) LD1Q = V(21)
      if ( LDTYPE(1) == 5 ) LD1C = V(21)
      if ( LDTYPE(2) == 2 ) LD2L = V(22)
      if ( LDTYPE(2) == 3 ) LD2S = V(22)
      if ( LDTYPE(2) == 4 ) LD2Q = V(22)
      if ( LDTYPE(2) == 5 ) LD2C = V(22)
      if ( LDTYPE(2) == 0 ) then
        LD2U = LD1U
        LD2L = LD1L
        LD2S = LD1S
        LD2Q = LD1Q
        LD2C = LD1C
      end if

      if ( NSINE > 0 ) then
        do i = 1,NSINE
          SINT = V(20+i*3)     
! sine reference time
          SINP = V(21+i*3)     
! sine period
          SINA = V(22+i*3)     
! sine amplitude
!           TWOPI =  atan(1.0d0) * 8.0d0
          SINTERM = SINA * sin( TWOPI * (HJD-SINT) / SINP )
          if ( PSINE(i) ==  1 )  BS = BS * (1.0d0+SINTERM)
          if ( PSINE(i) ==  2 )  RP = RP * (1.0d0+SINTERM)
          if ( PSINE(i) ==  3 )  RS = RS * (1.0d0+SINTERM)
          if ( PSINE(i) ==  6 )  FI = FI + SINTERM
          if ( PSINE(i) == 15 )  EL = EL * (1.0d0+SINTERM)
          if ( PSINE(i) == 17 )  SFACT = SFACT + SINTERM
          if ( PSINE(i) == -1 )  LPMULT = LPMULT * (1.0d0+SINTERM)
          if ( PSINE(i) == -2 )  LSMULT = LSMULT * (1.0d0+SINTERM)
        end do
      end if

      if ( NPOLY > 0 ) then
        do i = 1,NPOLY
          j = 32 + (i*6)
          SINTERM = V(j+1)*(HJD-V(j)) + V(j+2)*((HJD-V(j))**2) + V(j+3)*
     &  ((HJD-V(j))**3) + V(j+4)*((HJD-V(j))**4)+ V(j+5)*((HJD-V(j))**5)
          if ( PPOLY(i) ==  1 )  BS = BS + SINTERM  
! Yes this is POLY-
          if ( PPOLY(i) ==  2 )  RP = RP + SINTERM  
! TERM but's called
          if ( PPOLY(i) ==  3 )  RS = RS + SINTERM  
! SINETERM to avoid
          if ( PPOLY(i) ==  6 )  FI = FI + SINTERM   
! proliferation of
          if ( PPOLY(i) == 15 )  EL = EL + SINTERM          
! variables
          if ( PPOLY(i) == 17 )  SFACT = SFACT + SINTERM
          if ( PPOLY(i) == -1 )  LPMULT = LPMULT + SINTERM
          if ( PPOLY(i) == -2 )  LSMULT = LSMULT + SINTERM
        end do
      end if

      PHASE = GETPHASE(HJD,V(19),V(20))

!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!         RS=RP*RATIO
      if ( Q <= 0.0d0 ) then
        CALL BIAX (RP,0.0d0,RPA,RPB,EP)
        CALL BIAX (RS,0.0d0,RSA,RSB,ES)
      else
        CALL BIAX (RP,Q,RPA,RPB,EP)
        CALL BIAX (RS,1.0d0/Q,RSA,RSB,ES)
      end if

C
C        CORRECT THE OBSERVED PHASE FOR ANY EPOCH ERROR IN EPHEMERIS
C
      THETA=PHASE+DPH
C
      SINI  = SIN(FI*RAD)
      SINI2 = SINI*SINI
      COSI2 = 1.0d0  - SINI2
C
C        TRANSLATE TIDAL LEAD/LAG ANGLE TO RADIANS
      TANGR=TANGL*RAD
C
C     EQUATION 9
C        CONVERT PHASE TO RADIANS
      FMN=THETA*TWOPI
C
C        GET CURRENT VALUES OF E, AND W
      CALL GETEW (ECOSW,ESINW,E,W)
C
C        TEST FOR CIRCULAR ORBIT
      IF (E)   17,20,17
   20 COSVW=COS(FMN)
      SINVW=SIN(FMN)
      RV=1.0d0
      GO TO 25
C
C        SOLUTION OF KEPLER'S EQUATION BY DIFFERENTIAL CORRECTIONS
C        (NON-ZERO ECCENTRICITY ONLY . . . )
C
C     EQUATION 6
C
   17 OMEGA = 450.0d0  - W
   23 IF (OMEGA - 360.0d0)         22,21,21
   21 OMEGA = OMEGA - 360.0d0
      GO TO 23
   22 OMEGA = OMEGA*RAD
C        SINE AND COSINE OF OMEGA
      COSW=COS(OMEGA)
      SINW=SIN(OMEGA)
C
C        COMPUTE MEAN ANOMALY CORRECTION TO PHASE
C        CORRESPONDING TO V=OMEGA=90-W
C        AT WHICH PHASE COS(V-OMEGA)=1
      E0=ATAN2(SQRT(1.0d0-E*E)*SINW,COSW+E)
C
C        MEAN ANOMALY OF MID-PRIMARY ECLIPSE
      FMA0=E0-E*SIN(E0)
C
C        MEAN ANOMALY
      FMA=FMN+FMA0
C     FIRST APPROXIMATION OF ECCENTRIC ANOMALY
      EA=FMA+E*SIN(FMA)
C
      DO 10 J=1,15
C        EVALUATE SINE AND COSINE OF ECCENTRIC ANOMALY
      SINE=SIN(EA)
      COSE=COS(EA)
      DENOM=1.0d0-E*COSE
      DISC=FMA-EA+E*SINE
      EA=EA+DISC/DENOM
C        TEST FOR CONVERGENCE
      IF (ABS(DISC) - 2.0d-5)     15,15,10
   10 CONTINUE
C
C
C        EVALUATE SINE AND COSINE OF TRUE ANOMALY
   15 COSV=(COSE-E)/DENOM
      SINV=SINE*SQRT(1.0d0-E*E)/DENOM
C
C        RADIUS VECTOR
      RV = (1.0d0-E*E)/(1.0d0+E*COSV)
C
C        THE PHOTOMETRIC PHASE ARGUMENT IN TERMS OF ORBIT PARAMETERS
C        VW = V-OMEGA
      COSVW=COSV*COSW+SINV*SINW
      SINVW=SINV*COSW-COSV*SINW
C
   25 COS2=COSVW*COSVW
      SIN2=1.0d0-COS2
C
      CSVWT=COS(TANGR)*COSVW-SIN(TANGR)*SINVW
C
C
C        PHOTOMETRIC EFFECTS
C
C

      FMAXP = 0.0d0
      FMAXS = 0.0d0
      DELTP = 0.0d0
      DELTS = 0.0d0
      SHORT = 0.0d0

!-----------------------------------------------------------------------
! Alvaro Gimenez and J Diaz-Cordoves have corrected the treatment of LD
! and stellar shapes.  This treatment can be used by putting GIMENEZ=2
! Their treatment for nonlinear LD can be used by putting GIMENEZ=3
!-----------------------------------------------------------------------
! This whole thing affects only the brightness normalisation of the two
! eclipsing stars: any problems here affect the radiative parameters
! but not the geometric parameters (radii, inclination etc).
!-----------------------------------------------------------------------

      if ( GIMENEZ==1 ) then                         
! LINEAR LD ONLY

!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))
! Original
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)               
! lines
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))
! if the
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)               
! stars
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP   
! are
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES   
! oblate
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0                                
! Original
!       FMAXS=1.0E0-US/3.0E0                                
! lines if
!       DELTP=0.0E0                                         
! the stars
!       DELTS=0.0E0                                         
! are
!       SHORT=0.0                                           
! spherical

        if ( Q >= 0.0d0 ) then
          FMAXP=((1.0d0-LD1U)+0.666666667d0*LD1U*(1.0d0+0.2d0*EP))
     1        *(1.0d0+3.0d0*YP*EP)/(1.0d0-EP)
          FMAXS=((1.0d0-LD2U)+0.666666667d0*LD2U*(1.0d0+0.2d0*ES))
     1        *(1.0d0+3.0d0*YS*ES)/(1.0d0-ES)
          DELTP=(15.0d0+LD1U)/(15.0d0-5.0d0*LD1U)*(1.0d0+YP)*EP
          DELTS=(15.0d0+LD2U)/(15.0d0-5.0d0*LD2U)*(1.0d0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0d0-LD1U/3.0d0
          FMAXS=1.0d0-LD2U/3.0d0
          DELTP=0.0d0
          DELTS=0.0d0
          SHORT=0.0
        end if
!-----------------------------------------------------------------------

      else if ( GIMENEZ==2 ) then                    
! LINEAR LD ONLY

!       FMAXP=(1.0E0-UP*(1.0E0-2.0E0/5.0E0*EP)/3.0E0+YP*EP  
! Original
!      1      *(3.0E0-13.0E0/15.0E0*UP))/(1.0E0-EP)         
! lines
!       FMAXS=(1.0E0-US*(1.0E0-2.0E0/5.0E0*ES)/3.0E0+YS*ES  
! if the
!      1      *(3.0E0-13.0E0/15.0E0*US))/(1.0E0-ES)         
! stars
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP   
! are
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES   
! oblate
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0                                
! Original
!       FMAXS=1.0E0-US/3.0E0                                
! lines if
!       DELTP=0.0E0                                         
! the stars
!       DELTS=0.0E0                                         
! are
!       SHORT=0.0                                           
! spherical

        if ( Q >= 0.0d0 ) then
          FMAXP=(1.0d0-LD1U*(1.0d0-2.0d0/5.0d0*EP)/3.0d0+YP*EP
     1          *(3.0d0-13.0d0/15.0d0*LD1U))/(1.0d0-EP)
          FMAXS=(1.0d0-LD2U*(1.0d0-2.0d0/5.0d0*ES)/3.0d0+YS*ES
     1          *(3.0d0-13.0d0/15.0d0*LD2U))/(1.0d0-ES)
          DELTP=(15.0d0+LD1U)/(15.0d0-5.0d0*LD1U)*(1.0d0+YP)*EP
          DELTS=(15.0d0+LD2U)/(15.0d0-5.0d0*LD2U)*(1.0d0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0d0-LD1U/3.0d0
          FMAXS=1.0d0-LD2U/3.0d0
          DELTP=0.0d0
          DELTS=0.0d0
          SHORT=0.0d0
        end if
!-----------------------------------------------------------------------
! And this is Gimenez's code for including nonlinear LD. He includes
! the linear (UP), quadratic (UP, U2P) and square-root (UP, U3P) laws.
!-----------------------------------------------------------------------

      else if ( GIMENEZ==3 ) then

!      FMAXP=1.0E0-UP*(1.0E0-2.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0-3.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0-4.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.5E0-13.0E0*UP/30.0E0-U2P/5.0E0-23.0E0*U3P/90.0E0)
!      FMAXP=FMAXP/(1.0E0-EP)
!      FMINP=1.0E0-UP*(1.0E0+4.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0+6.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0+8.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.0E0-7.0E0*UP/15.0E0-4.0E0*U2P/15.0E0-13.0E0*U3P/45.0E0)
!      FMINS=1.0E0-US*(1.0E0+4.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0+6.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0+8.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.0E0-7.0E0*US/15.0E0-4.0E0*U2S/15.0E0-13.0E0*U3S/45.0E0)
!      FMAXS=1.0E0-US*(1.0E0-2.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0-3.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0-4.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.5E0-13.0E0*US/30.0E0-U2S/5.0E0-23.0E0*U3S/90.0E0)
!      FMAXS=FMAXS/(1.0E0-ES)
!      DELTP=1.0E0-FMINP/FMAXP
!      DELTS=1.0E0-FMINS/FMAXS
!      SHORT=SINI2*CSVWT*CSVWT

!   26 FMAXP=1.0E0-UP/3.0E0-U2P/6.0E0-U3P/5.0E0
!      FMAXS=1.0E0-US/3.0E0-U2S/6.0E0-U3S/5.0E0
!      DELTP=0.0E0
!      DELTS=0.0E0
!      SHORT=0.0

        if ( Q >= 0.0d0 .or. LDTYPE(1)==1 .or. LDTYPE(1)==5 .or.
     &                         LDTYPE(2)==1 .or. LDTYPE(2)==5 ) then
          FMAXP=1.0d0-LD1U*(1.0d0-2.0d0*EP/5.0d0)/3.0d0-
     &          LD1Q*(1.0d0-3.0d0*EP/5.0d0)/6.0d0-
     &          LD1S*(1.0d0-4.0d0*EP/9.0d0)/5.0d0+2.0d0*YP*EP
     &         *(1.5E0-13.0d0*LD1U/30.0d0-LD1Q/5.0d0-23.0d0*LD1S/90.0d0)
          FMAXP=FMAXP/(1.0d0-EP)
          FMINP=1.0d0-LD1U*(1.0d0+4.0d0*EP/5.0d0)/3.0d0-
     &          LD1Q*(1.0d0+6.0d0*EP/5.0d0)/6.0d0-
     &          LD1S*(1.0d0+8.0d0*EP/9.0d0)/5.0d0+2.0d0*YP*EP
     &   *(1.0d0-7.0d0*LD1U/15.0d0-4.0d0*LD1Q/15.0d0-13.0d0*LD1S/45.0d0)
          FMINS=1.0d0-LD2U*(1.0d0+4.0d0*ES/5.0d0)/3.0d0-
     &          LD2Q*(1.0d0+6.0d0*ES/5.0d0)/6.0d0-
     &          LD2S*(1.0d0+8.0d0*ES/9.0d0)/5.0d0+2.0d0*YS*ES
     &   *(1.0d0-7.0d0*LD2U/15.0d0-4.0d0*LD2Q/15.0d0-13.0d0*LD2S/45.0d0)
          FMAXS=1.0d0-LD2U*(1.0d0-2.0d0*ES/5.0d0)/3.0d0-
     &          LD2Q*(1.0d0-3.0d0*ES/5.0d0)/6.0d0-
     &          LD2S*(1.0d0-4.0d0*ES/9.0d0)/5.0d0+2.0d0*YS*ES
     &         *(1.5E0-13.0d0*LD2U/30.0d0-LD2Q/5.0d0-23.0d0*LD2S/90.0d0)
          FMAXS=FMAXS/(1.0d0-ES)
          DELTP=1.0d0-FMINP/FMAXP
          DELTS=1.0d0-FMINS/FMAXS
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0-LD1U/3.0-LD1Q/6.0-LD1S/5.0+LD1L*2.0/9.0-LD1C/10.0
          FMAXS=1.0-LD2U/3.0-LD2Q/6.0-LD2S/5.0+LD2L*2.0/9.0-LD2C/10.0
          DELTP=0.0d0
          DELTS=0.0d0
          SHORT=0.0d0
        end if
!----------------------------------------------------------------------
      end if
!----------------------------------------------------------------------
! Complete original code before the above messing:
! C
! C
! C        PHOTOMETRIC EFFECTS
! C
! C
! C        TEST FOR SIMPLE CASE OF TWO SPHERICAL STARS
!       IF (EP .EQ. 0.  .AND.  ES .EQ. 0.)   GO TO 26
! C
! C        EITHER OR BOTH STARS ARE OBLATE
! C
!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
! C        CHANGE IN INTENSITY RATIO DUE TO OBLATENESS RELATED VARIABLES
! C        FROM QUADRATURE TO MINIMUM
! C        FACE ON TO END ON
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
! C        FORE-SHORTENING FUNCTION OF OBLATENESS
!       SHORT=SINI2*CSVWT*CSVWT
!       GO TO 27
! C
! C        BOTH STARS ARE SPHERICAL
! C
!    26 FMAXP=1.0E0-UP/3.0E0
!       FMAXS=1.0E0-US/3.0E0
!       DELTP=0.0E0
!       DELTS=0.0E0
!       SHORT=0.0
!----------------------------------------------------------------------

C
C        UN-NORMALIZED BRIGHTNESS OF STELLAR COMPONENTS AT QUADRATURE
   27 OP=PI*RPB*RPB*FMAXP
      OS=PI*RSB*RSB*FMAXS*BS
C        THE NORMALIZING FACTOR
      OTOT=OP+OS
C        BRIGHTNESS CONTRIBUTION FROM EACH COMPONENT
      LP=OP/OTOT*(1.0d0-DELTP*SHORT)
      LS=OS/OTOT*(1.0d0-DELTS*SHORT)
C
C        REFLECTION AND RERADIATION EQUATION
      IF (SP .EQ. 0.0d0  .AND.  SS .EQ. 0.0d0)   GO TO 28
      HEAT=SINI*COSVW
      HEAT2=0.5d0+0.5d0*HEAT*HEAT
      DLP=SP*(HEAT2+HEAT)
      DLS=SS*(HEAT2-HEAT)
      GO TO 29
   28 DLP=0.0d0
      DLS=0.0d0
C
C        WHICH ECLIPSE COULD THIS BE
   29 IF (COSVW)         40,40,30
C
C     PRIMARY ECLIPSE
C
   30 R1 = RP
      R2 = RS
!---------------------------------------------------------------------
!
! JKT mod (10/8/2006): the line these replaced was      UU = UP       
!
!---------------------------------------------------------------------
!
      LDU = LD1U                                                      
!
      LDL = LD1L                                                      
!
      LDS = LD1S                                                      
!
      LDQ = LD1Q                                                      
!
      LDC = LD1C                                                      
!
!---------------------------------------------------------------------
!
      LE=LP
      DLE=DLP
      GO TO 60
C
C
C     SECONDARY ECLIPSE
C
   40 R1 = RS
      R2 = RP
!-----------------------------------------------------------------------
! JKT mod (10/8/2006): the line these replaced was      UU = US       
!
!---------------------------------------------------------------------
!
      LDU = LD2U                                                      
!
      LDL = LD2L                                                      
!
      LDS = LD2S                                                      
!
      LDQ = LD2Q                                                      
!
      LDC = LD2C                                                      
!
!---------------------------------------------------------------------
!
      LE=LS
      DLE=DLS
C
   60 SUM = 0.0d0
      ALAST = 0.0d0
      AREA=0.0d0
C
C     EQUATION  5
C
      DD = SINVW*SINVW + COSVW*COSVW*COSI2
C      IF (DD .LE. 1.0d-6)  DD=0.0d0
      DD = DD*RV*RV
      D = SQRT(ABS(DD))
      R22 = R2*R2
C
C     EQUATION 17
C
      GAMN = 90.01d0*RAD
      DGAMA = DGAM*RAD
      DGM = DGAMA/2.0d0
      RK = 0.0d0
      GAM = 0.0d0
   50 GAM = GAM + DGAMA
C        HAS LIMIT OF INTEGRATION BEEN REACHED
      IF (GAM - GAMN)              48,48,49
C
   48 RR = R1*SIN(GAM)
      R12 = RR*RR
C
      AA = 0.0d0
C        ARE THE PROJECTED DISKS CONCENTRIC
      IF (D)                       405,406,405
  406 IF (RR - R2)                 230,230,403
  403 IF (RK - R2)                 404, 49, 49
  404 AA = PI*R22
      GO TO 215
C        TEST FOR NO ECLIPSE
  405 IF (D-R1-R2)                 240,216,216
  216 SUM = 0.0d0
      GO TO 49
C        DECIDE WHICH AREA EQUATIONS FOR NON-CONCENTRIC ECLIPSE
  240 IF (D-RR-R2)                 245,215,215
  245 IF (D-R2+RR)                 230,230,250
  250 IF (R1-R2)                   255,255,280
  255 IF (DD-R22+R12)              205,210,210
  280 IF (D-RR+R2)                 290,260,260
  260 IF (RR-R2)                   255,255,265
  265 IF (DD-R12+R22)              270,210,210
C
C     EQUATION 12
C
  270 S1 = ABS((R12 - R22 - DD)*0.5d0/D)
      A1 = ABS(R2-S1)
      B2 = ABS(RR-S1-D  )
      AA=PI*R22-(R22*ACOS((R2-A1)/R2)
     1   - (R2-A1)*SQRT(2.0d0*R2*A1-A1*A1))
     2   +R12*ACOS((RR-B2)/RR)-(RR-B2)*SQRT(2.0d0*RR*B2-B2*B2)
      GO TO 215
C
  290 IF (R1 - R2 - D)             260,260,295
  295 IF (RK - R2 - D)             300,215,215
  300 RR = R2 + D
      R12 = RR*RR
      GAMN = 0.0d0
      GO TO 260
C
  230 AA = PI*R12
      GO TO 215
C
C     EQUATION 10
C
  205 S = ABS((R12 - R22 + DD)*0.5d0/D)
      A = ABS(RR-S)
      B1 = ABS(R2-S-D)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0d0*RR*A - A*A)
      AB1 = R22*ACOS((R2-B1)/R2) - (R2-B1)*SQRT(2.0d0*R2*B1-B1*B1)
      AA = PI*R12 - A1 + AB1
      GO TO 215
C
C     EQUATION 1
C
  210 S = ABS((R12 - R22 + DD)*0.5d0/D)
      A = ABS(RR-S)
      B = ABS(S-D+R2)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0d0*RR*A - A*A)
      AA1 = R22*ACOS((R2-B)/R2) - (R2-B)*SQRT(2.0d0*R2*B - B*B)
      AA = A1 + AA1
C
  215 DAREA = AA - ALAST
!---------------------------------------------------------------------
!
! JKT modification (10/9/2006). The removed line was:                 
!
!     SUM = SUM + DAREA*(1.0d0  - UU + UU*COS(GAM-DGM))               
!
!---------------------------------------------------------------------
!
      COSGAM = cos(GAM-DGM)                                           
!
      SUM = SUM + DAREA*(1.0d0 - LDU*(1.0d0-COSGAM)                   
!
     &          - LDL*COSGAM*log(COSGAM) - LDS*(1.0d0-sqrt(COSGAM))   
!
     &         - LDQ*(1.0d0-COSGAM)**2 - LDC*(1.0d0-COSGAM)**3)       
!
!---------------------------------------------------------------------
!
      ALAST = AA
      AREA = AREA + DAREA
C
      RK = RR
      GO TO 50
C
C        LIGHT LOSS FROM ECLIPSE
C
   49 ADISK = PI*R1*R1
!---------------------------------------------------------------------
!
! JKT modification (10/9/2006).  See 1992A+A...259..227D for more info
!
! The removed line was:           ALPHA = SUM/(ADISK*(1.0-UU/3.0))    
!
!---------------------------------------------------------------------
!
      ALPHA = 1.0d0 - LDU/3.0d0 + LDL*2.0d0/9.0d0 -                   
!
     &          LDS/5.0d0 - LDQ/6.0d0 - LDC/10.0d0                    
!
      ALPHA = SUM/(ADISK*ALPHA)                                       
!
!---------------------------------------------------------------------
!
      LECL = ALPHA*LE
      AREA = AREA/ADISK
      REFL=DLP+DLS-AREA*DLE
C
C        THEORETICAL INTENSITY WITH THIRD LIGHT AND QUADRATURE
C        SCALE FACTOR APPLIED
C
!---------------------------------------------------------------------
!
! This is the original line from EBOP:
!---------------------------------------------------------------------
!
!      FLITE = ((LP+LS-LECL+REFL)*(1.0d0-EL)+EL)*SFACT
!---------------------------------------------------------------------
!

      LP = LP * LPMULT              
! sine/poly applied to star 1 light
      LS = LS * LSMULT              
! sine/poly applied to star 2 light
C      FLITE = ((LP+LS-LECL+REFL)*(1.0d0-EL)+EL)
      FLITE = (LP+LS-LECL+REFL)
C     R.F. DIAZ : return flux instead of mag
C      FMAG = -2.5d0 * log10(FLITE) + SFACT
      FMAG = FLITE

C     Lines commented to return undiluted LP and LS
C      LP = LP * (1.0d0-EL)          
! account for third light *AFTER*
C      LS = LS * (1.0d0-EL)          

! FLITE and FMAG have been found
C     J.M. Almenara : return LP, LS (take into account 
C     the 3rd light), REFL and GLECL
      
      GREFL = REFL
      GLECL = LECL    

      END
!=======================================================================
      SUBROUTINE BIAX (R,Q,A,B,EPS)
           
! EBOP subroutine to calculate biaxial ellipsoid dimensions
           
! and oblateness for each star after Chandrasekhar (1933).
      implicit none
      real*8 R,Q,A,B,EPS

      if ( Q <= 0.0d0 )  then
        A = R
        B = R
        EPS = 0.0d0
      else
        A = R * ( 1.0d0 + (1.0d0 + 7.0d0*Q)/6.0d0 * R**3.0d0)
        B = R * ( 1.0d0 + (1.0d0 - 2.0d0*Q)/6.0d0 * R**3.0d0)
        EPS = (A - B) / A
        B = ( (1.0d0 - EPS) * R**3.0d0) ** (1.0d0/3.0d0)
        A = B / (1.0d0 - EPS)
      end if

      END SUBROUTINE BIAX
!=======================================================================
      SUBROUTINE GETEW (ECOSW,ESINW,E,W)
           
! EBOP subroutine to calculate e and w from e(cos)w e(sin)w
      implicit none
      real*8 ECOSW,ESINW,E,W

      if ( ECOSW == 0.0d0  .and.  ESINW == 0.0d0 ) then
        E = 0.0d0
        W = 0.0d0
      else
        W = atan2( ESINW,ECOSW )
        E = sqrt( ESINW*ESINW + ECOSW*ECOSW )
        W = W * 180.0d0 / 3.1415926536d0
      end if

      END SUBROUTINE GETEW
!=======================================================================
      SUBROUTINE ECQUADPHASES (ECCIN,OMEGAIN,PSHIFT,PHASES)
           
! Calculates orbital phases of the eclipses and quadratures.
           
! PHASES(1) and PHASES(3) are the prim and sec eclipse times
           
! PHASES(2) and PHASES(4) are the *photometric( quadratures.
      implicit none
      real*8 ECCIN,OMEGAIN         
! IN: eccentricity and peri.long.
      real*8 PSHIFT                
! IN: time of primary eclipse
      real*8 PHASES(4)             
! OUT: phases of eclipses and quads
      real*8 PI,DEG2RAD            
! LOCAL: constants
      real*8 ECC,OMEGA             
! LOCAL: values of ecc and peri.long
      real*8 ECOSW,ESINW           
! LOCAL: values of combination terms
      real*8 EFAC,EFAC2            
! LOCAL: useful eccentricity factors
      real*8 TERM1,TERM2           
! LOCAL: calculation helper varibles
      real*8 PHASEDIFF             
! LOCAL: diff of prim and sec minima

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)

           
! Get actual eccentricity and periastron longitude values
           
! from the input, which could be (e+10,w) or (ecosw,esinw)

      if ( ECCIN > 10.0d0 ) then
        ECC = ECCIN - 10.0d0
        OMEGA = OMEGAIN / DEG2RAD
        ECOSW = ECC * cos(OMEGA)
        ESINW = ECC * sin(OMEGA)
      else
        ECC = sqrt(ECCIN**2 + OMEGAIN**2)
        OMEGA = atan2(OMEGAIN,ECCIN)
        ECOSW = ECCIN
        ESINW = OMEGAIN
      end if

      EFAC = sqrt( (1.0d0-ECC) / (1.0d0+ECC) )
      EFAC2 = sqrt(1.0d0 - ECC**2)

!     write(6,*)" "
!     write(6,'(a16,2(f13.8))')"pi,deg2rad      ",pi,deg2rad
!     write(6,'(a16,2(f13.8))')"eccin,omegain   ",eccin,omegain
!     write(6,'(a16,2(f13.8))')"ecc,omega       ",ecc,omega
!     write(6,'(a16,2(f13.8))')"ecosw,esinw     ",ecosw,esinw
!     write(6,'(a16,2(f13.8))')"efac,efac2      ",efac,efac2
!     write(6,'(a16,2(f13.8))')"pshift,omegadeg ",pshift,omega*deg2rad

! The equation for phase difference comes from Hilditch (2001) page 238
! equation 5.66, originally credited to the monograph by  Kopal (1959).

      TERM1 = 2.0d0 * atan( ECOSW / EFAC2 )
      TERM2 = 2.0d0 * ECOSW * EFAC2 / (1.0d0 - ESINW**2)
      PHASEDIFF = ( PI + TERM1 + TERM2 ) / ( 2.0d0 * PI )

      PHASES(1) = PSHIFT
      PHASES(2) = PSHIFT + PHASEDIFF/2.0d0
      PHASES(3) = PSHIFT + PHASEDIFF
      PHASES(4) = PSHIFT + 0.50d0 + PHASEDIFF/2.0d0

!     write(6,*)" "
!     write(6,'(A19,2(F14.8))')"term1,term2        ",term1,term2
!     write(6,'(A19,2(F14.8))')"timediff,phasediff ",phasediff
!     write(6,'(A19,2(F14.8))')"eclipse phases     ",phases(1),phases(3)
!     write(6,'(A19,2(F14.8))')"quadrature phases  ",phases(2),phases(4)

      END SUBROUTINE ECQUADPHASES
!=======================================================================
