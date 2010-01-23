
!------------------------------------------------------------------------------
!
!  This code is used to determine the rigid-rotor harmonic oscillator 
!  frequencies and hindered internal rotor barrier heights (and 
!  their corresponding low-temoerature frequencies) for a given molecule.
!
!  It requires as input a file with the following format:
!       Cp(300) Cp(400) Cp(500) Cp(600) Cp(800) Cp(1000) Cp(1500)
!       Number of atoms
!       Number of internal rotors
!       Linearity (0 if linear, 1 if non-linear)
!       List of the number of each bond type, described below.
!
!  The frequency estimation takes place in two parts.
!  In the first part, the code looks at the number and type of bonds 
!  available and generates a list of RRHO frequencies that are typical
!  of those bonds.  These frequencies are stored in the vector:
!
!  Total_predicted_freq
!
!  And they are determined by calling the subroutine:
!
!  calc_predicted_freq
!
!  Once the code has returned the vector of predicted frequencies,
!  the contribution of these frequencies towards the heat capacity is
!  calculated, and this contribution is then subtracted from the total
!  heat capacity (which is supplied in the first line of the input file).
!  What remains is the heat capacity due to the unkown active degrees
!  (e.g. all internal rotors, as well as any RRHO frequencies not 
!  predicted by calc_predicted_freq).
!  The second part of the code takes the remaining heat capacity and uses
!  it to determine the barrier heights and unknown frequencies.  This is
!  accomplished by calling:
!
!  assign_cases
!
!  This subroutine will look at the number of unknown RRHO frequencies and 
!  the number of rotors, and it will select one of 26 possible cases to
!  fit these RRHO and HR parameters to the heat capacity data.

!  The fitting is done using the code DQED, available on netlib,
!  which is a bounded, contrained least squares solver for non-linear 
!  equations.  The DQED code must be supplied as a separate file, 
!  dqed.f90, at the time of compilation.
!
!
!  The vector Total_predicted_freq is copies to the vector
!  Total_harm_osc_freq, and the unknown RRHO frequencies determined by the 
!  fitting procedure are appended to the end of the vector.
!  The hindered-rotor parameters are stored in an array HR_params.
!
!---------------------------------------------------------------------------  
!
!  The code is broken into several parts.
!
!  The code begins with a large module heat_capacity_functions.
!  This module provides the heat capacity functions and their Jacobians
!
!  The first 500 lines or so of this module were cut and pasted from Argonne.
!  The Argonne code is used to calculated the modified Bessel functions 
!  Of the first kind and orders zero and one:  I_0 and I_1.
!  
!  Following the Argonne Code are four new functions.
!
!  The first two functions are: cv_hind_rot and d_cv_hind_rot
!
!  cv_hind_rot takes as input a barrier height, CV_temps, and harmonic 
!  Oscillator frequency and returns the corresponding heat capacity
!  d_cv_hind_rot is the derivative of cv_hind_rot with respect to the barrier
!  Height.  This function is required for the analytic Jacobian.
!
!  The final two functions are:  cv_harm_osc and d_cv_harm_osc
!
!  cv_harm_osc takes as input the harmonic oscillator frequency and CV_temps
!  And returns the corresponding Cv value.
!
!  d_cv_harm_osc is the derivative of cv_harm_osc with respect to the frequency
!  Which is required for the analytic Jacobian
!
!  The module heat_capacity_functions ends with d_cv_harm_osc.
!
!  The next module is frequencies.  It contains all the subroutines needed 
!  to estimate the RRHO frequencies based on the input structure.
!
!  The third module is cases.  It contains the 26 subroutines used in the
!  fitting procedure.
!
!  The fourth module is open_stuff.  It contains the subroutines that are used 
!  to read the input file and extract the relevant info.
!  Most importantly, it contains the main subroutine
!
!  Calculate_RRHO_HR_params
!
!  Which is the master program.  This subroutine will be called by the main
!  program.  It requires as input the name of the file containing the
!  input data and the name of the file to which it should write the
!  frequencies.  The format of the output file is designed to be read by the 
!  separate code, calc_rho.f90, which uses this info to calculate the 
!  density of states.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  BEGIN THE MODULE FOR CALCULATING HEAT CAPACITIES
!------------------------------------------------------------------------------

MODULE heat_capacity_functions

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

  PRIVATE
  PUBLIC :: calci0, besi0, calci1, besi1
  PUBLIC :: cv_hind_rot, d_cv_hind_rot, d_cv_hind_rot_nu
  PUBLIC :: cv_harm_osc, d_cv_harm_osc

  CONTAINS

! The I0 Bessel functions

SUBROUTINE calci0(arg, result, jint)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!--------------------------------------------------------------------

! This packet computes modified Bessel functions of the first kind
!   and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
!   arguments X.  It contains two function type subprograms, BESI0
!   and BESEI0, and one subroutine type subprogram, CALCI0.
!   The calling statements for the primary entries are

!                   Y=BESI0(X)
!   and
!                   Y=BESEI0(X)

!   where the entry points correspond to the functions I0(X) and
!   EXP(-ABS(X))*I0(X), respectively.  The routine CALCI0 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCI0 with the statement
!          CALL CALCI0(ARG,RESULT,JINT)
!   where the parameter usage is as follows

!      Function                     Parameters for CALCI0
!       Call              ARG                  RESULT          JINT

!     BESI0(ARG)    ABS(ARG) .LE. XMAX        I0(ARG)           1
!     BESEI0(ARG)    any real ARG        EXP(-ABS(ARG))*I0(ARG) 2

!   The main computation evaluates slightly modified forms of
!   minimax approximations generated by Blair and Edwards, Chalk
!   River (Atomic Energy of Canada Limited) Report AECL-4928,
!   October, 1974.  This transportable program is patterned after
!   the machine-dependent FUNPACK packet NATSI0, but cannot match
!   that version for efficiency or accuracy.  This version uses
!   rational functions that theoretically approximate I-SUB-0(X)
!   to at least 18 significant decimal digits.  The accuracy
!   achieved depends on the arithmetic system, the compiler, the
!   intrinsic functions, and proper selection of the machine-
!   dependent constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   beta   = Radix for the floating-point system
!   maxexp = Smallest power of beta that overflows

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XSMALL = Positive argument such that 1.0 - X = 1.0 to
!            machine precision for all ABS(X) .LE. XSMALL.
!   XINF =   Largest positive machine number; approximately beta**maxexp
!   XMAX =   Largest argument acceptable to BESI0;  Solution to
!            equation:
!               W(X) * (1+1/(8*X)+9/(128*X**2) = beta**maxexp
!            where  W(X) = EXP(X)/SQRT(2*PI*X)


!     Approximate values for some important machines are:

!                          beta       maxexp       XSMALL

! CRAY-1        (S.P.)       2         8191       3.55E-15
! Cyber 180/855
!   under NOS   (S.P.)       2         1070       3.55E-15
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)       2          128       2.98E-8
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)       2         1024       5.55D-17
! IBM 3033      (D.P.)      16           63       6.95D-18
! VAX           (S.P.)       2          127       2.98E-8
! VAX D-Format  (D.P.)       2          127       6.95D-18
! VAX G-Format  (D.P.)       2         1023       5.55D-17


!                               XINF          XMAX

! CRAY-1        (S.P.)       5.45E+2465     5682.810
! Cyber 180/855
!   under NOS   (S.P.)       1.26E+322       745.893
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)       3.40E+38         91.900
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)       1.79D+308       713.986
! IBM 3033      (D.P.)       7.23D+75        178.182
! VAX           (S.P.)       1.70D+38         91.203
! VAX D-Format  (D.P.)       1.70D+38         91.203
! VAX G-Format  (D.P.)       8.98D+307       713.293

!*******************************************************************

! Error returns

!  The program returns XINF for BESI0 for ABS(ARG) .GT. XMAX.


!  Intrinsic functions required are:

!     ABS, SQRT, EXP


!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439

!  Latest modification: March 12, 1992

!--------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (8), INTENT(IN)   :: arg
REAL (8), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: jint

! Local variables

INTEGER    :: i
REAL (8)  :: a, b, sump, sumq, x, xx
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: one = 1.0_dp, one5 = 15.0_dp, forty = 40.0_dp,  &
                         exp40 = 2.353852668370199854D17, two25 = 225.0_dp, &
                         rec15 = 6.6666666666666666666D-2
!--------------------------------------------------------------------
!  Machine-dependent constants
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: XSMALL = 5.55D-17, XINF = 1.79D308, XMAX = 713.986D0
!--------------------------------------------------------------------
!  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: p(15) =  &
        (/ -5.2487866627945699800D-18, -1.5982226675653184646D-14,  &
           -2.6843448573468483278D-11, -3.0517226450451067446D-08,  &
           -2.5172644670688975051D-05, -1.5453977791786851041D-02,  &
           -7.0935347449210549190D+00, -2.4125195876041896775D+03,  &
           -5.9545626019847898221D+05, -1.0313066708737980747D+08,  &
           -1.1912746104985237192D+10, -8.4925101247114157499D+11,  &
           -3.2940087627407749166D+13, -5.5050369673018427753D+14,  &
           -2.2335582639474375249D+15 /)
REAL (8), PARAMETER  :: Q(5) =  &
        (/ -3.7277560179962773046D+03, 6.5158506418655165707D+06,  &
           -6.5626560740833869295D+09, 3.7604188704092954661D+12,  &
           -9.7087946179594019126D+14 /)
!--------------------------------------------------------------------
!  Coefficients for 15.0 .LE. ABS(ARG)
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: PP(8) =  &
        (/ -3.9843750000000000000D-01, 2.9205384596336793945D+00,  &
           -2.4708469169133954315D+00, 4.7914889422856814203D-01,  &
           -3.7384991926068969150D-03,-2.6801520353328635310D-03,  &
            9.9168777670983678974D-05,-2.1877128189032726730D-06 /)
REAL (8), PARAMETER  :: QQ(7) =  &
        (/ -3.1446690275135491500D+01, 8.5539563258012929600D+01,  &
           -6.0228002066743340583D+01, 1.3982595353892851542D+01,  &
           -1.1151759188741312645D+00, 3.2547697594819615062D-02,  &
           -5.5194330231005480228D-04 /)
!--------------------------------------------------------------------

x = ABS(arg)
IF (x < xsmall) THEN
  result = one
ELSE IF (x < one5) THEN
!--------------------------------------------------------------------
!  XSMALL .LE.  ABS(ARG)  .LT. 15.0
!--------------------------------------------------------------------
  xx = x * x
  sump = p(1)
  DO  i = 2, 15
    sump = sump * xx + p(i)
  END DO
  xx = xx - two25
  sumq = ((((xx+q(1))*xx+q(2))*xx+q(3))*xx+q(4)) * xx + q(5)
  result = sump / sumq
  IF (jint == 2) result = result * EXP(-x)
ELSE IF (x >= one5) THEN
  IF ((jint == 1).AND.(x > xmax)) THEN
    result = xinf
  ELSE
!--------------------------------------------------------------------
!  15.0  .LE.  ABS(ARG)
!--------------------------------------------------------------------
    xx = one / x - rec15
    sump = ((((((pp(1)*xx+pp(2))*xx+pp(3))*xx+pp(4))*xx+pp(5))*xx+  &
        pp(6))*xx+pp(7)) * xx + pp(8)
    sumq = ((((((xx+qq(1))*xx+qq(2))*xx+qq(3))*xx+qq(4))*xx+  &
        qq(5))*xx+qq(6)) * xx + qq(7)
    result = sump / sumq
    IF (jint == 2) THEN
      result = (result-pp(1)) / SQRT(x)
    ELSE
!--------------------------------------------------------------------
!  Calculation reformulated to avoid premature overflow
!--------------------------------------------------------------------
      IF (x <= (xmax-one5)) THEN
        a = EXP(x)
        b = one
      ELSE
        a = EXP(x-forty)
        b = exp40
      END IF
      result = ((result*a-pp(1)*a)/SQRT(x)) * b
    END IF
  END IF
END IF
!--------------------------------------------------------------------
!  Return for ABS(ARG) .LT. XSMALL
!--------------------------------------------------------------------
RETURN
!----------- Last line of CALCI0 -----------
END SUBROUTINE calci0



FUNCTION BESI0(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This long precision subprogram computes approximate values for
!   modified Bessel functions of the first kind of order zero for
!   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI0).

!--------------------------------------------------------------------

REAL (8), INTENT(IN)  :: x
REAL (8)              :: fn_val

! Local variable
INTEGER  :: jint

!--------------------------------------------------------------------
jint = 1
CALL calci0(x, fn_val, jint)
RETURN
!---------- Last line of BESI0 ----------
END FUNCTION BESI0


SUBROUTINE calci1(arg, result, jint)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!--------------------------------------------------------------------

! This packet computes modified Bessel functioons of the first kind
!   and order one, I1(X) and EXP(-ABS(X))*I1(X), for real arguments X.
!   It contains two function type subprograms, BESI1 and BESEI1,
!   and one subroutine type subprogram, CALCI1.
!   The calling statements for the primary entries are

!                   Y=BESI1(X)
!   and
!                   Y=BESEI1(X)

!   where the entry points correspond to the functions I1(X) and
!   EXP(-ABS(X))*I1(X), respectively.  The routine CALCI1 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCI1 with the statement
!          CALL CALCI1(ARG,RESULT,JINT)
!   where the parameter usage is as follows

!      Function                     Parameters for CALCI1
!       Call              ARG                  RESULT          JINT

!     BESI1(ARG)    ABS(ARG) .LE. XMAX        I1(ARG)           1
!     BESEI1(ARG)    any real ARG        EXP(-ABS(ARG))*I1(ARG) 2

!   The main computation evaluates slightly modified forms of minimax
!   approximations generated by Blair and Edwards, Chalk River (Atomic
!   Energy of Canada Limited) Report AECL-4928, October, 1974.
!   This transportable program is patterned after the machine-dependent
!   FUNPACK packet NATSI1, but cannot match that version for efficiency or
!   accuracy.  This version uses rational functions that theoretically
!   approximate I-SUB-1(X) to at least 18 significant decimal digits.
!   The accuracy achieved depends on the arithmetic system, the compiler, the
!   intrinsic functions, and proper selection of the machine-dependent
!   constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   beta   = Radix for the floating-point system
!   maxexp = Smallest power of beta that overflows

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XSMALL = Positive argument such that 1.0 - X = 1.0 to
!            machine precision for all ABS(X) .LE. XSMALL.
!   XINF =   Largest positive machine number; approximately beta**maxexp
!   XMAX =   Largest argument acceptable to BESI1;  Solution to equation:
!               EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp


!     Approximate values for some important machines are:

!                          beta       maxexp       XSMALL

! CRAY-1        (S.P.)       2         8191       3.55E-15
! Cyber 180/855
!   under NOS   (S.P.)       2         1070       3.55E-15
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)       2          128       2.98E-8
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)       2         1024       5.55D-17
! IBM 3033      (D.P.)      16           63       6.95D-18
! VAX           (S.P.)       2          127       2.98E-8
! VAX D-Format  (D.P.)       2          127       6.95D-18
! VAX G-Format  (D.P.)       2         1023       5.55D-17


!                               XINF          XMAX

! CRAY-1        (S.P.)       5.45E+2465     5682.810
! Cyber 180/855
!   under NOS   (S.P.)       1.26E+322       745.894
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)       3.40E+38         91.906
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)       1.79D+308       713.987
! IBM 3033      (D.P.)       7.23D+75        178.185
! VAX           (S.P.)       1.70D+38         91.209
! VAX D-Format  (D.P.)       1.70D+38         91.209
! VAX G-Format  (D.P.)       8.98D+307       713.293

!*******************************************************************

! Error returns

!  The program returns the value XINF for ABS(ARG) .GT. XMAX.


! Intrinsic functions required are:

!     ABS, SQRT, EXP


!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL  60439

!  Latest modification: March 13, 1992

!--------------------------------------------------------------------

REAL (8), INTENT(IN)   :: arg
REAL (8), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: jint

INTEGER    :: j
REAL (8)  :: a, b, sump, sumq, x, xx
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: one = 1.0_dp, one5 = 15.0_dp, forty = 40.0_dp,  &
                         exp40 = 2.353852668370199854D17, two25 = 225.0_dp, &
                         rec15 = 6.6666666666666666666D-2, half = 0.5_dp,  &
                         zero = 0.0_dp
!--------------------------------------------------------------------
!  Machine-dependent constants
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: XSMALL = 5.55D-17, XINF = 1.79D308, XMAX = 713.986D0
!--------------------------------------------------------------------
!  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: P(15) =  &
           (/ -1.9705291802535139930D-19, -6.5245515583151902910D-16,  &
              -1.1928788903603238754D-12, -1.4831904935994647675D-09,  &
              -1.3466829827635152875D-06, -9.1746443287817501309D-04,  &
              -4.7207090827310162436D-01, -1.8225946631657315931D+02,  &
              -5.1894091982308017540D+04, -1.0588550724769347106D+07,  &
              -1.4828267606612366099D+09, -1.3357437682275493024D+11,  &
              -6.9876779648010090070D+12, -1.7732037840791591320D+14,  &
              -1.4577180278143463643D+15 /)
REAL (8), PARAMETER  ::  Q(5) =   &
           (/ -4.0076864679904189921D+03, 7.4810580356655069138D+06,  &
              -8.0059518998619764991D+09, 4.8544714258273622913D+12,  &
              -1.3218168307321442305D+15 /)
!--------------------------------------------------------------------
!  Coefficients for 15.0 .LE. ABS(ARG)
!--------------------------------------------------------------------
REAL (8), PARAMETER  :: PP(8) =  &
           (/ -6.0437159056137600000D-02, 4.5748122901933459000D-01,  &
              -4.2843766903304806403D-01, 9.7356000150886612134D-02,  &
              -3.2457723974465568321D-03, -3.6395264712121795296D-04, &
               1.6258661867440836395D-05, -3.6347578404608223492D-07 /)
REAL (8), PARAMETER  :: QQ(6) =  &
           (/ -3.8806586721556593450D+00, 3.2593714889036996297D+00,  &
              -8.5017476463217924408D-01, 7.4212010813186530069D-02,  &
              -2.2835624489492512649D-03, 3.7510433111922824643D-05 /)
REAL (8), PARAMETER  :: PBAR = 3.98437500D-01
!--------------------------------------------------------------------

x = ABS(arg)
IF (x < xsmall) THEN
!--------------------------------------------------------------------
!  Return for ABS(ARG) .LT. XSMALL
!--------------------------------------------------------------------
  result = half * x
ELSE IF (x < one5) THEN
!--------------------------------------------------------------------
!  XSMALL .LE. ABS(ARG) .LT. 15.0
!--------------------------------------------------------------------
  xx = x * x
  sump = p(1)
  DO  j = 2, 15
    sump = sump * xx + p(j)
  END DO
  xx = xx - two25
  sumq = ((((xx+q(1))*xx+q(2))*xx+q(3))*xx+q(4)) * xx + q(5)
  result = (sump/sumq) * x
  IF (jint == 2) result = result * EXP(-x)
ELSE IF ((jint == 1) .AND. (x > xmax)) THEN
  result = xinf
ELSE
!--------------------------------------------------------------------
!  15.0 .LE. ABS(ARG)
!--------------------------------------------------------------------
  xx = one / x - rec15
  sump = ((((((pp(1)*xx+pp(2))*xx+pp(3))*xx+pp(4))*xx+pp(5))*xx+  &
      pp(6))*xx+pp(7)) * xx + pp(8)
  sumq = (((((xx+qq(1))*xx+qq(2))*xx+qq(3))*xx+qq(4))*xx+qq(5)) * xx + qq(6)
  result = sump / sumq
  IF (jint /= 1) THEN
    result = (result+pbar) / SQRT(x)
  ELSE
!--------------------------------------------------------------------
!  Calculation reformulated to avoid premature overflow
!--------------------------------------------------------------------
    IF (x > xmax-one5) THEN
      a = EXP(x-forty)
      b = exp40
    ELSE
      a = EXP(x)
      b = one
    END IF
    result = ((result*a+pbar*a)/SQRT(x)) * b
!--------------------------------------------------------------------
!  Error return for ABS(ARG) .GT. XMAX
!--------------------------------------------------------------------
  END IF
END IF
IF (arg < zero) result = -result
RETURN
!----------- Last line of CALCI1 -----------
END SUBROUTINE calci1

!------------------------------------------------------------------------------

FUNCTION BESI1(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This long precision subprogram computes approximate values for
!   modified Bessel functions of the first kind of order one for
!   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI1).

!--------------------------------------------------------------------

REAL (8), INTENT(IN)  :: x
REAL (8)              :: fn_val

! Local variable

INTEGER  :: jint
!--------------------------------------------------------------------
jint = 1
CALL calci1(x, fn_val, jint)
RETURN
!---------- Last line of BESI1 ----------
END FUNCTION BESI1

!------------------------------------------------------------------------------
!  BEGIN THE NEW FUNCTIONS FOR CALCULATING HEAT CAPACITY
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This function calculates the heat capacity for hindered rotors.
! It is called by the "dqed" solver in the section for "diff"
!------------------------------------------------------------------------------
FUNCTION cv_hind_rot(V,T,nu)  RESULT(fn_val)

  IMPLICIT NONE
  REAL(8) :: fn_val                  ! Unitless
  REAL(8), INTENT(IN) :: T           ! Kelvin
  REAL(8), INTENT(IN) :: V           ! cm^-1
  REAL(8), INTENT(IN) :: nu          ! cm^-1
  REAL(8) :: kB = 0.695039         ! cm^-1 / K
  REAL(8) :: u, a
  REAL(8) :: BB
  REAL(8) :: expo, in_expo

  u = V/(2. * kB * T)
  a = nu/V
  expo = exp(-2*u*a)
  in_expo = (1.0 / (1. - expo))
  BB =  BESI1(u) / BESI0(u)

  fn_val = - 0.5 + u * ( u - BB - u * BB*BB )                                  &
                 + 4 * u*u * a*a * expo  * in_expo * in_expo 

END FUNCTION cv_hind_rot

!------------------------------------------------------------------------------
!  This function is the partial derivative of the hindered-rotor heat capacity
!  With respect to the barrier height V.
! It is called by the "dqed" solver in the section for the jacobian.
!------------------------------------------------------------------------------
FUNCTION d_cv_hind_rot(V,T,nu)  RESULT(fn_val)
  
  IMPLICIT NONE
  REAL(8) :: fn_val                  ! Unitless
  REAL(8), INTENT(IN) :: T           ! Kelvin
  REAL(8), INTENT(IN) :: V           ! cm^-1
  REAL(8), INTENT(IN) :: nu          ! cm^-1
  REAL(8) :: kB = 0.695039         ! cm^-1 / K
  REAL(8) :: u, a
  REAL(8) :: BB
  REAL(8) :: expo, in_expo

  u = V/(2. * kB * T)
  a = nu/V
  expo = exp(-2.*u*a)
  in_expo = 1.0/(expo-1.)
  BB =  BESI1(u) / BESI0(u)

  fn_val = (u - BB - u * BB*BB) / (2. * kB * T)                 &
         + u * ( 1./ (2. * kB * T) - (1. - BB/u)/ (2. * kB * T) &
         - u * (BB - BB*BB/u) / ( kB * T)               & 
         + u * BB*BB*BB / (kB * T)  )

END FUNCTION d_cv_hind_rot


!------------------------------------------------------------------------------
!  This function is the partial derivative of the hindered-rotor heat capacity
!  With respect to the frequency nu.
! It is called by the "dqed" solver in the section for the jacobian.
!------------------------------------------------------------------------------
FUNCTION d_cv_hind_rot_nu(V,T,nu)  RESULT(fn_val)
  
  IMPLICIT NONE
  REAL(8) :: fn_val                  ! Unitless
  REAL(8), INTENT(IN) :: T           ! Kelvin
  REAL(8), INTENT(IN) :: V           ! cm^-1
  REAL(8), INTENT(IN) :: nu          ! cm^-1
  REAL(8) :: kB = 0.695039         ! cm^-1 / K
  REAL(8) :: u, a
  REAL(8) :: BB
  REAL(8) :: expo, in_expo

  a = nu/(kB * T)
  expo = exp(-a)
  in_expo = 1.0/(expo-1.)
 

  fn_val = a /(kB * T) * expo * in_expo*in_expo * (2 - a + 2.*a * expo * in_expo )

END FUNCTION d_cv_hind_rot_nu

!------------------------------------------------------------------------------
!  This function calculates the heat capacity for harmonic oscillators.
! It is called by the "dqed" solver in the section for "diff"
!------------------------------------------------------------------------------
Function cv_harm_osc(nu, T)

  IMPLICIT NONE
  REAL(8) :: cv_harm_osc
  REAL(8), INTENT(IN) :: nu        ! cm^-1
  REAL(8), INTENT(IN) :: T         ! K
  REAL(8) :: kB = 0.695039         ! cm^-1 / K
  REAL(8) :: rat, expo, in_expo


  rat = nu / (kB * T)
  expo = exp(rat)
  in_expo = (1.0 / (expo - 1.))

  cv_harm_osc = expo * (rat * in_expo)*(rat * in_expo)

END FUNCTION cv_harm_osc

!------------------------------------------------------------------------------
!  This function is the partial derivative of the harmonic oscillator heat
!  Capacity with respect to the frequency nu.
! It is called by the "dqed" solver in the section for the jacobian.
!------------------------------------------------------------------------------
FUNCTION d_Cv_harm_osc(nu,T)

  IMPLICIT NONE
  REAL(8) :: d_Cv_harm_osc
  REAL(8), INTENT(IN) :: nu        ! cm^-1
  REAL(8), INTENT(IN) :: T         ! K
  REAL(8) :: kB = 0.695039         ! cm^-1 / K
  REAL(8) :: rat, expo, in_expo

  rat = nu / (kB * T) 
  expo = exp( rat )
  in_expo = (1.0 / (expo - 1.))


  d_Cv_harm_osc = expo * rat * in_expo * in_expo /  (kB * T)          &
                * (2. + rat - 2. * rat * expo * in_expo)

END FUNCTION d_Cv_harm_osc

END MODULE heat_capacity_functions 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------



MODULE frequencies

CONTAINS
!------------------------------------------------------------------------------

SUBROUTINE calc_predicted_freq(bond_info, Bond_degeneracy, Total_predicted_freq)
!------------------------------------------------------------------------------
! This code contains the characteristic or typical frequencies for both
! bond stretching and for bond bending.
!
! The data are combined with the number of a particular type of bond
! (e.g. the number of hydrogen bonded with a carbon that shares a double bond
! with a second carbon), which is taken from the connectivity table.
!
! The connectivity table and list of characteristic frequencies are used to 
! calculate as many of the characteristic bending and stretching frequencies as
! possible.  These characteristic frequencies are then used to calculate the 
! contribution of these bonds to the overall heat capacity using the usual 
! Stat-mech formulas.
!
! ALL FREQUENCIES SHOULD BE GIVEN IN UNITS OF cm^-1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  BEGIN THE LIST OF ALL THE CONSTANTS
!------------------------------------------------------------------------------

  IMPLICIT NONE

! The array which contains all the bond info
  INTEGER, DIMENSION(:), INTENT(IN) :: bond_info

! Number of atoms, number of hindered internal rotors, linearity
  INTEGER :: N_atom, N_rot, linearity

! All the integers corresponding to the number of particular bonds
  INTEGER :: RsCH3, RdCH2, CtCH, RSCH2sR, CdCHsR, Aldehyde,Cumulene
  INTEGER :: Ketene, CtCsR, RsCHsR2, CdCsR2, Ketone, RsCsR3
  INTEGER :: RsCH2r, RdCHr, RsCHrsR, CdCrsR, OdCrsR, RsCrsR2
  INTEGER :: Alcohol, Ether, ROOH, ROOR, Peroxy
  INTEGER :: Rings

! All the arrays corresponding to the characteristic frequencies for each bond
  REAL(8), ALLOCATABLE, DIMENSION(:) :: RsCH3_pred_freq, RdCH2_pred_freq 
  REAL(8), ALLOCATABLE, DIMENSION(:) :: CtCH_pred_freq, RsCH2sR_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: CdCHsR_pred_freq, Aldehyde_pred_freq 
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Cumulene_pred_freq, Ketene_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: CtCsR_pred_freq, RsCHsR2_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: CdCsR2_pred_freq, Ketone_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: RsCsR3_pred_freq, RsCH2r_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: RdCHr_pred_freq, RsCHrsR_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: CdCrsR_pred_freq, OdCrsR_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: RsCrsR2_pred_freq, Alcohol_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Ether_pred_freq, ROOH_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: ROOR_pred_freq, peroxy_pred_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Rings_pred_freq  


! The total number of bonds specified by the connectivity table
  INTEGER, INTENT(IN) :: Bond_degeneracy

! The array corresponding to all the characteristic frequencies
  REAL(8), INTENT(OUT), DIMENSION(Bond_degeneracy) :: Total_predicted_freq

!------------------------------------------------------------------------------
!  END THE LIST OF ALL THE CONSTANTS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! BEGIN THE LIST OF TYPICAL FREQUENCIES 
!------------------------------------------------------------------------------


  REAL(8), DIMENSION(10) :: RsCH3_typ_freq    = (/ 2750., 2850.0,  &
                                                1350.0, 1500.0,   &
                                                 700.0, 800.0,    &
                                                1000.0, 1100.0,   &
                                                1350.0, 1400.0 /)

  REAL(8), DIMENSION(8)  :: RdCH2_typ_freq    = (/ 2950.0, 3100.0,  &
                                                   1330.0, 1430.0,  &
                                                    900.0, 1050.0,  &
                                                   1000.0, 1050.0  /)

  REAL(8), DIMENSION(4)  :: CtCH_typ_freq     = (/  750.0, 770.0,  &
                                                  3350.0, 3450.0   /)

  REAL(8), DIMENSION(12) :: RsCH2sR_typ_freq  = (/ 2750.0, 2850.0,  &
                                                   1425.0, 1450.0,  &
                                                   1225.0, 1275.0,  &
                                                   1270.0, 1340.0,  &
                                                    700.0,  800.0,  &
                                                    300.0,  400.0   /)

  REAL(8), DIMENSION(10) :: CdCHsR_typ_freq   = (/ 2995.0, 3025.0,  &
                                                    975.0, 1000.0,  &
                                                   1300.0, 1375.0,  &
                                                    400.0,  500.0,  &
                                                   1630.0, 1680.0   /)

  REAL(8), DIMENSION(10) :: Aldehyde_typ_freq = (/ 2695.0, 2870.0,  &
                                                    700.0,  800.0,  &
                                                   1380.0, 1410.0,  &
                                                    450.0,  500.0,  &
                                                   1750.0, 1800.0   /)
       
  REAL(8), DIMENSION(4)  :: Cumulene_typ_freq = (/  540.0,  610.0,  &
                                                   1970.0, 2140.0   /)

  REAL(8), DIMENSION(6)  :: Ketene_typ_freq   = (/ 2110.0, 2130.0,  &
                                                    495.0,  530.0,  &
                                                    650.0,  925.0   /)

  REAL(8), DIMENSION(4)  :: CtCsR_typ_freq    = (/ 2100.0, 2250.0,  &
                                                    500.0,  550.0   /)
       
  REAL(8), DIMENSION(8) :: RsCHsR2_typ_freq  = (/  1380.0, 1390.0,  &
                                                    370.0,  380.0,  &
                                                   2800.0, 3000.0,  &
                                                    430.0,  440.0   /)

  REAL(8), DIMENSION(8) :: CdCsR2_typ_freq   = (/   325.0,  375.0,  &
                                                    415.0,  465.0,  &
                                                    420.0,  450.0,  &
                                                   1700.0, 1750.0   /)

  REAL(8), DIMENSION(8) :: Ketone_typ_freq   = (/   365.0,  385.0,  &
                                                    505.0,  600.0,  &
                                                    445.0,  480.0,  &
                                                   1700.0, 1720.0   /)

  REAL(8), DIMENSION(6) :: RsCsR3_typ_freq   = (/   350.0,  400.0,  &
                                                   1190.0, 1240.0,  &
                                                    400.0,  500.0   /)

  REAL(8), DIMENSION(8)  :: RsCH2r_typ_freq   = (/ 3000.0, 3100.0,  &
                                                    415.0,  465.0,  &
                                                    780.0,  850.0,  &
                                                   1435.0, 1475.0   /)

  REAL(8), DIMENSION(6)  :: RdCHr_typ_freq    = (/ 3115.0, 3125.0,  &
                                                    620.0,  680.0,  &
                                                    785.0,  800.0   /)

  REAL(8), DIMENSION(8) :: RsCHrsR_typ_freq  = (/  3000.0, 3050.0,  &
                                                    390.0,  425.0,  &
                                                   1340.0, 1360.0,  &
                                                    335.0,  370.0   /)

  REAL(8), DIMENSION(4)  :: CdCrsR_typ_freq   = (/ 1670.0, 1700.0,  &
                                                    300.0,  440.0   /)

  REAL(8), DIMENSION(4)  :: OdCrsR_typ_freq   = (/ 1850.0, 1860.0,  &
                                                    440.0,  470.0   /)

  REAL(8), DIMENSION(4)  :: RsCrsR2_typ_freq  = (/  360.0,  370.0,  &
                                                    300.0,  400.0   /)

  REAL(8), DIMENSION(4)  :: Alcohol_typ_freq  = (/ 3580.0, 3650.0,  &
                                                   1210.0, 1345.0   /)

  REAL(8), DIMENSION(2)  :: Ether_typ_freq    = (/  350.0,  500.0   /)

  REAL(8), DIMENSION(8) :: ROOH_typ_freq     = (/ 3580.0, 3650.0,  &
                                                   1300.0, 1320.0,  &
                                                    350.0,  425.0,  &
                                                    825.0,  875.0   /)

  REAL(8), DIMENSION(4)  :: ROOR_typ_freq     = (/  350.0,  500.0,  &
                                                    795.0,  815.0   /)

  REAL(8), DIMENSION(4)  :: peroxy_typ_freq   = (/  470.0,  515.0,  &
                                                   1100.0, 1170.0   /)

  REAL(8), DIMENSION(4) :: Rings_typ_freq    = (/ 2750.0, 3150.0,   &
                                                   900.0, 1100.0    /)



!------------------------------------------------------------------------------
! END THE LIST OF TYPICAL FREQUENCIES 
!------------------------------------------------------------------------------

! Unpack molecular information

! non radicals
! terminal groups
  RsCH3     = bond_info(1)
  RdCH2     = bond_info(2)
  CtCH      = bond_info(3)

! non-terminal groups (2 heavy atoms)
  RsCH2sR   = bond_info(4)
  CdCHsR    = bond_info(5)
  Aldehyde  = bond_info(6)
  Cumulene  = bond_info(7)
  Ketene    = bond_info(8)
  CtCsR     = bond_info(9)
  
  ! non-terminal groups (3 heavy atoms)
  RsCHsR2   = bond_info(10)
  CdCsR2    = bond_info(11)
  Ketone    = bond_info(12)

  ! non-terminal groups (4 heavy atoms)
  RsCsR3    = bond_info(13)

! Single radicals
! terminal groups
  RsCH2r    = bond_info(14)
  RdCHr     = bond_info(15)

! non-terminal groups
  RsCHrsR   = bond_info(16)
  CdCrsR    = bond_info(17)
  OdCrsR    = bond_info(18)
  RsCrsR2   = bond_info(19)

! oxygen-specific groups
  Alcohol   = bond_info(20)
  Ether     = bond_info(21)
  ROOH      = bond_info(22)
  ROOR      = bond_info(23)
  Peroxy    = bond_info(24)

! ring-specific groups
  Rings     = bond_info(25)

! Useful for debugging
 WRITE(*,*) ' RsCH3 = ', RsCH3
 WRITE(*,*) ' RdCH2 = ', RdCH2
 WRITE(*,*) ' CtCH = ', CtCH
 WRITE(*,*) ' RSCH2sR = ', RSCH2sR
 WRITE(*,*) ' CdCHsR = ', CdCHsR
 WRITE(*,*) ' Aldehyde = ', Aldehyde
 WRITE(*,*) ' Cumulene = ', Cumulene
 WRITE(*,*) ' Ketene = ', Ketene
 WRITE(*,*) ' CtCsR = ', CtCsR
 WRITE(*,*) ' RsCHsR2 = ', RsCHsR2
 WRITE(*,*) ' CdCsR2 = ', CdCsR2
 WRITE(*,*) ' Ketone = ', Ketone
 WRITE(*,*) ' RsCsR3 = ', RsCsR3
 WRITE(*,*) ' RsCH2r = ', RsCH2r
 WRITE(*,*) ' RdCHr = ', RdCHr
 WRITE(*,*) ' RsCHrsR = ', RsCHrsR
 WRITE(*,*) ' CdCrsR = ', CdCrsR
 WRITE(*,*) ' OdCrsR = ', OdCrsR
 WRITE(*,*) ' RsCrsR2 = ', RsCrsR2
 WRITE(*,*) ' Alcohol = ', Alcohol
 WRITE(*,*) ' Ether = ', Ether
 WRITE(*,*) ' ROOH = ', ROOH
 WRITE(*,*) ' ROOR = ', ROOR
 WRITE(*,*) ' Peroxy = ', Peroxy
 WRITE(*,*) ' Rings = ', Rings


!------------------------------------------------------------------------------
!  BEGIN ALLOCATING ALL THE VECTORS FOR PREDICTED FREQUENCIES
!  The length of each vector is the number of instances of a particular type
!  times the number frequencies that correspond to that type
!------------------------------------------------------------------------------

  ALLOCATE(  RsCH3_pred_freq(    RsCH3    * 8 ) )
  ALLOCATE(  RdCH2_pred_freq(    RdCH2    * 5 ) )
  ALLOCATE(  CtCH_pred_freq(     CtCH     * 3 ) )
  ALLOCATE(  RsCH2sR_pred_freq(  RsCH2sR  * 7 ) )
  ALLOCATE(  CdCHsR_pred_freq(   CdCHsR   * 5 ) )
  ALLOCATE(  Aldehyde_pred_freq( Aldehyde * 5 ) )
  ALLOCATE(  Cumulene_pred_freq( Cumulene * 3 ) )

  ALLOCATE(  Ketene_pred_freq(   Ketene   * 3 ) )
  ALLOCATE(  CtCsR_pred_freq(    CtCsR    * 2 ) )
  ALLOCATE(  RsCHsR2_pred_freq(  RsCHsR2  * 6 ) )
  ALLOCATE(  CdCsR2_pred_freq(   CdCsR2   * 4 ) )
  ALLOCATE(  Ketone_pred_freq(   Ketone   * 4 ) )
  ALLOCATE(  RsCsR3_pred_freq(   RsCsR3   * 5 ) )
  ALLOCATE(  RsCH2r_pred_freq(   RsCH2r   * 5 ) )
  ALLOCATE(  RdCHr_pred_freq(    RdCHr    * 3 ) )

  ALLOCATE(  RsCHrsR_pred_freq(  RsCHrsR  * 4 ) )
  ALLOCATE(  CdCrsR_pred_freq(   CdCrsR   * 2 ) )
  ALLOCATE(  OdCrsR_pred_freq(   OdCrsR   * 2 ) )
  ALLOCATE(  RsCrsR2_pred_freq(  RsCrsR2  * 3 ) )
  ALLOCATE(  Alcohol_pred_freq(  Alcohol  * 2 ) )
  ALLOCATE(  Ether_pred_freq(    Ether    * 1 ) )
  ALLOCATE(  ROOH_pred_freq(     ROOH     * 4 ) )
  ALLOCATE(  ROOR_pred_freq(     ROOR     * 2 ) )
  ALLOCATE(  peroxy_pred_freq(   peroxy   * 2 ) )

  ALLOCATE(  Rings_pred_freq(    Rings    * 2 ) )

!------------------------------------------------------------------------------
!  END ALLOCATING ALL THE VECTORS FOR PREDICTED FREQUENCIES
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  BEGIN calculating the characteristic frequencies
!------------------------------------------------------------------------------

  CALL  Bond_avg_8x3_2_1(RsCH3,  RsCH3_typ_freq, RsCH3_pred_freq)
  CALL  Bond_avg_5x2_1(RdCH2,  RdCH2_typ_freq, RdCH2_pred_freq)
  CALL  Bond_avg_3x2_1(CtCH,  CtCH_typ_freq, CtCH_pred_freq)
  CALL  Bond_avg_7x2_1(RsCH2sR,  RsCH2sR_typ_freq, RsCH2sR_pred_freq)
  CALL  Bond_avg_5x1(CdCHsR,  CdCHsR_typ_freq, CdCHsR_pred_freq)
  CALL  Bond_avg_5x1(Aldehyde,  Aldehyde_typ_freq, Aldehyde_pred_freq)
  CALL  Bond_avg_3x2_1(Cumulene,  Cumulene_typ_freq, Cumulene_pred_freq)
  CALL  Bond_avg_3x1(Ketene, Ketene_typ_freq, Ketene_pred_freq)
  CALL  Bond_avg_2x1(CtCsR, CtCsR_typ_freq, CtCsR_pred_freq)
  CALL  Bond_avg_6x2_2_1(RsCHsR2, RsCHsR2_typ_freq, RsCHsR2_pred_freq)
  CALL  Bond_avg_4x1(CdCsR2,  CdCsR2_typ_freq, CdCsR2_pred_freq)
  CALL  Bond_avg_4x1(Ketone,  Ketone_typ_freq, Ketone_pred_freq)
  CALL  Bond_avg_5x2_2_1(RsCsR3, RsCsR3_typ_freq, RsCsR3_pred_freq)
  CALL  Bond_avg_5x2_1(RsCH2r,  RsCH2r_typ_freq, RsCH2r_pred_freq)
  CALL  Bond_avg_3x1(RdCHr,  RdCHr_typ_freq, RdCHr_pred_freq)
  CALL  Bond_avg_4x1(RsCHrsR,  RsCHrsR_typ_freq, RsCHrsR_pred_freq)
  CALL  Bond_avg_2x1(CdCrsR,  CdCrsR_typ_freq, CdCrsR_pred_freq)
  CALL  Bond_avg_2x1(OdCrsR,  OdCrsR_typ_freq, OdCrsR_pred_freq)
  CALL  Bond_avg_3x2_1(RsCrsR2, RsCrsR2_typ_freq, RsCrsR2_pred_freq)
  CALL  Bond_avg_2x1(Alcohol,  Alcohol_typ_freq, Alcohol_pred_freq)
  CALL  Bond_avg_1x1(Ether,  Ether_typ_freq, Ether_pred_freq)
  CALL  Bond_avg_4x1(ROOH,  ROOH_typ_freq, ROOH_pred_freq)
  CALL  Bond_avg_2x1(ROOR,  ROOR_typ_freq, ROOR_pred_freq)
  CALL  Bond_avg_2x1(peroxy,  peroxy_typ_freq, peroxy_pred_freq)
  CALL  Bond_avg_2x1(Rings,  Rings_typ_freq, Rings_pred_freq)

!------------------------------------------------------------------------------
!  END calculating the bond degeneracy, sump cp bonds, and char. freqs.
!------------------------------------------------------------------------------

! --- CALCULATE THE TOTAL CHARACTERISTIC FREQUENCIES
! append all the individual predicted frequencies to one long vector


  Total_predicted_freq = (/  RsCH3_pred_freq, RdCH2_pred_freq,     &
                             CtCH_pred_freq, RsCH2sR_pred_freq,    &
                             CdCHsR_pred_freq, Aldehyde_pred_freq, &
                             Cumulene_pred_freq, Ketene_pred_freq, & 
                             CtCsR_pred_freq, RsCHsR2_pred_freq,   &
                             CdCsR2_pred_freq, Ketone_pred_freq,   &
                             RsCsR3_pred_freq, RsCH2r_pred_freq,   &
                             RdCHr_pred_freq, RsCHrsR_pred_freq,   &
                             CdCrsR_pred_freq, OdCrsR_pred_freq,   &
                             RsCrsR2_pred_freq, Alcohol_pred_freq, &
                             Ether_pred_freq, ROOH_pred_freq,      &
                             ROOR_pred_freq, peroxy_pred_freq,     &
                             Rings_pred_freq                       /)

 !------------------------------------------------------------------------------
!   Deallocate all the characteristic frequency arrays
!------------------------------------------------------------------------------

  DEALLOCATE(  RsCH3_pred_freq )
  DEALLOCATE(  RdCH2_pred_freq )
  DEALLOCATE(  CtCH_pred_freq )
  DEALLOCATE(  RsCH2sR_pred_freq )
  DEALLOCATE(  CdCHsR_pred_freq )
  DEALLOCATE(  Aldehyde_pred_freq )
  DEALLOCATE(  Cumulene_pred_freq )
  DEALLOCATE(  Ketene_pred_freq  )
  DEALLOCATE(  CtCsR_pred_freq )
  DEALLOCATE(  RsCHsR2_pred_freq )
  DEALLOCATE(  CdCsR2_pred_freq )
  DEALLOCATE(  Ketone_pred_freq )
  DEALLOCATE(  RsCsR3_pred_freq )
  DEALLOCATE(  RsCH2r_pred_freq )
  DEALLOCATE(  RdCHr_pred_freq )
  DEALLOCATE(  RsCHrsR_pred_freq )
  DEALLOCATE(  CdCrsR_pred_freq )
  DEALLOCATE(  OdCrsR_pred_freq )
  DEALLOCATE(  RsCrsR2_pred_freq )
  DEALLOCATE(  Alcohol_pred_freq )
  DEALLOCATE(  Ether_pred_freq )
  DEALLOCATE(  ROOH_pred_freq )
  DEALLOCATE(  ROOR_pred_freq )
  DEALLOCATE(  peroxy_pred_freq )
  DEALLOCATE(  Rings_pred_freq ) 

END SUBROUTINE calc_predicted_freq
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! The following subroutines are used to interpolate the correct number of 
! bonds.  They are all pretty much the same in nature.
! For each vibrational mode, the vector typical_frequencies will have an
! upper and lower bound.  If the there is one instance of that frequency, then
! the code will return the average between the upper and lower bounds.
! If there are two instances of the frequency, it will return the two bounds.
! If there are more than two, then it will call a subroutine, linspace, 
! (like the function linspace in MATLAB), which returns a vector from low to 
! high in evenly spaced increments.
!
!  The naming convention is as follows:
!  The first number (to the left of the x) is the number of frequencies the
!  subroutine will return per bond type (e.g. Bond_avg_8x will return 8 times
!  Number_of_bonds frequencies.
!  
!  The numbers to the right of the x indicate whether a particular frequency
!  range will return more than one frequency within that range.
!  For example:  x3_2_1  will return 3 frequencies in the first range,
!  2 frequencies in the second range, and all the remaining ranges will return
!  only 1 frequency per range.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_8x3_2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(8*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          3, predicted_frequencies(1:3 ) )

     predicted_frequencies(4) = typical_frequencies(3)
     predicted_frequencies(5) = typical_frequencies(4)

     predicted_frequencies(6) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0
     predicted_frequencies(7) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0     
     predicted_frequencies(8) = ( typical_frequencies(9) + &
                                  typical_frequencies(10) ) / 2.0

  ELSE
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          3*Number_of_bonds, predicted_frequencies( 1:3*Number_of_bonds ) )
     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          2*Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(5*number_of_bonds) ) )

     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (5*Number_of_bonds+1):(6*number_of_bonds) ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (6*Number_of_bonds+1):(7*number_of_bonds) ) )
     call linspace(typical_frequencies(9),typical_frequencies(10), &
          Number_of_bonds, predicted_frequencies( (7*Number_of_bonds+1):(8*number_of_bonds) ) )

     
  ENDIF

END SUBROUTINE Bond_avg_8x3_2_1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_5x2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(5*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN
     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)

     predicted_frequencies(3) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0
     predicted_frequencies(4) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0     
     predicted_frequencies(5) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0

  ELSE
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(4*number_of_bonds) ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (4*Number_of_bonds+1):(5*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_5x2_1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_3x2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(3*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN
     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)

     predicted_frequencies(3) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0
   
  ELSE
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_3x2_1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_7x2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(7*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN
     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)

     predicted_frequencies(3) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0
     predicted_frequencies(4) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0     
     predicted_frequencies(5) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0
     predicted_frequencies(6) = ( typical_frequencies(9) + &
                                  typical_frequencies(10) ) / 2.0
     predicted_frequencies(7) = ( typical_frequencies(11) + &
                                  typical_frequencies(12) ) / 2.0     
     
  ELSE
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(4*number_of_bonds) ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (4*Number_of_bonds+1):(5*number_of_bonds) ) )
     call linspace(typical_frequencies(9),typical_frequencies(10), &
          Number_of_bonds, predicted_frequencies( (5*Number_of_bonds+1):(6*number_of_bonds) ) )
     call linspace(typical_frequencies(11),typical_frequencies(12), &
          Number_of_bonds, predicted_frequencies( (6*Number_of_bonds+1):(7*number_of_bonds) ) )

  ENDIF

END SUBROUTINE Bond_avg_7x2_1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_5x1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(5*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = ( typical_frequencies(1) + &
                                  typical_frequencies(2) ) / 2.0
     predicted_frequencies(2) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0     
     predicted_frequencies(3) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0
     predicted_frequencies(4) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0
     predicted_frequencies(5) = ( typical_frequencies(9) + &
                                  typical_frequencies(10) ) / 2.0     

    
  ELSE

     call linspace(typical_frequencies(1),typical_frequencies(2), &
          Number_of_bonds, predicted_frequencies( 1:number_of_bonds ) )     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (Number_of_bonds+1):(2*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(4*number_of_bonds) ) )
     call linspace(typical_frequencies(9),typical_frequencies(10), &
          Number_of_bonds, predicted_frequencies( (4*Number_of_bonds+1):(5*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_5x1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_3x1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(3*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = ( typical_frequencies(1) + &
                                  typical_frequencies(2) ) / 2.0
     predicted_frequencies(2) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0     
     predicted_frequencies(3) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0

  ELSE

     call linspace(typical_frequencies(1),typical_frequencies(2), &
          Number_of_bonds, predicted_frequencies( 1:number_of_bonds ) )     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (  Number_of_bonds+1):(2*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_3x1
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

SUBROUTINE Bond_avg_6x2_2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(6*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
!     WRITE (*,*) 'N = 0'
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
 ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)
     predicted_frequencies(3) = typical_frequencies(3)
     predicted_frequencies(4) = typical_frequencies(4)

     predicted_frequencies(5) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0     
     predicted_frequencies(6) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0
   

  ELSE   
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          2*Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):4*Number_of_bonds ) )
    
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (4*Number_of_bonds+1):5*Number_of_bonds ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (5*Number_of_bonds+1):(6*number_of_bonds) ) )

ENDIF

END SUBROUTINE Bond_avg_6x2_2_1
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_4x1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(4*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = ( typical_frequencies(1) + &
                                  typical_frequencies(2) ) / 2.0
     predicted_frequencies(2) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0     
     predicted_frequencies(3) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0
     predicted_frequencies(4) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0
     

    
  ELSE

     call linspace(typical_frequencies(1),typical_frequencies(2), &
          Number_of_bonds, predicted_frequencies( 1:number_of_bonds ) )     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (Number_of_bonds+1):(2*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(4*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_4x1
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

SUBROUTINE Bond_avg_5x2_2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(5*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
!     WRITE (*,*) 'N = 0'
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
 ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)
     predicted_frequencies(3) = typical_frequencies(3)
     predicted_frequencies(4) = typical_frequencies(4)

     predicted_frequencies(5) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0     
   

  ELSE   
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          2*Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):4*Number_of_bonds ) )
    
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (4*Number_of_bonds+1):5*Number_of_bonds ) )

ENDIF

END SUBROUTINE Bond_avg_5x2_2_1
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_6x2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(6*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN
     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)

     predicted_frequencies(3) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0
     predicted_frequencies(4) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0     
     predicted_frequencies(5) = ( typical_frequencies(7) + &
                                  typical_frequencies(8) ) / 2.0
     predicted_frequencies(6) = ( typical_frequencies(9) + &
                                  typical_frequencies(10) ) / 2.0

  ELSE
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(4*number_of_bonds) ) )
     call linspace(typical_frequencies(7),typical_frequencies(8), &
          Number_of_bonds, predicted_frequencies( (4*Number_of_bonds+1):(5*number_of_bonds) ) )
     call linspace(typical_frequencies(9),typical_frequencies(10), &
          Number_of_bonds, predicted_frequencies( (5*Number_of_bonds+1):(6*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_6x2_1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_2x1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(2*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = ( typical_frequencies(1) + &
                                  typical_frequencies(2) ) / 2.0
     predicted_frequencies(2) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0     

  ELSE

     call linspace(typical_frequencies(1),typical_frequencies(2), &
          Number_of_bonds, predicted_frequencies( 1:number_of_bonds ) )     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (  Number_of_bonds+1):(2*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_2x1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_1x1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(2), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(1*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN

     predicted_frequencies(1) = ( typical_frequencies(1) + &
                                  typical_frequencies(2) ) / 2.0

  ELSE

     call linspace(typical_frequencies(1),typical_frequencies(2), &
          Number_of_bonds, predicted_frequencies( 1:number_of_bonds ) )     

  ENDIF

END SUBROUTINE Bond_avg_1x1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Bond_avg_4x2_1 (Number_of_bonds, & 
                              typical_frequencies, predicted_frequencies)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number_of_bonds
  REAL(8), DIMENSION(:), INTENT(IN) :: typical_frequencies
  REAL(8),DIMENSION(4*Number_of_bonds), INTENT(OUT) :: predicted_frequencies
  INTEGER                                      :: i,j, k

  IF (Number_of_bonds == 0 ) THEN
     predicted_frequencies = (/ 0.0, 0.0, 0.0, 0.0 /) 
  
  ELSEIF (Number_of_bonds == 1) THEN
     predicted_frequencies(1) = typical_frequencies(1)
     predicted_frequencies(2) = typical_frequencies(2)

     predicted_frequencies(3) = ( typical_frequencies(3) + &
                                  typical_frequencies(4) ) / 2.0
     predicted_frequencies(4) = ( typical_frequencies(5) + &
                                  typical_frequencies(6) ) / 2.0     

  ELSE
     call linspace(typical_frequencies(1),typical_frequencies(2), &
          2*Number_of_bonds, predicted_frequencies( 1:2*Number_of_bonds ) )
     
     call linspace(typical_frequencies(3),typical_frequencies(4), &
          Number_of_bonds, predicted_frequencies( (2*Number_of_bonds+1):(3*number_of_bonds) ) )
     call linspace(typical_frequencies(5),typical_frequencies(6), &
          Number_of_bonds, predicted_frequencies( (3*Number_of_bonds+1):(4*number_of_bonds) ) )
     
  ENDIF

END SUBROUTINE Bond_avg_4x2_1
!-------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE linspace(low,high,n,out)
  IMPLICIT NONE
  REAL(8), INTENT(in) :: low
  REAL(8), INTENT(in) :: high
  INTEGER, INTENT(in) :: n
  REAL(8), DIMENSION(n) :: out
  INTEGER :: i
  REAL(8) :: increment

  if (n==2) THEN
     out = (/low, high/)
  else
     out(1) = low
     increment = (high - low) /( n-1)
     DO i = 1, (n-1)
        out(i+1) = out(i) + increment
     END DO
  END IF
  
END SUBROUTINE linspace

!------------------------------------------------------------------------------
   
SUBROUTINE Cv_vib(frequency, CV_temps, Cv)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: frequency   ! cm^-1
  REAL(8), INTENT(IN) :: CV_temps ! K
  REAL(8) :: kB = 0.695039           ! cm^-1 / K
  REAL(8) :: theta
  REAL(8), INTENT(OUT) :: Cv

  theta = ( frequency / (kB* CV_temps) )
  Cv = exp(theta) * ( theta / (exp(theta) - 1) )**2

END SUBROUTINE Cv_vib

!------------------------------------------------------------------------------
END MODULE frequencies

MODULE Cases

CONTAINS


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Begin Solver section
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! This file contains all the cases concerning the number of unknown harmonic
! Oscillator frequencies and hindered internal rotors
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 0:  N_vib = 0, N_rot = 0
!------------------------------------------------------------------------------

SUBROUTINE Case_0

END SUBROUTINE Case_0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 1:  N_vib = 1, N_rot = 0
!------------------------------------------------------------------------------

SUBROUTINE Case_1(Total_char_freq, Total_harm_osc_freq, HR_params  )
  USE heat_capacity_functions
 
  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 1
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_1_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1


!  Define the initial guess for the solution
  x(1) = 800.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_1_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 1 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 = results(1)
  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1


END SUBROUTINE Case_1
!------------------------------------------------------------------------------
SUBROUTINE Case_1_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 1
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i))                                 &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 

    END DO

  END IF

  RETURN

END SUBROUTINE Case_1_hd
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
! CASE 2:  N_vib = 2, N_rot = 0
!------------------------------------------------------------------------------

SUBROUTINE Case_2(Total_char_freq, Total_harm_osc_freq, HR_params  )
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 2
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_2_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 300.0 !Harm. osc. frequency
  x(2) = 2000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_2_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 2 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '  
  results = x

  Theta_1 = results(1)
  Theta_2 = results(2)
  Total_harm_osc_freq( size(Total_char_freq ) +1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) +2 ) = Theta_2

END SUBROUTINE Case_2
!------------------------------------------------------------------------------
SUBROUTINE Case_2_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 2
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i)) + cv_harm_osc(x(2),t(i))                                &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i)) 

    END DO

  END IF

  RETURN

END SUBROUTINE Case_2_hd

!------------------------------------------------------------------------------
! CASE 3:  N_vib = 3, N_rot = 0
!------------------------------------------------------------------------------
SUBROUTINE Case_3(Total_char_freq, Total_harm_osc_freq, HR_params )
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 3
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_3_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the third harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 300.0 !Harm. osc. frequency
  x(2) = 1000.0 !Harm. osc. frequency
  x(3) = 2000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_3_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 3 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 = results(1)
  Theta_2 = results(2)
  Theta_3 = results(3)
  Total_harm_osc_freq( size(Total_char_freq ) +1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) +2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) +3 ) = Theta_3

END SUBROUTINE Case_3
!------------------------------------------------------------------------------
SUBROUTINE Case_3_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 3
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i)) &
               + cv_harm_osc(x(2),t(i)) &
               + cv_harm_osc(x(3),t(i)) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i)) 
      fj(mcon+i,3) = 2. * diff(i) * d_Cv_harm_osc(x(3),t(i)) 

    END DO

  END IF

  RETURN

END SUBROUTINE Case_3_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 4:  N_vib = 4, N_rot = 0
!------------------------------------------------------------------------------
SUBROUTINE Case_4(Total_char_freq, Total_harm_osc_freq, HR_params  )
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  REAL(8) ::  Theta_4
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 4
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
 ! EXTERNAL Case_4_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the third harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the fourth harmonic oscillator pseudo-frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 180.     ! Lower bound is 180 cm^-1
  bu(4)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 300.0 !Harm. osc. frequency
  x(2) = 1000.0 !Harm. osc. frequency
  x(3) = 2000.0 !Harm. osc. frequency
  x(4) = 3000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_4_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 4 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 = results(1)
  Theta_2 = results(2)
  Theta_3 = results(3)
  Theta_4 = results(4)
  Total_harm_osc_freq( size(Total_char_freq ) +1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) +2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) +3 ) = Theta_3
  Total_harm_osc_freq( size(Total_char_freq ) +4 ) = Theta_4
  
END SUBROUTINE Case_4
!------------------------------------------------------------------------------
SUBROUTINE Case_4_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 4
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i)) &
               + cv_harm_osc(x(2),t(i)) &
               + cv_harm_osc(x(3),t(i)) &
               + cv_harm_osc(x(4),t(i)) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i)) 
      fj(mcon+i,3) = 2. * diff(i) * d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) * d_Cv_harm_osc(x(4),t(i)) 

    END DO

  END IF

  RETURN

END SUBROUTINE Case_4_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 5:  N_vib = 5, N_rot = 0
!------------------------------------------------------------------------------
SUBROUTINE Case_5(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  REAL(8) ::  Theta_4
  REAL(8) ::  Theta_5
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 5
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_5_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the third harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the fourth harmonic oscillator pseudo-frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 180.     ! Lower bound is 180 cm^-1
  bu(4)  = 4000.    ! Upper bound is 4000 cm^-1

! x(5) corresponds to the fifth harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 250.0 !Harm. osc. frequency
  x(2) = 1000.0 !Harm. osc. frequency
  x(3) = 2000.0 !Harm. osc. frequency
  x(4) = 3000.0 !Harm. osc. frequency
  x(5) = 4000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_5_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 5 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 = results(1)
  Theta_2 = results(2)
  Theta_3 = results(3)
  Theta_4 = results(4)
  Theta_5 = results(5)
  Total_harm_osc_freq( size(Total_char_freq ) +1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) +2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) +3 ) = Theta_3
  Total_harm_osc_freq( size(Total_char_freq ) +4 ) = Theta_4
  Total_harm_osc_freq( size(Total_char_freq ) +5 ) = Theta_5

END SUBROUTINE Case_5
!------------------------------------------------------------------------------
SUBROUTINE Case_5_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 5
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i)) &
               + cv_harm_osc(x(2),t(i)) &
               + cv_harm_osc(x(3),t(i)) &
               + cv_harm_osc(x(4),t(i)) &
               + cv_harm_osc(x(5),t(i)) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i)) 
      fj(mcon+i,3) = 2. * diff(i) * d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) * d_Cv_harm_osc(x(4),t(i)) 
      fj(mcon+i,5) = 2. * diff(i) * d_Cv_harm_osc(x(5),t(i)) 

    END DO

  END IF

  RETURN

END SUBROUTINE Case_5_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 6:  N_vib = 6, N_rot = 0
!------------------------------------------------------------------------------
SUBROUTINE Case_6(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   
  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  REAL(8) ::  Theta_4
  REAL(8) ::  Theta_5
  REAL(8) ::  Theta_6
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_6_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the third harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the fourth harmonic oscillator pseudo-frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 180.     ! Lower bound is 180 cm^-1
  bu(4)  = 4000.    ! Upper bound is 4000 cm^-1

! x(5) corresponds to the fifth harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

! x(6) corresponds to the sixth harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 250.0 !Harm. osc. frequency
  x(2) = 800.0 !Harm. osc. frequency
  x(3) = 1200.0 !Harm. osc. frequency
  x(4) = 2000.0 !Harm. osc. frequency
  x(5) = 2500.0 !Harm. osc. frequency
  x(5) = 4000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_6_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 6 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 = results(1)
  Theta_2 = results(2)
  Theta_3 = results(3)
  Theta_4 = results(4)
  Theta_5 = results(5)
  Theta_6 = results(6)
  Total_harm_osc_freq( size(Total_char_freq ) +1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) +2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) +3 ) = Theta_3
  Total_harm_osc_freq( size(Total_char_freq ) +4 ) = Theta_4
  Total_harm_osc_freq( size(Total_char_freq ) +5 ) = Theta_5
  Total_harm_osc_freq( size(Total_char_freq ) +6 ) = Theta_6

END SUBROUTINE Case_6
!------------------------------------------------------------------------------
SUBROUTINE Case_6_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i)) &
               + cv_harm_osc(x(2),t(i)) &
               + cv_harm_osc(x(3),t(i)) &
               + cv_harm_osc(x(4),t(i)) &
               + cv_harm_osc(x(5),t(i)) &
               + cv_harm_osc(x(6),t(i)) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i)) 
      fj(mcon+i,3) = 2. * diff(i) * d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) * d_Cv_harm_osc(x(4),t(i)) 
      fj(mcon+i,5) = 2. * diff(i) * d_Cv_harm_osc(x(5),t(i)) 
      fj(mcon+i,6) = 2. * diff(i) * d_Cv_harm_osc(x(6),t(i))

    END DO

  END IF

  RETURN

END SUBROUTINE Case_6_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 7:  N_vib >= 7, N_rot = 0
!------------------------------------------------------------------------------
SUBROUTINE Case_7(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(INOUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  REAL(8) ::  Theta_4
  INTEGER :: g_1
  INTEGER :: g_2
  INTEGER :: g_3
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_7_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the first degeneracy
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 0.     ! Lower bound is 0
  bu(2)  = REAL(N_vib) - 1.    ! Upper bound is N_vib - 1

! x(3) corresponds to the second harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the second degeneracy
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 0.     ! Lower bound is 0
  bu(4)  = REAL(N_vib) - 1. -  bu(2)    ! Upper bound is N_vib - 1 - g_1


! x(5) corresponds to the third harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

! x(6) corresponds to the fourth harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 250.0 !Harm. osc. frequency
  x(2) = FLOOR(REAL(N_vib)/3.) !Degeneracy
  x(3) = 800.0 !Harm. osc. frequency
  x(4) = FLOOR(REAL(N_vib)/3.) !Degeneracy
  x(5) = 1600.0 !Harm. osc. frequency
  x(5) = 4000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_7_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 7 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 =      results(1)
  g_1     = NINT(results(2))
  Theta_2 =      results(3)
  g_2     = NINT(results(4))
  Theta_3 =      results(5)
  Theta_4 =      results(6)
  g_3     = (N_vib - 1 - g_1 - g_2)
  Total_harm_osc_freq( size(Total_char_freq ) +1 ) = Theta_1
  DO i = 1, g_1
     Total_harm_osc_freq( size(Total_char_freq ) +1 + i ) = Theta_2
  ENDDO
  DO i = 1, g_2
     Total_harm_osc_freq( size(Total_char_freq ) +1 + g_1 + i ) = Theta_3
  ENDDO
  DO i = 1, g_3
     Total_harm_osc_freq( size(Total_char_freq ) &
          + 1 + g_1 + g_2 + i ) = Theta_4
  ENDDO

 
END SUBROUTINE Case_7
!------------------------------------------------------------------------------
SUBROUTINE Case_7_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i)) &
               + x(2) * cv_harm_osc(x(3),t(i)) &
               + x(4) * cv_harm_osc(x(5),t(i)) &
               + (REAL(N_vib) -1. - x(2) - x(4)) * cv_harm_osc(x(6),t(i)) &
               - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
      fj(mcon+i,2) = 2. * diff(i) * ( Cv_harm_osc(x(3),t(i))  &
                   -  Cv_harm_osc(x(6),t(i)) )
      fj(mcon+i,3) = 2. * diff(i) * x(2) * d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) * ( Cv_harm_osc(x(5),t(i))  &
                   -  Cv_harm_osc(x(6),t(i)) )
      fj(mcon+i,5) = 2. * diff(i) * x(4) * d_Cv_harm_osc(x(5),t(i)) 
      fj(mcon+i,6) = 2. * diff(i) * (REAL(N_vib) -1. - x(2) - x(4)) * &
                     d_Cv_harm_osc(x(6),t(i))

    END DO

  END IF

  RETURN

END SUBROUTINE Case_7_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 8:  N_vib = 0, N_rot = 1
!------------------------------------------------------------------------------

SUBROUTINE Case_8(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 2
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_8_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 50.0 !Hind. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_8_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 8 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1  = results(1)
  nu_1 = results(2)

  HR_params(1,:) = (/ nu_1, V_1 /)

END SUBROUTINE Case_8
!------------------------------------------------------------------------------
SUBROUTINE Case_8_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 2
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2))                               &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))

    END DO

  END IF

  RETURN

END SUBROUTINE Case_8_hd
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 9:  N_vib = 1, N_rot = 1
!------------------------------------------------------------------------------
SUBROUTINE Case_9(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  Theta_1
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 3
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_9_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 150.0 !Hind freq
  x(3) = 1200.0 !Harm. Osc. Freq


!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_9_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 9 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
results = x
V_1     = results(1)
nu_1    = results(2)
Theta_1 = results(3)

Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1

  HR_params(1,:) = (/ nu_1, V_1 /)

END SUBROUTINE Case_9
!------------------------------------------------------------------------------
SUBROUTINE Case_9_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 3
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2)) &
             + cv_harm_osc(x(3),t(i)) & 
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_Cv_harm_osc(x(3),t(i)) 

   END DO

  END IF

  RETURN

END SUBROUTINE Case_9_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 10:  N_vib = 2, N_rot = 1
!------------------------------------------------------------------------------
SUBROUTINE Case_10(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 4
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_10_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the first harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the second harmonic oscillator pseudo-frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 180.     ! Lower bound is 180 cm^-1
  bu(4)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 100.0 !Hind freq
  x(3) = 300.0 !Harm. Osc. Freq
  x(4) = 2000.0 !Harm. Osc. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_10_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 10 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     = results(1)
  nu_1    = results(2)
  Theta_1 = results(3)
  Theta_2 = results(4)

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2

  HR_params(1,:) = (/ nu_1, V_1 /)

END SUBROUTINE Case_10
!------------------------------------------------------------------------------
SUBROUTINE Case_10_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 4
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2)) &
             + cv_harm_osc(x(3),t(i)) & 
             + cv_harm_osc(x(4),t(i)) & 
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) *  d_Cv_harm_osc(x(4),t(i)) 
   END DO

  END IF

  RETURN

END SUBROUTINE Case_10_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 11:  N_vib = 3, N_rot = 1
!------------------------------------------------------------------------------
SUBROUTINE Case_11(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 5
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_11_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the first harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the second harmonic oscillator pseudo-frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 180.     ! Lower bound is 180 cm^-1
  bu(4)  = 4000.    ! Upper bound is 4000 cm^-1

! x(5) corresponds to the third harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 150.0 !Hind freq
  x(3) = 300.0 !Harm. Osc. Freq
  x(4) = 1000.0 !Harm. Osc. Freq
  x(5) = 2000.0 !Harm. Osc. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_11_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 11 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     = results(1)
  nu_1    = results(2)
  Theta_1 = results(3)
  Theta_2 = results(4)
  Theta_3 = results(5)

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) + 3 ) = Theta_3

  HR_params(1,:) = (/ nu_1, V_1 /)

END SUBROUTINE Case_11
!------------------------------------------------------------------------------
SUBROUTINE Case_11_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 5
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2)) &
             + cv_harm_osc(x(3),t(i)) & 
             + cv_harm_osc(x(4),t(i)) & 
             + cv_harm_osc(x(5),t(i)) & 
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) *  d_Cv_harm_osc(x(4),t(i)) 
      fj(mcon+i,5) = 2. * diff(i) *  d_Cv_harm_osc(x(5),t(i)) 
   END DO

  END IF

  RETURN

END SUBROUTINE Case_11_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 12:  N_vib = 4, N_rot = 1
!------------------------------------------------------------------------------
SUBROUTINE Case_12(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  REAL(8) ::  Theta_4
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_12_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the first harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the second harmonic oscillator pseudo-frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 180.     ! Lower bound is 180 cm^-1
  bu(4)  = 4000.    ! Upper bound is 4000 cm^-1

! x(5) corresponds to the third harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

! x(6) corresponds to the fourth harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 150.0 !Hind freq
  x(3) = 300.0 !Harm. Osc. Freq
  x(4) = 1000.0 !Harm. Osc. Freq
  x(5) = 2000.0 !Harm. Osc. Freq
  x(6) = 3000.0 !Harm. Osc. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_12_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 12 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     = results(1)
  nu_1    = results(2)
  Theta_1 = results(3)
  Theta_2 = results(4)
  Theta_3 = results(5)
  Theta_4 = results(6)

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) + 3 ) = Theta_3
  Total_harm_osc_freq( size(Total_char_freq ) + 4 ) = Theta_4

  HR_params(1,:) = (/ nu_1, V_1 /)
          
END SUBROUTINE Case_12
!------------------------------------------------------------------------------
SUBROUTINE Case_12_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2)) &
             + cv_harm_osc(x(3),t(i)) & 
             + cv_harm_osc(x(4),t(i)) & 
             + cv_harm_osc(x(5),t(i)) & 
             + cv_harm_osc(x(6),t(i)) & 
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_Cv_harm_osc(x(3),t(i)) 
      fj(mcon+i,4) = 2. * diff(i) *  d_Cv_harm_osc(x(4),t(i)) 
      fj(mcon+i,5) = 2. * diff(i) *  d_Cv_harm_osc(x(5),t(i)) 
      fj(mcon+i,6) = 2. * diff(i) *  d_Cv_harm_osc(x(6),t(i)) 

   END DO

  END IF

  RETURN

END SUBROUTINE Case_12_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 13:  N_vib >= 5, N_rot = 1
!------------------------------------------------------------------------------
SUBROUTINE Case_13(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  INTEGER ::  g_1
  INTEGER ::  g_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_13_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the first harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the harmonic oscillator degeneracy
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 0.     ! Lower bound is 0
  bu(4)  = REAL(N_vib) -1.    ! Upper bound is N_vib - 1

! x(5) corresponds to the second harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

! x(6) corresponds to the third harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 150.0 !Hind freq
  x(3) = 300.0 !Harm. Osc. Freq
  x(4) = FLOOR(REAL(N_vib)/2.) !Degeneracy
  x(5) = 1000.0 !Harm. Osc. Freq
  x(6) = 3000.0 !Harm. Osc. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_13_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 13 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     =      results(1)
  nu_1    =      results(2)
  Theta_1 =      results(3)
  g_1     = NINT(results(4))
  Theta_2 =      results(5)
  Theta_3 =      results(6)
  g_2     = (N_vib - 1 - g_1)

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  DO i = 1, g_1
     Total_harm_osc_freq( size(Total_char_freq ) +1 + i ) = Theta_2
  ENDDO
  DO i = 1, g_2
     Total_harm_osc_freq( size(Total_char_freq ) +1 + g_1 + i ) = Theta_3
  ENDDO

  HR_params(1,:) = (/ nu_1, V_1 /)

END SUBROUTINE Case_13
!------------------------------------------------------------------------------
SUBROUTINE Case_13_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2)) &
             + cv_harm_osc(x(3),t(i)) & 
             + x(4) * cv_harm_osc(x(5),t(i)) & 
             + (REAL(N_vib) - 1 - x(4) ) * cv_harm_osc(x(6),t(i)) & 
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_Cv_harm_osc(x(3),t(i))
      fj(mcon+i,4) = 2. * diff(i) *  (Cv_harm_osc(x(5),t(i)) & 
                   -  Cv_harm_osc(x(6),t(i)) )
      fj(mcon+i,5) = 2. * diff(i) * x(4) * d_Cv_harm_osc(x(5),t(i)) 
      fj(mcon+i,6) = 2. * diff(i) *  (REAL(N_vib) - 1 - x(4) ) & 
                   * d_Cv_harm_osc(x(6),t(i)) 

   END DO

  END IF

  RETURN

END SUBROUTINE Case_13_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 14:  N_vib = 0, N_rot = 2
!------------------------------------------------------------------------------
SUBROUTINE Case_14(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(INOUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_2
   INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 4
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_14_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the first hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the second hindered-rotor barrier height
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 10.     ! Lower bound is 10 cm^-1
  bu(3)  = 10000.    ! Upper bound is 10000 cm^-1

! x(4) corresponds to the second hindered-rotor frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 40.     ! Lower bound is 40 cm^-1
  bu(4)  = 600.    ! Upper bound is 600 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 80.0 !Hind. Freq
  x(3) = 5000.0 !Barrier Height
  x(4) = 160.0 !Hind. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_14_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 14 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1  = results(1)
  nu_1 = results(2)
  V_2  = results(3)
  nu_2 = results(4)

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)

END SUBROUTINE Case_14
!------------------------------------------------------------------------------
SUBROUTINE Case_14_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 4
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2))      &
             +     cv_hind_rot(x(3),t(i), x(4))      &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_cv_hind_rot(x(3),t(i), x(4))
      fj(mcon+i,4) = 2. * diff(i) *  d_cv_hind_rot_nu(x(3),t(i), x(4))
    END DO

  END IF

  RETURN

END SUBROUTINE Case_14_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 15:  N_vib = 1, N_rot = 2
!------------------------------------------------------------------------------
SUBROUTINE Case_15(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_2
  REAL(8) ::  Theta_1
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 5
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_15_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the first hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the second hindered-rotor barrier height
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 10.     ! Lower bound is 10 cm^-1
  bu(3)  = 10000.    ! Upper bound is 10000 cm^-1

! x(4) corresponds to the second hindered-rotor frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 40.     ! Lower bound is 40 cm^-1
  bu(4)  = 600.    ! Upper bound is 600 cm^-1

! x(5) corresponds to the harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 80.0 !Hind. Freq
  x(3) = 5000.0 !Barrier Height
  x(4) = 160.0 !Hind. Freq
  x(5) = 1200.0 !Harm Osc

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_15_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 15 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     = results(1)
  nu_1    = results(2)
  V_2     = results(3)
  nu_2    = results(4)
  Theta_1 = results(5)

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)

END SUBROUTINE Case_15
!------------------------------------------------------------------------------
SUBROUTINE Case_15_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 5
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2))      &
                 + cv_hind_rot(x(3),t(i), x(4))      &
                 + cv_harm_osc(x(5),t(i)) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_cv_hind_rot(x(3),t(i), x(4))
      fj(mcon+i,4) = 2. * diff(i) *  d_cv_hind_rot_nu(x(3),t(i), x(4))
      fj(mcon+i,5) = 2. * diff(i) *  d_Cv_harm_osc(x(5),t(i))

    END DO

  END IF

  RETURN

END SUBROUTINE Case_15_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 16:  N_vib = 2, N_rot = 2
!------------------------------------------------------------------------------
SUBROUTINE Case_16(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_2
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_16_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the first hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the second hindered-rotor barrier height
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 10.     ! Lower bound is 10 cm^-1
  bu(3)  = 10000.    ! Upper bound is 10000 cm^-1

! x(4) corresponds to the second hindered-rotor frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 40.     ! Lower bound is 40 cm^-1
  bu(4)  = 600.    ! Upper bound is 600 cm^-1

! x(5) corresponds to the first harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

! x(6) corresponds to the second harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 1200.0 !Barrier Height
  x(2) = 80.0 !Hind. Freq
  x(3) = 5000.0 !Barrier Height
  x(4) = 160.0 !Hind. Freq
  x(5) = 300.0 !Harm Osc
  x(6) = 2000.0 !Harm Osc

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_16_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 16 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     = results(1)
  nu_1    = results(2)
  V_2     = results(3)
  nu_2    = results(4)
  Theta_1 = results(5)
  Theta_2 = results(6)

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)

END SUBROUTINE Case_16
!------------------------------------------------------------------------------
SUBROUTINE Case_16_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2))      &
                 + cv_hind_rot(x(3),t(i), x(4))      &
                 + cv_harm_osc(x(5),t(i)) &
                 + cv_harm_osc(x(6),t(i)) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_cv_hind_rot(x(3),t(i), x(4))
      fj(mcon+i,4) = 2. * diff(i) *  d_cv_hind_rot_nu(x(3),t(i), x(4))
      fj(mcon+i,5) = 2. * diff(i) *  d_Cv_harm_osc(x(5),t(i))
      fj(mcon+i,6) = 2. * diff(i) *  d_Cv_harm_osc(x(6),t(i))
    END DO

  END IF

  RETURN

END SUBROUTINE Case_16_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 17:  N_vib = 3, N_rot = 2
!------------------------------------------------------------------------------
SUBROUTINE Case_17(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_2
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_17_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the third harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the first hindered rotor pseudo-barrier
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 10.      ! Has lower bound of 10 cm^-1
  bu(4)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(5) corresponds to the first hindered rotor oscillator frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 40.      ! Has lower bound of 40 cm^-1
  bu(5)  = 180.   ! Has upper bound of 180 cm^-1

! x(6) corresponds to the second hindered rotor pseudo-barrier
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 10.0     ! Has lower bound of 10 cm^-1
  bu(6)  = 10000.   ! Has upper bound of 10,000 cm^-1

!  Define the initial guess for the solution
!  Note that the two pseudo-frequencies should be different so that the
!  Jacobian of x(1) is not automatically zero
  x(1) = 300.0      ! First harm. osc. pseudo-frequency: 500 cm^-1
  x(2) = 1500.0      ! Second harm. osc. pseudo-frequency: 500 cm^-1
  x(3) = 3000.0     ! Third harm. osc. pseudo-frequency: 1500 cm^-1
  x(4) = 1000.0     ! First hind. rot. pseudo-barrier: 500 cm^-1
  x(5) = 100.0      ! First hind. rot. osc. freq. : 100 cm^-1
  x(6) = 2000.0     ! Second hind. rot. pseudo-barrier: 1500 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_17_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )

! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 17'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 = results(1)
  Theta_2 = results(2)
  Theta_3 = results(3)
  V_1     = results(4)
  nu_1    = results(5)
  V_2     = results(6)
  nu_2    = nu_mid

  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) + 3 ) = Theta_3

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)

END SUBROUTINE Case_17
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This subroutine supplies the function to be minimized 
SUBROUTINE Case_17_hd( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference 
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid 
  INTEGER ::  N_vib, N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( Cv_harm_osc(x(1),t(i))                                &
             +  Cv_harm_osc(x(2),t(i))                      &
             +  Cv_harm_osc(x(3),t(i))                      &
             + cv_hind_rot(x(4),t(i),x(5) )                &
             + cv_hind_rot(x(6),t(i), nu_mid)                &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i))
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i))
!     The Jacobian for the third pseudo-frequency
      fj(mcon+i,3) = 2. * diff(i) * d_Cv_harm_osc(x(3),t(i)) 
!     The Jacobian for the hind. rot. degeneracy 
      fj(mcon+i,4) = 2. * diff(i) * d_cv_hind_rot(x(4),t(i), x(5))
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,5) = 2. * diff(i) *  d_cv_hind_rot_nu(x(4),t(i), x(5))
!     The Jacobian for the second pseudo-barrier height
      fj(mcon+i,6) = 2. * diff(i) * d_cv_hind_rot(x(6),t(i), nu_mid)

    END DO

  END IF

  RETURN
END SUBROUTINE Case_17_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 18:  N_vib >= 4, N_rot = 2
!------------------------------------------------------------------------------
SUBROUTINE Case_18(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params

  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_2
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER ::  g_1
  INTEGER ::  g_2
  INTEGER ::  i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_18_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the degeneracy of the first harmonic oscillator
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 0.0      ! Lower bound is zero
  bu(1)  = REAL(N_vib) ! Upper bound is the number of unknown harm. osc.

! x(2) corresponds to the first harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the second harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the first hindered rotor pseudo-barrier
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 10.      ! Has lower bound of 10 cm^-1
  bu(4)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(5) corresponds to the first hindered rotor oscillator frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 40.      ! Has lower bound of 40 cm^-1
  bu(5)  = 180.   ! Has upper bound of 180 cm^-1

! x(6) corresponds to the second hindered rotor pseudo-barrier
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 10.0     ! Has lower bound of 10 cm^-1
  bu(6)  = 10000.   ! Has upper bound of 10,000 cm^-1

!  Define the initial guess for the solution
!  Note that the two pseudo-frequencies should be different so that the
!  Jacobian of x(1) is not automatically zero
  x(1) = floor(REAL(N_vib)/2.) !Harm. osc. degen.:  1/2 the # of unknown osc
  x(2) = 300.0      ! First harm. osc. pseudo-frequency: 500 cm^-1
  x(3) = 3000.0     ! Second harm. osc. pseudo-frequency: 1500 cm^-1
  x(4) = 1000.0     ! First hind. rot. pseudo-barrier: 500 cm^-1
  x(5) = 100.0      ! First hind. rot. osc. freq. : 100 cm^-1
  x(6) = 3000.0     ! Second hind. rot. pseudo-barrier: 1500 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_18_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )

! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 18'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  g_1     = NINT(results(1))
  Theta_1 =      results(2)
  Theta_2 =      results(3)
  V_1     =      results(4)
  nu_1    =      results(5)
  V_2     =      results(6)
  nu_2    = nu_mid
  g_2     = (N_vib - g_1)

  DO i = 1,g_1
     Total_harm_osc_freq( size(Total_char_freq ) + i ) = Theta_1
  ENDDO
  DO i = 1, g_2
     Total_harm_osc_freq( size(Total_char_freq ) + g_1 + i ) = Theta_2
  ENDDO

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)

END SUBROUTINE Case_18
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This subroutine supplies the function to be minimized 
SUBROUTINE Case_18_hd( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference 
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib, N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( x(1) * cv_harm_osc(x(2),t(i))                                &
             + (N_vib -x(1) ) *   Cv_harm_osc(x(3),t(i))                      &
             + cv_hind_rot(x(4),t(i),x(5) )                &
             + cv_hind_rot(x(6),t(i), nu_mid)                &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *                                           &
                   ( Cv_harm_osc(x(2),t(i)) - Cv_harm_osc(x(3),t(i)) )
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) * ( x(1) * d_Cv_harm_osc(x(2),t(i)) )
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,3) =  2. * diff(i) *                                          &
                   ( (N_vib - x(1) ) * d_Cv_harm_osc(x(3),t(i)) )
!     The Jacobian for the hind. rot. degeneracy 
      fj(mcon+i,4) = 2. * diff(i) * d_cv_hind_rot(x(4),t(i), x(5))
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,5) = 2. * diff(i) *  d_cv_hind_rot_nu(x(4),t(i), x(5))
!     The Jacobian for the second pseudo-barrier height
      fj(mcon+i,6) = 2. * diff(i) * d_cv_hind_rot(x(6),t(i), nu_mid)

    END DO

  END IF

  RETURN
END SUBROUTINE Case_18_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 19:  N_vib = 0, N_rot = 3
!------------------------------------------------------------------------------
SUBROUTINE Case_19(Total_char_freq, Total_harm_osc_freq, HR_params )
 
  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(INOUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  V_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_2
  REAL(8) ::  nu_3
  REAL(8) ::  V_3
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_19_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first hindered-rotor barrier height
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.     ! Lower bound is 10 cm^-1
  bu(1)  = 10000.    ! Upper bound is 10000 cm^-1

! x(2) corresponds to the first hindered-rotor frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the second hindered-rotor barrier height
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 10.     ! Lower bound is 10 cm^-1
  bu(3)  = 10000.    ! Upper bound is 10000 cm^-1

! x(4) corresponds to the second hindered-rotor frequency
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 40.     ! Lower bound is 40 cm^-1
  bu(4)  = 600.    ! Upper bound is 600 cm^-1

! x(5) corresponds to the third hindered-rotor barrier height
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 10.     ! Lower bound is 10 cm^-1
  bu(5)  = 10000.    ! Upper bound is 10000 cm^-1

! x(6) corresponds to the third hindered-rotor frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 40.     ! Lower bound is 40 cm^-1
  bu(6)  = 600.    ! Upper bound is 600 cm^-1

!  Define the initial guess for the solution
  x(1) = 3000.0 !Barrier Height
  x(2) = 80.0 !Hind. Freq
  x(3) = 1000.0 !Barrier Height
  x(4) = 150.0 !Hind. Freq
  x(5) = 3000.0 !Barrier Height
  x(6) = 300.0 !Hind. Freq

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed

!  Call the program

  CALL dqed ( Case_19_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 19 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1  = results(1)
  nu_1 = results(2)
  V_2  = results(3)
  nu_2 = results(4)
  V_3  = results(5)
  nu_3 = results(6)

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)
  HR_params(3,:) = (/ nu_3, V_3 /)

END SUBROUTINE Case_19
!------------------------------------------------------------------------------
SUBROUTINE Case_19_hd( x, fj, ldfj, igo, iopt, ropt )


  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = (   cv_hind_rot(x(1),t(i), x(2))      &
             +     cv_hind_rot(x(3),t(i), x(4))      &
             +     cv_hind_rot(x(5),t(i), x(6))      &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *  d_cv_hind_rot(x(1),t(i), x(2))
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(1),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_cv_hind_rot(x(3),t(i), x(4))
      fj(mcon+i,4) = 2. * diff(i) *  d_cv_hind_rot_nu(x(3),t(i), x(4))
      fj(mcon+i,5) = 2. * diff(i) *  d_cv_hind_rot(x(5),t(i), x(6))
      fj(mcon+i,6) = 2. * diff(i) *  d_cv_hind_rot_nu(x(5),t(i), x(6))
    END DO

  END IF

  RETURN

END SUBROUTINE Case_19_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 20:  N_vib = 1, N_rot = 3
!------------------------------------------------------------------------------
SUBROUTINE Case_20(Total_char_freq, Total_harm_osc_freq, HR_params  )

  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: HR_params

  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  nu_3
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  V_3
  REAL(8) ::  Theta_1
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_20_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the harm. osc. freq. of  the first hindered rotor
  ind(2) = 3     ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40
  bu(2)  = 600.    ! Upper bound is 600

! x(3) corresponds to the first pseudo-barrier height
  ind(3) = 3       ! Has both an upper and lower bound
  bl(3)  = 10.     ! Lower bound is 10 cm^-1
  bu(3)  = 10000.  ! Upper bound is 10000 cm^-1

! x(4) corresponds to the harm. osc. freq. of  the second hindered rotor
  ind(4) = 3     ! Has both an upper and lower bound
  bl(4)  = 40.     ! Lower bound is 40
  bu(4)  = 600.    ! Upper bound is 600

! x(5) corresponds to the second pseudo-barrier height
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 10.0     ! Lower bound is 10 cm^-1
  bu(5)  = 10000.0  ! Upper bound is 10000 cm^-1

! x(6) corresponds to the third pseudo-barrier height
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 10.0     ! Lower bound is 10 cm^-1
  bu(6)  = 10000.0  ! Upper bound is 10000 cm^-1

!  Define the initial guess for the solution
  x(1) = 300.0 !Harm. osc. frequency
  x(2) = 80.0 !Hind. Freq
  x(3) = 1000.0 !Barrier Height
  x(4) = 250.0 !Hind. Freq
  x(5) = 2000.0 !Barrier Height
  x(6) = 3000.0 !Barrier Height


!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_20_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 20 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x

  Theta_1  = results(1)
  nu_1     = results(2)
  V_1      = results(3)
  nu_2     = results(4)
  V_2      = results(5)
  V_3      = results(6)
  nu_3     = nu_mid

  Total_harm_osc_freq( size(Total_char_freq ) + 1) = Theta_1

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)
  HR_params(3,:) = (/ nu_3, V_3 /)

END SUBROUTINE Case_20
!------------------------------------------------------------------------------
SUBROUTINE Case_20_hd ( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i))                             &
             + cv_hind_rot(x(3),t(i), x(2))                  &
             + cv_hind_rot(x(5),t(i), x(4)) &
             + cv_hind_rot(x(6),t(i), nu_mid) &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) *  d_cv_hind_rot_nu(x(3),t(i), x(2))
      fj(mcon+i,3) = 2. * diff(i) *  d_cv_hind_rot(x(3),t(i), x(2))
      fj(mcon+i,4) = 2. * diff(i) *  d_cv_hind_rot_nu(x(5),t(i), x(4))
      fj(mcon+i,5) = 2. * diff(i) *  d_cv_hind_rot(x(5),t(i), x(4))
      fj(mcon+i,6) = 2. * diff(i) *  d_cv_hind_rot(x(6),t(i), nu_mid)

    END DO

  END IF

  RETURN
END SUBROUTINE Case_20_hd

!------------------------------------------------------------------------------
! CASE 21:  N_vib = 2, N_rot = 3
!------------------------------------------------------------------------------
SUBROUTINE Case_21(Total_char_freq, Total_harm_osc_freq, HR_params )

  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  nu_3
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  V_3
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_21_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the harm. osc. freq. of the first hindered rotor
  ind(3) = 3      ! Has both an upper and lower bound
  bl(3)  = 40.     ! Lower bound is 40
  bu(3)  = 600.    ! Upper bound is 600

! x(4) corresponds to the first pseudo-barrier height
  ind(4) = 3       ! Has both an upper and lower bound
  bl(4)  = 10.     ! Lower bound is 10 cm^-1
  bu(4)  = 10000.  ! Upper bound is 10000 cm^-1

! x(5) corresponds to the second pseudo-barrier height
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 10.0     ! Lower bound is 10 cm^-1
  bu(5)  = 10000.0  ! Upper bound is 10000 cm^-1

! x(6) corresponds to the third harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 10.     ! Lower bound is 180 cm^-1
  bu(6)  = 10000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 300.0 !Harm. osc. frequency
  x(2) = 3000.0 !Harm. osc. frequency
  x(3) = 150. !Hind Rot. Harm osc. freq. 
  x(4) = 1000. ! First pseudo-barrier height: 500 cm^-1
  x(5) = 2000. ! Second pseudo-barrier height: 3000 cm^-1
  x(6) = 3000. ! Third pseudo-barrier height: 3000 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_21_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 21 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1  = results(1)
  Theta_2  = results(2)
  nu_1     = results(3)
  V_1      = results(4)
  V_2      = results(5)
  V_3      = results(6)
  nu_2     = nu_low
  nu_3     = nu_high

  Total_harm_osc_freq( size(Total_char_freq ) + 1) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2) = Theta_2

  HR_params(1,:) = (/ nu_1, V_1 /)
  HR_params(2,:) = (/ nu_2, V_2 /)
  HR_params(3,:) = (/ nu_3, V_3 /)
  
END SUBROUTINE Case_21
!------------------------------------------------------------------------------
SUBROUTINE Case_21_hd ( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high,nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i))                             &
               + cv_harm_osc(x(2),t(i))                             &
               + cv_hind_rot(x(4),t(i), x(3))                       &
               + cv_hind_rot(x(5),t(i), nu_low)                       &
               + cv_hind_rot(x(6),t(i), nu_high)                       &
             - cp_difference(i) ) 

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the first harm. osc. freq 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
!     The Jacobian for the second harm. osc. freq 
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i))
!     The Jacobian for the first hind rot harm. osc. freq
      fj(mcon+i,3) = 2. * diff(i) *  d_cv_hind_rot_nu(x(4),t(i), x(3))
!     The Jacobian for the first hind rot barrier height
      fj(mcon+i,4) = 2. * diff(i) *  d_cv_hind_rot(x(4),t(i), x(3))
!     The Jacobian for the second hind rot barrier height
      fj(mcon+i,5) = 2. * diff(i) *  d_cv_hind_rot(x(5),t(i), nu_low)
!     The Jacobian for the third hind rot barrier height
      fj(mcon+i,6) = 2. * diff(i) *  d_cv_hind_rot(x(6),t(i), nu_high)


    END DO

  END IF

  RETURN
END SUBROUTINE Case_21_hd
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! CASE 22:  N_vib = 3, N_rot >= 3
!------------------------------------------------------------------------------
SUBROUTINE Case_22(Total_char_freq, Total_harm_osc_freq, HR_params )

  USE heat_capacity_functions
   

  IMPLICIT NONE
  
! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  nu_3
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  V_3
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  REAL(8) ::  Theta_3
  INTEGER ::  g_1
  INTEGER ::  g_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_22_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot
 
 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the degeneracy of the first hindered rotor
  ind(2) = 3      ! Has both an upper and lower bound
  bl(2)  = 0.     ! Lower bound is 0
  bu(2)  = REAL(N_rot)    ! Upper bound is N_rot

! x(3) corresponds to the first pseudo-barrier height
  ind(3) = 3       ! Has both an upper and lower bound
  bl(3)  = 10.     ! Lower bound is 10 cm^-1
  bu(3)  = 10000.  ! Upper bound is 10000 cm^-1

! x(4) corresponds to the second pseudo-barrier height
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 10.0     ! Lower bound is 10 cm^-1
  bu(4)  = 10000.0  ! Upper bound is 10000 cm^-1

! x(5) corresponds to the second harmonic oscillator pseudo-frequency
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 180.     ! Lower bound is 180 cm^-1
  bu(5)  = 4000.    ! Upper bound is 4000 cm^-1

! x(6) corresponds to the third harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
  x(1) = 300.0 !Harm. osc. frequency
  x(2) = FLOOR(REAL(N_rot)/2.) !Degeneracy of first pseudo-frequency: 1/2 N_rot
  x(3) = 500.0 ! First pseudo-barrier height: 500 cm^-1
  x(4) = 3000. ! Second pseudo-barrier height: 3000 cm^-1
  x(5) = 1500.0 !Harm. osc. frequency
  x(6) = 3000.0 !Harm. osc. frequency

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_22_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 22 Baby Yea!'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1  =      results(1)
  g_1      = NINT(results(2) )
  V_1      =      results(3)
  V_2      =      results(4)
  Theta_2  =      results(5)
  Theta_3  =      results(6)
  nu_1     =      nu_low
  nu_2     =      nu_high
  g_2      = (N_rot - g_1)
  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2
  Total_harm_osc_freq( size(Total_char_freq ) + 3 ) = Theta_3

DO i = 1,g_1
   HR_params(i,:) = (/ nu_1, V_1 /)
ENDDO

DO i = (g_1+1),N_rot
   HR_params(i,:) = (/ nu_2, V_2 /)
ENDDO

END SUBROUTINE Case_22
!------------------------------------------------------------------------------
SUBROUTINE Case_22_hd ( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
 REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib
  INTEGER ::  N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_harm_osc(x(1),t(i))                             &
             + x(2) * cv_hind_rot(x(3),t(i), nu_low)                  &
             + (REAL(N_rot) - x(2) ) * cv_hind_rot(x(4),t(i), nu_high) &
             + cv_harm_osc(x(5),t(i))                               &
             + cv_harm_osc(x(6),t(i))                               &
             - cp_difference(i) ) 

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i)) 
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i)                                      &
                   * ( cv_hind_rot(x(3),t(i), nu_low)                    & 
                   -   cv_hind_rot(x(4),t(i), nu_high) )
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,3) = 2. * diff(i) * ( x(2) * d_cv_hind_rot(x(3),t(i), nu_low) )
!     The Jacobian for the second pseudo-barrier height
      fj(mcon+i,4) = 2. * diff(i) *  &
                        (REAL(N_rot) - x(2) ) * d_cv_hind_rot(x(4),t(i), nu_high)
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,5) = 2. * diff(i) * d_Cv_harm_osc(x(5),t(i)) 
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,6) = 2. * diff(i) * d_Cv_harm_osc(x(6),t(i)) 
    END DO

  END IF

  RETURN
END SUBROUTINE Case_22_hd
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! CASE 23:  N_vib >= 4, N_rot >= 3
!------------------------------------------------------------------------------
SUBROUTINE Case_23(Total_char_freq, Total_harm_osc_freq, HR_params )

  USE heat_capacity_functions
   

  IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER ::  g_1
  INTEGER ::  g_2
  INTEGER ::  g_3
  INTEGER ::  g_4
  INTEGER ::  i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_23_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the degeneracy of the first harmonic oscillator
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 0.0      ! Lower bound is zero
  bu(1)  = REAL(N_vib) ! Upper bound is the number of unknown harm. osc.

! x(2) corresponds to the first harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the second harmonic oscillator pseudo-frequency
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 180.     ! Lower bound is 180 cm^-1
  bu(3)  = 4000.    ! Upper bound is 4000 cm^-1

! x(4) corresponds to the degeneracy of the first hindered rotor
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 0.0      ! Lower bound is zero
  bu(4)  = REAL(N_rot)      ! Upper bound is the number of hindered rotors

! x(5) corresponds to the first hindered rotor pseudo-barrier
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 10.      ! Has lower bound of 10 cm^-1
  bu(5)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(6) corresponds to the second hindered rotor pseudo-barrier
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 10.0     ! Has lower bound of 10 cm^-1
  bu(6)  = 10000.   ! Has upper bound of 10,000 cm^-1

!  Define the initial guess for the solution
!  Note that the two pseudo-frequencies should be different so that the
!  Jacobian of x(1) is not automatically zero
  x(1) = floor(REAL(N_vib)/2.) !Harm. osc. degen.:  1/2 the # of unknown osc
  x(2) = 300.0      ! First harm. osc. pseudo-frequency: 500 cm^-1
  x(3) = 3000.0     ! Second harm. osc. pseudo-frequency: 1500 cm^-1
  x(4) = floor(REAL(N_rot)/2.) ! Hind. rot. degen.:  1/2 the # of hind. rot.
  x(5) = 500.0      ! First hind. rot. pseudo-barrier: 500 cm^-1
  x(6) = 4000.0     ! Second hind. rot. pseudo-barrier: 1500 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_23_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )

! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 23'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  g_1      = NINT(results(1) )
  Theta_1  =      results(2)
  Theta_2  =      results(3)
  g_2      = NINT(results(4) )
  V_1      =      results(5)
  V_2      =      results(6)
  nu_1     =      nu_low
  nu_2     =      nu_high
  g_3      = (N_vib - g_1)
  g_4      = (N_rot - g_2)

  DO i = 1,g_1
     Total_harm_osc_freq( size(Total_char_freq ) + i ) = Theta_1
  ENDDO
  DO i = 1, g_3
     Total_harm_osc_freq( size(Total_char_freq ) + g_1 + i ) = Theta_2
  ENDDO

  DO i = 1,g_2
     HR_params(i,:) = (/ nu_1, V_1 /)
  ENDDO

  DO i = (g_2+1),N_rot
     HR_params(i,:) = (/ nu_2, V_2 /)
  ENDDO

END SUBROUTINE Case_23
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This subroutine supplies the function to be minimized 
SUBROUTINE Case_23_hd( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference 
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib, N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( x(1) * cv_harm_osc(x(2),t(i))                                &
             + (REAL(N_vib) -x(1) ) *   Cv_harm_osc(x(3),t(i))                &
             + x(4) * cv_hind_rot(x(5),t(i), nu_low)                            &
             + (REAL(N_rot) - x(4) ) * cv_hind_rot(x(6),t(i), nu_high)           &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) *                                           &
                   ( Cv_harm_osc(x(2),t(i)) - Cv_harm_osc(x(3),t(i)) )
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) * ( x(1) * d_Cv_harm_osc(x(2),t(i)) )
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,3) =  2. * diff(i) *                                          &
                   ( (REAL(N_vib) - x(1) ) * d_Cv_harm_osc(x(3),t(i)) )
!     The Jacobian for the hind. rot. degeneracy 
      fj(mcon+i,4) = 2. * diff(i) *                                           &
                   ( cv_hind_rot(x(5),t(i), nu_low) -                           &
                     cv_hind_rot(x(6),t(i), nu_high) )
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,5) = 2. * diff(i) * ( x(4) *d_cv_hind_rot( x(5),t(i),  nu_low) )
!     The Jacobian for the second pseudo-barrier height
      fj(mcon+i,6) = 2. * diff(i) *                                           &
                   ( (REAL(N_rot) - x(4) ) * d_cv_hind_rot(x(6),t(i), nu_high) )

    END DO

  END IF

  RETURN
END SUBROUTINE Case_23_hd

!------------------------------------------------------------------------------
! CASE 24:  N_vib = 0, N_rot >= 4
!------------------------------------------------------------------------------
SUBROUTINE Case_24(Total_char_freq, Total_harm_osc_freq, HR_params )

  USE heat_capacity_functions
   

  IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(INOUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  nu_3
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  V_3
  INTEGER ::  g_1
  INTEGER ::  g_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_24_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first hindered rotor pseudo-barrier
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.      ! Has lower bound of 10 cm^-1
  bu(1)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(2) corresponds to the first hindered rotor oscillator freq
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the degeneracy of the first hindered rotor
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 0.0      ! Lower bound is zero
  bu(3)  = REAL(N_rot) - 1.0 ! Upper bound is the number of hind. rot - 1

! x(4) corresponds to the second hindered rotor pseudo-barrier
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 10.0     ! Has lower bound of 10 cm^-1
  bu(4)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(5) corresponds to the third hindered rotor pseudo-barrier
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 10.      ! Has lower bound of 10 cm^-1
  bu(5)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(6) corresponds to the second hindered rotor oscillator freq
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 40.     ! Lower bound is 40 cm^-1
  bu(6)  = 600.    ! Upper bound is 600 cm^-1

!  Define the initial guess for the solution
!  Note that the two pseudo-frequencies should be different so that the
!  Jacobian of x(1) is not automatically zero
  x(1) = 3000.0      ! First hind. rot. pseudo-barrier: 1000 cm^-1
  x(2) = 100.0      ! First hind. rot osc. frequency: 150 cm^-1
  x(3) = floor(REAL(N_rot)/2.) -1.! Hind. rot. degen.:  1/2 the # of hind. rot.
  x(4) = 500.0      ! second hind. rot. pseudo-barrier: 500 cm^-1
  x(5) = 1500.0     ! third hind. rot. pseudo-barrier: 1500 cm^-1
  x(2) = 100.0      ! First hind. rot osc. frequency: 150 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_24_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )

! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 24'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1  =      results(1)
  nu_1 =      results(2)
  g_1  = NINT(results(3) )
  V_2  =      results(4)
  V_3  =      results(5)
  nu_2 =      results(6)
  nu_3 = nu_mid
  g_2  = (N_rot - 1 - g_1)

  IF (g_1==1) THEN
!     HR_rho_1_total = HR_rho_2
  ELSE
!     call convolution_rho(HR_rho_2, HR_rho_2, (g_1-1), delta_nu, HR_rho_1_total)
  ENDIF
  IF (g_2==1) THEN
!     HR_rho_2_total = HR_rho_3
  ELSE
!     call convolution_rho(HR_rho_3, HR_rho_3, (g_2-1), delta_nu, HR_rho_2_total)
  ENDIF


  HR_params(1,:) = (/ nu_1, V_1 /)
  
  DO i = 2,(g_1+1)
     HR_params(i,:) = (/ nu_2, V_2 /)
  ENDDO

  DO i = (g_1+2),N_rot
     HR_params(i,:) = (/ nu_3, V_3 /)
  ENDDO


END SUBROUTINE Case_24
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This subroutine supplies the function to be minimized 
SUBROUTINE Case_24_hd( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference 
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib, N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_hind_rot(x(1),t(i), x(2))  &
             + x(3) * cv_hind_rot(x(4),t(i), x(6))  &
             + (REAL(N_rot) - 1. - x(3) ) * cv_hind_rot(x(5),t(i), nu_mid)  &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_cv_hind_rot( x(1),t(i),  x(2)) 
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) * d_cv_hind_rot_nu( x(1),t(i),  x(2)) 
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,3) =  2. * diff(i) *                                   &
         ( cv_hind_rot( x(4),t(i),  x(6)) - cv_hind_rot( x(5),t(i),  nu_mid) ) 
!     The Jacobian for the hind. rot. degeneracy 
      fj(mcon+i,4) = 2. * diff(i) * x(3) *  d_cv_hind_rot( x(4),t(i), x(6)) 
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,5) = 2. * diff(i) *  (REAL(N_rot) - 1. - x(3) )    &
                   *  d_cv_hind_rot( x(5),t(i), nu_mid)
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,6) = 2. * diff(i) * d_cv_hind_rot_nu( x(4),t(i),  x(6)) 

    END DO

  END IF

  RETURN
END SUBROUTINE Case_24_hd


!------------------------------------------------------------------------------
! CASE 25:  N_vib = 1, N_rot >= 4
!------------------------------------------------------------------------------
SUBROUTINE Case_25(Total_char_freq, Total_harm_osc_freq, HR_params )

  USE heat_capacity_functions
   

  IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  nu_3
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  V_3
  REAL(8) ::  Theta_1
  INTEGER ::  g_1
  INTEGER ::  g_2
  INTEGER :: i

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_25_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first hindered rotor pseudo-barrier
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 10.      ! Has lower bound of 10 cm^-1
  bu(1)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(2) corresponds to the first hindered rotor oscillator freq
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 40.     ! Lower bound is 40 cm^-1
  bu(2)  = 600.    ! Upper bound is 600 cm^-1

! x(3) corresponds to the degeneracy of the first hindered rotor
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 0.0      ! Lower bound is zero
  bu(3)  = REAL(N_rot) - 1.0 ! Upper bound is the number of hind. rot - 1

! x(4) corresponds to the second hindered rotor pseudo-barrier
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 10.0     ! Has lower bound of 10 cm^-1
  bu(4)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(5) corresponds to the third hindered rotor pseudo-barrier
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 10.      ! Has lower bound of 10 cm^-1
  bu(5)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(6) corresponds to the first harmonic oscillator pseudo-frequency
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 180.     ! Lower bound is 180 cm^-1
  bu(6)  = 4000.    ! Upper bound is 4000 cm^-1

!  Define the initial guess for the solution
!  Note that the two pseudo-frequencies should be different so that the
!  Jacobian of x(1) is not automatically zero
  x(1) = 3000.0      ! First hind. rot. pseudo-barrier: 1000 cm^-1
  x(2) = 150.0      ! First hind. rot osc. frequency: 150 cm^-1
  x(3) = floor(REAL(N_rot)/2.) -1.! Hind. rot. degen.:  1/2 the # of hind. rot.
  x(4) = 500.0      ! second hind. rot. pseudo-barrier: 500 cm^-1
  x(5) = 1500.0     ! third hind. rot. pseudo-barrier: 1500 cm^-1
  x(6) = 300.0      ! First harm. osc.: 1200 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_25_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )

! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 25'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  V_1     =      results(1)
  nu_1    =      results(2)
  g_1     = NINT(results(3) )
  V_2     =      results(4)
  V_3     =      results(5)
  Theta_1 =      results(6)
  nu_2    = nu_low
  nu_3    = nu_high
  g_2     = (N_rot - 1 - g_1)
  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1

  HR_params(1,:) = (/ nu_1, V_1 /)
  
  DO i = 2,(g_1+1)
     HR_params(i,:) = (/ nu_2, V_2 /)
  ENDDO

  DO i = (g_1+2),N_rot
     HR_params(i,:) = (/ nu_3, V_3 /)
  ENDDO

END SUBROUTINE Case_25
!------------------------------------------------------------------------------
! This subroutine supplies the function to be minimized 
SUBROUTINE Case_25_hd( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference 
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib, N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( cv_hind_rot(x(1),t(i), x(2))                &
             +   x(3) * cv_hind_rot(x(4),t(i), nu_low)         &
             +   (REAL(N_rot) - 1. - x(3) ) * cv_hind_rot(x(5),t(i), nu_high)  &
             +   Cv_harm_osc(x(6),t(i))                    &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_cv_hind_rot( x(1),t(i),  x(2)) 
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) * d_cv_hind_rot_nu( x(1),t(i),  x(2)) 
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,3) =  2. * diff(i) *                                   &
         ( cv_hind_rot( x(4),t(i),  nu_low) - cv_hind_rot( x(5),t(i),  nu_high) ) 
!     The Jacobian for the hind. rot. degeneracy 
      fj(mcon+i,4) = 2. * diff(i) * x(3) *  d_cv_hind_rot( x(4),t(i), nu_low) 
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,5) = 2. * diff(i) *  (REAL(N_rot) - 1. - x(3) )    &
                   *  d_cv_hind_rot( x(5),t(i), nu_high)
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,6) = 2. * diff(i) *  d_Cv_harm_osc(x(6),t(i))   
    END DO

  END IF

  RETURN
END SUBROUTINE Case_25_hd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! CASE 26:  N_vib = 2, N_rot >= 4
!------------------------------------------------------------------------------
SUBROUTINE Case_26(Total_char_freq, Total_harm_osc_freq, HR_params )

  USE heat_capacity_functions

  IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference

! These variables are specific to this case
  REAL(8), INTENT(IN), DIMENSION(:) :: Total_char_freq
  REAL(8), INTENT(OUT), DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), INTENT(OUT), DIMENSION(:,:) :: HR_params
  REAL(8) ::  nu_1
  REAL(8) ::  nu_2
  REAL(8) ::  V_1
  REAL(8) ::  V_2
  REAL(8) ::  Theta_1
  REAL(8) ::  Theta_2
  INTEGER ::  g_1
  INTEGER ::  g_2

! These variables are required by the nonlinear solver
  INTEGER, PARAMETER :: liwork = 103
  INTEGER, PARAMETER :: lwork = 785
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  INTEGER, PARAMETER :: ldfj = mcon + mequa
  REAL(8) :: bl(nvars+mcon)
  REAL(8) :: bu(nvars+mcon)
  REAL(8) :: fj(ldfj,nvars+1)
  REAL(8) :: fnorm
  INTEGER :: igo
  INTEGER :: ind(nvars+mcon)
  INTEGER :: iopt(24)
  INTEGER :: iwork(liwork)
  REAL(8) :: ropt(1)
  REAL(8) :: work(lwork)
  REAL(8) :: x(nvars)
  REAL(8) :: results(nvars)
!  EXTERNAL Case_26_hd

! These variables are specific to the molecule
  REAL(8) ::  nu_low
  REAL(8) ::  nu_high
  REAL(8) ::  nu_mid
  INTEGER :: N_vib, N_rot, i

 COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps
! Define the contraints on the target variables 

! x(1) corresponds to the first harmonic oscillator pseudo-frequency
  ind(1) = 3        ! Has both an upper and lower bound
  bl(1)  = 180.     ! Lower bound is 180 cm^-1
  bu(1)  = 4000.    ! Upper bound is 4000 cm^-1

! x(2) corresponds to the second harmonic oscillator pseudo-frequency
  ind(2) = 3        ! Has both an upper and lower bound
  bl(2)  = 180.     ! Lower bound is 180 cm^-1
  bu(2)  = 4000.    ! Upper bound is 4000 cm^-1

! x(3) corresponds to the degeneracy of the first hindered rotor
  ind(3) = 3        ! Has both an upper and lower bound
  bl(3)  = 0.0      ! Lower bound is zero
  bu(3)  = REAL(N_rot)! Upper bound is the number of hind. rot

! x(4) corresponds to the first hindered rotor pseudo-barrier
  ind(4) = 3        ! Has both an upper and lower bound
  bl(4)  = 10.0     ! Has lower bound of 10 cm^-1
  bu(4)  = 10000.   ! Has upper bound of 10,000 cm^-1

! x(5) corresponds to the first hindered rotor pseudo-barrier\osc. freq.
  ind(5) = 3        ! Has both an upper and lower bound
  bl(5)  = 40.      ! Has lower bound of 40 cm^-1
  bu(5)  = 600.   ! Has upper bound of 600 cm^-1

! x(6) corresponds to the second harmonic oscillator pseudo-barrier
  ind(6) = 3        ! Has both an upper and lower bound
  bl(6)  = 10.     ! Lower bound is 10 cm^-1
  bu(6)  = 10000.    ! Upper bound is 10000 cm^-1

!  Define the initial guess for the solution
!  Note that the two pseudo-frequencies should be different so that the
!  Jacobian of x(1) is not automatically zero
  x(1) = 300.0      !First harm osc. frequency: 800 cm^-1
  x(2) = 3000.0      ! Second harm osc. frequency: 1200 cm^-1
  x(3) = floor(REAL(N_rot)/2.)! Hind. rot. degen.:  1/2 the # of hind. rot.
  x(4) = 1000.0      ! First hind. rot. pseudo-barrier: 1000 cm^-1
  x(5) = 100.0     ! First hind. rot.osc. freq: 100 cm^-1
  x(6) = 4000.0      ! Second hind. rot. pseudo-barrier: 2000 cm^-1

!  Tell how much storage we gave the solver.
  iwork(1) = lwork
  iwork(2) = liwork

!  Additional solver options
  iopt(1)=4         ! Set the option to change the value of TOLF
  iopt(2)=1         ! Where in ROPT to look for the new value
  ropt(1)=1.E-9     ! New value for TOLF
  iopt(3)=2         ! Change the number of interations
  iopt(4)=100000    ! Maximum number of iterations
  iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
  iopt(6)=1         ! Forces a full model step
  iopt(7)=99        ! No further options are changed


!  Call the program

  CALL dqed ( Case_26_hd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, &
    fnorm, igo, iopt, ropt, iwork, work )

! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) ' CASE 26'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,i6)' ) '  Output flag from DQED, IGO = ', igo
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) '  Computed X:'
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(g14.8)' ) x(1:nvars)
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a,g14.6)' ) '  L2 norm of the residual, FNORM = ', fnorm
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'Expected X:'
 
! WRITE ( *, '(a)' ) ' '
! WRITE ( *, '(a)' ) 'DQED_PRB2'
! WRITE ( *, '(a)' ) '  Normal end of execution.'

! WRITE ( *, '(a)' ) ' '
  results = x
  Theta_1 =      results(1)
  Theta_2 =      results(2)
  g_1     = NINT(results(3) )
  V_1     =      results(4)
  nu_1    =      results(5)
  V_2     =      results(6)
  nu_2    = nu_mid
  g_2     = (N_rot  - g_1)
  Total_harm_osc_freq( size(Total_char_freq ) + 1 ) = Theta_1
  Total_harm_osc_freq( size(Total_char_freq ) + 2 ) = Theta_2

 IF (g_1==1) THEN
!     HR_rho_1_total = HR_rho_1
  ELSE
!     call convolution_rho(HR_rho_1, HR_rho_1, (g_1-1), delta_nu, HR_rho_1_total)
  ENDIF
  IF (g_2==1) THEN
!     HR_rho_2_total = HR_rho_2
  ELSE
!     call convolution_rho(HR_rho_2, HR_rho_2, (g_2-1), delta_nu, HR_rho_2_total)
  ENDIF
  
  DO i = 1,(g_1)
     HR_params(i,:) = (/ nu_1, V_1 /)
  ENDDO

  DO i = (g_1+1),N_rot
     HR_params(i,:) = (/ nu_2, V_2 /)
  ENDDO


END SUBROUTINE Case_26
!------------------------------------------------------------------------------
! This subroutine supplies the function to be minimized 
SUBROUTINE Case_26_hd( x, fj, ldfj, igo, iopt, ropt )

  USE heat_capacity_functions 
  IMPLICIT NONE

  INTEGER ldfj
  INTEGER, PARAMETER :: mcon = 0
  INTEGER, PARAMETER :: mequa = 7
  INTEGER, PARAMETER :: nvars = 6
  REAL(8), DIMENSION(7) :: CV_temps
  REAL(8), DIMENSION(7) :: cp_difference 
  REAL(8) ::  fj(ldfj,nvars+1)
  INTEGER :: i
  INTEGER :: igo
  INTEGER :: iopt(*)
  REAL(8) :: ropt(*)
  REAL(8), SAVE, DIMENSION ( mequa ) :: t
  REAL(8) :: x(nvars)
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_high 
  REAL(8) ::  nu_mid
  INTEGER ::  N_vib, N_rot
  REAL(8), DIMENSION(7) :: diff

COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

t = CV_temps

!  For each of the number of equations (i.e. seven)
  DO i = 1, mequa
!    calculate the difference between the functions and the cp data
     diff(i) = ( Cv_harm_osc(x(1),t(i))                    &
             +   Cv_harm_osc(x(2),t(i))                    &
             +   x(3) * cv_hind_rot(x(4),t(i), x(5))         &
             +   (REAL(N_rot) -  x(3) ) * cv_hind_rot(x(6),t(i), nu_mid)  &
             - cp_difference(i) )

!   Square this difference.  This is the value of the residual function
    fj(mcon+i,nvars+1) = diff(i) * diff(i)

 END DO

!  If IGO is nonzero, compute the derivatives.
  IF ( igo /= 0 ) THEN

    DO i = 1, mequa
      
!     The Jacobian for the harm. osc. degeneracy 
      fj(mcon+i,1) = 2. * diff(i) * d_Cv_harm_osc(x(1),t(i))  
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,2) = 2. * diff(i) * d_Cv_harm_osc(x(2),t(i)) 
!     The Jacobian for the second pseudo-frequency
      fj(mcon+i,3) =  2. * diff(i) *                                   &
         ( cv_hind_rot( x(4),t(i),  x(5)) - cv_hind_rot( x(6),t(i),  nu_mid) ) 
!     The Jacobian for the hind. rot. degeneracy 
      fj(mcon+i,4) = 2. * diff(i) * x(3) *  d_cv_hind_rot( x(4),t(i), x(5)) 
!     The Jacobian for the first pseudo-barrier height
      fj(mcon+i,5) = 2. * diff(i) * x(3) *  d_cv_hind_rot_nu( x(4),t(i), x(5))
!     The Jacobian for the first pseudo-frequency
      fj(mcon+i,6) = 2. * diff(i) * d_cv_hind_rot( x(6),t(i), nu_mid) 
    END DO

  END IF

  RETURN
END SUBROUTINE Case_26_hd
!------------------------------------------------------------------------------

SUBROUTINE assign_cases(N_rot, N_vib, Total_char_freq, Total_harm_osc_freq, HR_params )

IMPLICIT NONE

 INTEGER, INTENT(IN):: N_vib
 INTEGER, INTENT(IN) :: N_rot

 INTEGER :: i

 REAL(8), DIMENSION(:), INTENT(IN) :: Total_char_freq
 REAL(8), DIMENSION(:), INTENT(INOUT) :: Total_harm_osc_freq
 REAL(8), DIMENSION(:,:), INTENT(OUT) :: HR_params

 IF (N_rot<0) THEN
     WRITE(*,*) 'ERROR!  Number of rotors is negative!'
   
  ! IF THERE ARE NO ROTORS, THEN CALL CASES 0 THROUGH 7
  ELSE IF (N_rot==0) THEN

     IF (N_vib<0) THEN
       WRITE(*,*) 'ERROR!  Number of unknown harmonic oscillator frequencies &
                    &is negative!'
        
     ELSE IF (N_Vib==0) THEN
!       WRITE(*,*) 'CASE 0'
        
     ELSE IF (N_vib==1) THEN
!       WRITE(*,*) 'CASE 1'
        CALL Case_1(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==2) THEN
!       WRITE(*,*) 'CASE 2'
        CALL Case_2(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==3) THEN
!       WRITE(*,*) 'CASE 3'
        CALL Case_3(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==4) THEN
!       WRITE(*,*) 'CASE 4'
        CALL Case_4(Total_char_freq, Total_harm_osc_freq, HR_params )
       
     ELSE IF (N_vib==5) THEN
!       WRITE(*,*) 'CASE 5'
        CALL Case_5(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==6) THEN
!       WRITE(*,*) 'CASE 6'
        CALL Case_6(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib>=7) THEN
!       WRITE(*,*) 'CASE 7'
        CALL Case_7(Total_char_freq, Total_harm_osc_freq, HR_params )
  
     ELSE
       WRITE(*,*) 'Something is wrong!'
        
     END IF

  ! IF THERE IS ONE ROTOR, THEN CALL CASES 8 THROUGH 13   
  ELSE IF (N_rot==1) THEN

     IF (N_vib<0) THEN
       WRITE(*,*) 'ERROR!  Number of unknown harmonic oscillator frequencies &
                    &is negative!'
        
     ELSE IF (N_Vib==0) THEN
!       WRITE(*,*) 'CASE 8'
        CALL Case_8(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==1) THEN
!       WRITE(*,*) 'CASE 9'
        CALL Case_9(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==2) THEN
!       WRITE(*,*) 'CASE 10'
        CALL Case_10(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==3) THEN
!       WRITE(*,*) 'CASE 11'
        CALL Case_11(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==4) THEN
!       WRITE(*,*) 'CASE 12'
        CALL Case_12(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib>=5) THEN
!       WRITE(*,*) 'CASE 13'
        CALL Case_13(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE
       WRITE(*,*) 'Something is wrong!'
        
     END IF
  
  ! IF THERE ARE TWO ROTORS, THEN CALL CASES 14 THROUGH 18   
  ELSE IF (N_rot==2) THEN

     IF (N_vib<0) THEN
       WRITE(*,*) 'ERROR!  Number of unknown harmonic oscillator frequencies &
                    &is negative!'
        
     ELSE IF (N_Vib==0) THEN
!       WRITE(*,*) 'CASE 14'
        CALL Case_14(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib==1) THEN
!       WRITE(*,*) 'CASE 15'
        CALL Case_15(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib==2) THEN
!       WRITE(*,*) 'CASE 16'
        CALL Case_16(Total_char_freq, Total_harm_osc_freq, HR_params )
       
     ELSE IF (N_vib==3) THEN
!       WRITE(*,*) 'CASE 17'
        CALL Case_17(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib>=4) THEN
!       WRITE(*,*) 'CASE 18'
        CALL Case_18(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE
       WRITE(*,*) 'Something is wrong!'
        
     END IF

  ! IF THERE ARE THREE ROTORS, THEN CALL CASES 19 THROUGH 23   
  ELSE IF (N_rot==3) THEN

     IF (N_vib<0) THEN
       WRITE(*,*) 'ERROR!  Number of unknown harmonic oscillator frequencies &
                    &is negative!'
        
     ELSE IF (N_Vib==0) THEN
!       WRITE(*,*) 'CASE 19'
        CALL Case_19(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib==1) THEN
!       WRITE(*,*) 'CASE 20'
        CALL Case_20(Total_char_freq, Total_harm_osc_freq, HR_params )
               
     ELSE IF (N_vib==2) THEN
!       WRITE(*,*) 'CASE 21'
        CALL Case_21(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib==3) THEN
!       WRITE(*,*) 'CASE 22'
        CALL Case_22(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib>=4) THEN
!       WRITE(*,*) 'CASE 23'
        CALL Case_23(Total_char_freq, Total_harm_osc_freq, HR_params )
    
     ELSE
       WRITE(*,*) 'Something is wrong!'
        
     END IF

  ! IF THERE ARE FOUR OR MORE ROTORS, THEN CALL CASES 24,25,26,22,23      
  ELSE IF (N_rot>=4) THEN

     IF (N_vib<0) THEN
       WRITE(*,*) 'ERROR!  Number of unknown harmonic oscillator frequencies &
                    &is negative!'
        
     ELSE IF (N_Vib==0) THEN
!       WRITE(*,*) 'CASE 24'
        CALL Case_24(Total_char_freq, Total_harm_osc_freq, HR_params )

     ELSE IF (N_vib==1) THEN
!       WRITE(*,*) 'CASE 25'
        CALL Case_25(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==2) THEN
!       WRITE(*,*) 'CASE 26'
        CALL Case_26(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib==3) THEN
!       WRITE(*,*) 'CASE 22'
        CALL Case_22(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE IF (N_vib>=4) THEN
!       WRITE(*,*) 'CASE 23'
        CALL Case_23(Total_char_freq, Total_harm_osc_freq, HR_params )
        
     ELSE
       WRITE(*,*) 'Something is wrong!'
        
     END IF

  ELSE 
     WRITE(*,*) 'Something is wrong!'
  END IF





END SUBROUTINE assign_cases


END MODULE Cases
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
MODULE Open_stuff

CONTAINS

SUBROUTINE read_bonds(data, bond_info, bond_degeneracy)

  INTEGER :: RsCH3, RdCH2, CtCH, RSCH2sR, CdCHsR, Aldehyde,Cumulene
  INTEGER :: Ketene, CtCsR, RsCHsR2, CdCsR2, Ketone, RsCsR3
  INTEGER :: RsCH2r, RdCHr, RsCHrsR, CdCrsR, OdCrsR, RsCrsR2
  INTEGER :: Alcohol, Ether, ROOH, ROOR, Peroxy
  INTEGER :: Rings
  
  INTEGER :: N_atom, N_rot, linearity
  INTEGER, INTENT(IN), DIMENSION(:) :: data
  INTEGER, INTENT(OUT), DIMENSION(:) :: bond_info
  INTEGER, INTENT(OUT) :: bond_degeneracy

  N_atom = data(1)
  N_rot =  data(2)
  linearity =  data(3)

! non radicals
! terminal groups
  RsCH3     = data(4)
  RdCH2     = data(5)
  CtCH      = data(6)

! non-terminal groups (2 heavy atoms)
  RsCH2sR   = data(7)
  CdCHsR    = data(8)
  Aldehyde  = data(9)
  Cumulene  = data(10)
  Ketene    = data(11)
  CtCsR     = data(12)
  
  ! non-terminal groups (3 heavy atoms)
  RsCHsR2   = data(13)
  CdCsR2    = data(14)
  Ketone    = data(15)

  ! non-terminal groups (4 heavy atoms)
  RsCsR3    = data(16)

! Single radicals
! terminal groups
  RsCH2r    = data(17)
  RdCHr     = data(18)

! non-terminal groups
  RsCHrsR   = data(19)
  CdCrsR    = data(20)
  OdCrsR    = data(21)
  RsCrsR2   = data(22)

! oxygen-specific groups
  Alcohol   = data(23)
  Ether     = data(24)
  ROOH      = data(25)
  ROOR      = data(26)
  Peroxy    = data(27)

! cyclic species C-H stretch
  Rings     = data(28)

  ! RMG will probably double-count the number of ROOR bonds.
  ! Assuming that is the case, divide by 2:
  IF (ROOR/2==0) THEN
     !it is an odd number, so do nothing
  ELSE
     ROOR = ROOR / 2
  ENDIF

Bond_degeneracy =  8*RsCH3 + 5*RdCH2 + 3*CtCH + 7*RSCH2sR + 5*CdCHsR + &
                   5*Aldehyde + 3*Cumulene + 3*Ketene + 2*CtCsR + & 
                   6*RsCHsR2 + 4*CdCsR2 + 4*Ketone + 5*RsCsR3 + 5*RsCH2r + & 
                   3*RdCHr + 4*RsCHrsR + 2*CdCrsR + 2*OdCrsR + 3*RsCrsR2 + & 
                   2*Alcohol + 1*Ether + 4*ROOH + 2*ROOR + 2*Peroxy + &
                   2*Rings

WRITE(*,*) 'degeneracy = ', Bond_degeneracy
  

  bond_info(1)  = RsCH3
  bond_info(2)  = RdCH2
  bond_info(3)  = CtCH
  bond_info(4)  = RsCH2sR
  bond_info(5)  = CdCHsR
  bond_info(6)  = Aldehyde
  bond_info(7)  = Cumulene
  bond_info(8)  = Ketene
  bond_info(9)  = CtCsR
  bond_info(10) = RsCHsR2
  bond_info(11) = CdCsR2
  bond_info(12) = Ketone
  bond_info(13) = RsCsR3
  bond_info(14) = RsCH2r
  bond_info(15) = RdCHr
  bond_info(16) = RsCHrsR
  bond_info(17) = CdCrsR
  bond_info(18) = OdCrsR
  bond_info(19) = RsCrsR2
  bond_info(20) = Alcohol
  bond_info(21) = Ether
  bond_info(22) = ROOH
  bond_info(23) = ROOR 
  bond_info(24) = Peroxy
  bond_info(25) = Rings


END SUBROUTINE read_bonds
!-----------------------------------------------------------------------------

SUBROUTINE Calculate_RRHO_HR_params(input_file, output_file )

  USE heat_capacity_functions
  USE Cases
  USE frequencies

IMPLICIT NONE

! Global Variables
  REAL(8), DIMENSION(7) :: CV_temps 
  REAL(8) :: R 

! These variables are specific to the molecule and are used by the solver
  INTEGER :: N_atoms
  INTEGER :: linearity
  INTEGER :: degeneracy
  INTEGER :: N_vib
  INTEGER :: N_rot
  INTEGER, DIMENSION(28) :: data
  INTEGER, DIMENSION(25) :: bond_info ! must be data - 3 
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Total_predicted_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: HR_params
  REAL(8), DIMENSION(7) :: cp_data
  REAL(8), DIMENSION(7) :: sum_cp_bonds
  REAL(8), DIMENSION(7) :: cp_difference
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Cv
  INTEGER :: i, j
  REAL(8) :: trans_rot_cp_to_cv

  ! These three variables are the "fixed" frequencies for use in the 
  ! Hindered rotor fitting.
  REAL(8) ::  nu_low 
  REAL(8) ::  nu_mid 
  REAL(8) ::  nu_high

 ! These variables are used to open the file
  INTEGER :: OpenStatus
  CHARACTER(20) :: input_file
  CHARACTER(20) :: output_file


  ! these variables are set to the common block.
  ! They are required by the fitting routine within DQED
  COMMON N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

  !  Set the common block varibles:
  CV_temps = (/ 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0 /)
  nu_low = 50.0
  nu_mid = 150.0
  nu_high = 300.0

  ! Open in the input file, read everything, close it back up.
  open (UNIT = 12, FILE = input_file, STATUS = 'OLD', & 
       ACTION = 'READ', IOSTAT = OpenStatus)
  READ (12,*) cp_data, data
  CLOSE(12)
  
  ! call the subroutine which separates the data info and calculates
  ! how many RRHO frequencies it can determine from the structure
  call read_bonds(data, bond_info,  degeneracy)
  
  ! Normalize the cp data by the ideal gas constant
  R = 8.3145 ! Gas Constant in J/mol-K
  R = 1.9859 ! Gas Constant in cal/mol-K
  cp_data = cp_data/R
  
  ! Specify the number of atoms, hind. rotors, unknown harm. osc., & linearity
  N_atoms   = data(1)
  N_rot     = data(2)
  linearity = data(3)
  
  ! Depending upon the linearity of the molecule:
  ! Detemine whether to subtract 3.5 or 4 from the supplied heat capacity
  ! data;  the justification for 3.5 vs 4 is described below in the comments
  ! below the subroutine calc_predicted_freq.
  ! Determine whether the total number of vibrational models is 3N-5 or 3N-6
  IF (linearity<0) THEN
     WRITE(*,*) 'ERROR!  Linearity is less than zero!'
  ELSE IF (linearity==0) THEN
     WRITE(*,*) 'Linear molecule'
     trans_rot_cp_to_cv = 3.5
     N_vib = 3 * (N_atoms) - 5 - (N_rot) - (degeneracy)
     IF (N_vib < 0) THEN
!       WRITE(*,*) 'Something is wrong.  The system is over specified.'
!       WRITE(*,*) 'Either reduce the number of rotors or reduce the number of bond types.'
!       WRITE(*,*) 'Program will now exit to avoid crashing.'
        STOP
     ENDIF
     ALLOCATE( Total_harm_osc_freq(3 * (N_atoms) - 5 - (N_rot) ) )
  ELSE IF (linearity==1) THEN
     WRITE(*,*) 'Nonlinear molecule'
     trans_rot_cp_to_cv = 4.0
     N_vib = 3 * (N_atoms) - 6 - (N_rot) - (degeneracy)
	 IF (N_vib < 0) THEN
!       WRITE(*,*) 'Something is wrong.  The system is over specified.'
!       WRITE(*,*) 'Either reduce the number of rotors or reduce the number of bond types.'
!       WRITE(*,*) 'Program will now exit to avoid crashing.'
        STOP
     ENDIF
     ALLOCATE( Total_harm_osc_freq(3 * (N_atoms) - 6 - (N_rot) ) )
  ELSE
     WRITE(*,*) 'ERROR!  Linearity greater than one!'
  END IF
  

  ! Allocate the total number of characteristic frequencies
  ALLOCATE(Total_predicted_freq( degeneracy ) )
  
  ! Call the program which calculates the total number of characteristic freq.s
  CALL calc_predicted_freq(bond_info, degeneracy, Total_predicted_freq )
  
  ! Calculate the heat capacity at each CV_temps due to the char. freq.'s
  ! To isolate the contribution towards the heat capacity of the unknown 
  ! Harmonic oscillator frequencies and hindered rotors, we must subtract
  ! The translational, rotational, and known harm. osc. heat capacity from
  ! The supplied cp data
  ! cp_differences = 
  ! Cp - 3/2 (for translation)
  !    - 3/2 (for rotation) unless it is linear, in which case -1
  !    - 1 (for Cp - Cv = R)
  !    - sum_cp_bonds (the contribution of the characteristic freq.'s
  ALLOCATE(Cv( degeneracy ) )
  ! For each temperature, cycle through the i's
  DO i = 1, size(CV_temps)
  !  For each frequency, cycle through the j's
     DO j = 1, size(Total_predicted_freq)
        
        ! calculate the heat capacity for that frequency at that temperature
        CALL cv_vib( Total_predicted_freq(j), CV_temps(i), Cv(j))
        
     END DO
     ! sum over all the frequencies at that temperature
     sum_cp_bonds(i) = SUM(Cv)
     ! subtract everything.  What is left is the contribution from the 
     ! unknown active degrees of freedom.
     cp_difference(i) = cp_data(i) - sum_cp_bonds(i) - trans_rot_cp_to_cv 
  END DO
  DEALLOCATE(Cv)
  
  ! Set the first 3N-6 (or 3N-5) total harmonic oscillator frequencies equal to
  ! The frequencies calculated by the subroutine "calc_total_freq"
  DO i = 1,size(Total_predicted_freq)
     Total_harm_osc_freq(i) = Total_predicted_freq(i)
  ENDDO

  ! BEGIN THE SECTION WHICH DETERMINES WHICH SOLVER TO CALL

  ! Allocate the number of hindered rotor parameters
  ALLOCATE( HR_params(N_rot,2) )
  HR_params = 0.0
  
  ! Call the subroutine which, based upon N_rot and N_vib, solves the unknonwn
  ! cp data and returns the total harmonic osc. frequencies and the hindered
  ! rotor density of state.
  CALL assign_cases(N_rot, N_vib, Total_predicted_freq, Total_harm_osc_freq, HR_params )


  ! WRITE EVERYTHING TO THE OUTPUT FILE
  ! open the file
  open (UNIT = 22, FILE = output_file)

  ! Write the number of atoms, internal rotors, and the linearity
  WRITE(22,*) N_atoms
  WRITE(22,*) N_rot
  WRITE(22,*) linearity
  
  ! Write the characteristic RRHO frequencies determined by the structure
  WRITE(22,*) ''
  DO i = 1,size(Total_predicted_freq)
     WRITE(22,*) Total_predicted_freq(i)
  ENDDO
  WRITE(22,*) ''

  ! Write the additional RRHO frequencies determined by the heat capacity
  WRITE(22,*) ''
  DO i = (size(Total_predicted_freq)+1), size(Total_harm_osc_freq)
     WRITE(22,*) Total_harm_osc_freq(i)
  ENDDO
  WRITE(22,*) ''

  ! Write the hindered rotor frequencies and barrier heights
  DO i = 1,N_rot
     WRITE(22,*) HR_params(i,:)
  ENDDO
  WRITE(22,*) ''
  
  ! Close the output file
  CLOSE(22)

  ! Deallocate what's left open
  DEALLOCATE( Total_harm_osc_freq )
  DEALLOCATE(Total_predicted_freq )
  DEALLOCATE( HR_params )

END SUBROUTINE Calculate_RRHO_HR_params


END MODULE Open_stuff
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

PROGRAM main

  USE heat_capacity_functions
  USE Cases
  USE Open_stuff
  USE frequencies

  IMPLICIT NONE

  INTEGER :: i

  ! These variables are for the input and output filenames
  CHARACTER(20) :: input_file
  CHARACTER(20) :: output_file
  
  ! Fluffy stuff  FEEL FREE TO REMOVE
  CHARACTER(10) :: time
  CHARACTER(8) :: date
  CHARACTER(5) :: zone
  INTEGER, DIMENSION(8) :: start_value, end_value


! You can change this so that it will read the names from the command like.
! It is probably more useful that way.
!  READ(*,*) input_file, output_file
  input_file = 'dat'
  output_file = 'rho_input'

! Feel free to cut this.  I use it to determine how long the code runs.  
  CALL  date_and_time(date, time, zone, start_value)
  
! This is the command that calls the main program.
  CALL Calculate_RRHO_HR_params(input_file, output_file)

! Ditto for the next five lines.  Cut them.
  CALL  date_and_time(date, time, zone, end_value)
  WRITE(*,*) ''
  WRITE(*,70) (end_value(5) - start_value(5)),  (end_value(6) - start_value(6)), &
       (end_value(7) - start_value(7)),  (end_value(8) - start_value(8))
70 FORMAT('Total computation time is: ',I2.2,' hrs, ',I2.2, ' min, ',I2.2, ' s, ',I4.4 ,' ms')  

  STOP

END PROGRAM main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fitSpectralDataToHeatCapacity(Cv, Tlist, Ntemp, Nvib, Nhind, vib, hind)

	use cases

	integer, intent(in) :: Ntemp
	integer, intent(in) :: Nvib
	integer, intent(in) :: Nhind
	real(8), dimension(1:Ntemp), intent(in) :: Cv
	real(8), dimension(1:Ntemp), intent(in) :: Tlist
	real(8), dimension(1:Nvib), intent(out) :: vib
	real(8), dimension(1:Nhind,1:2), intent(out) :: hind

	real(8), dimension(1:0) :: totalCharFreq

	integer :: N_vib, N_rot
	real(8) :: nu_low, nu_mid, nu_high
	real(8), dimension(7) :: cp_difference, CV_temps

	! These variables are set to the common block
	! They are required by the fitting routine within DQED
	common N_vib, N_rot, nu_low, nu_high, nu_mid, cp_difference, CV_temps

	! Initialize output arrays to zero
	vib = 0.0
	hind = 0.0

	! Set common variables
	N_vib = Nvib
	N_rot = Nhind
	nu_low = 100.0
	nu_mid = 150.0
	nu_high = 300.0
	cp_difference = Cv
	CV_temps = Tlist

	! Fit the harmonic oscillator and hindered rotor modes to the heat capacity
	! The function assign_cases is in calc_freq_code.f90
	call assign_cases(Nhind, Nvib, totalCharFreq, vib, hind)

end subroutine
