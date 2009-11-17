subroutine difcen ( fj, func, fx, iopt, ldfj, mcon, mequa, nvars, &
  ropt, x )

!***********************************************************************
!
!! DIFCEN estimates a jacobian using central differences.
!
!  Modified:
!
!    18 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJ(LDFJ,NVARS), the estimated jacobian.
!
!    Input, external FUNC, the name of the user written
!    function evaluation routine.  FUNC should have the form:
!
!      subroutine func ( fx, iopt, mcon, mequa, nvars, ropt, x )
!
!    and should accept X as input, and return in FX the value
!    of the MEQUA+MCON functions.
!
!    Workspace, real ( kind = 8 ) FX(MEQUA+MCON).
!
!    Throughput, integer IOPT(*), parameters to be passed to FUNC.
!
!    Input, integer LDFJ, the leading dimension of FJ, which must
!    be at least MEQUA+MCON.
!
!    Input, integer MCON, the number of constraints.
!
!    Input, integer MEQUA, the number of nonlinear functions.
!
!    Input, integer NVARS, the number of variables.
!
!    Throughput, real ( kind = 8 ) ROPT(*), parameters to be passed to FUNC.
!
!    Input, real ( kind = 8 ) X(NVARS), the point at which the
!    jacobian should be evaluated.
!
  implicit none

  integer ldfj
  integer mcon
  integer mequa
  integer nvars

  real ( kind = 8 ) dxj
  real ( kind = 8 ) eps
  real ( kind = 8 ) fj(ldfj,nvars)
  external func
  real ( kind = 8 ) fx(mequa+mcon)
  integer iopt(*)
  integer j
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) x(nvars)
  real ( kind = 8 ) xsave
!
!  Get the square root of the machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  Consider each component X(J) of the set of variables.
!
  do j = 1, nvars
!
!  Set the appropriate increment DXJ to X(J).
!
    dxj = eps * ( abs ( x(j) ) + 1.0D+00 )
!
!  Make a copy XP of X, with X(J) incremented by DXJ.
!
    xsave = x(j)

    x(j) = xsave + dxj
!
!  Evaluate F(XP).
!
    call func ( fx, iopt, mcon, mequa, nvars, ropt, x )
!
!  Save F(XP).
!
    fj(1:mequa+mcon,j) = fx(1:mequa+mcon)
!
!  Make a copy XM of X, with X(J) decremented by DXJ.
!
    x(j) = xsave - dxj
!
!  Evaluate F(XM).
!
    call func ( fx, iopt, mcon, mequa, nvars, ropt, x )
!
!  Estimate the partial derivative d F/d X(J) by (F(XP)-F(XM))/(2*DXJ)
!
    fj(1:mequa+mcon,j) = ( fj(1:mequa+mcon,j) - fx(1:mequa+mcon) ) &
      / ( 2.0D+00 * dxj )
!
!  Restore the value of X(J).
!
    x(j) = xsave

  end do

  return
end
subroutine diffor ( fj, func, fx, iopt, ldfj, mcon, mequa, nvars, &
  ropt, x )

!***********************************************************************
!
!! DIFFOR estimates a jacobian using forward differences.
!
!  Modified:
!
!    18 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJ(LDFJ,NVARS), the estimated jacobian.
!
!    Input, external FUNC, the name of the user written
!    function evaluation routine.  FUNC should have the form:
!
!      subroutine func ( fx, iopt, mcon, mequa, nvars, ropt, x )
!
!    and should accept X as input, and return in FX the value
!    of the MEQUA+MCON functions.
!
!    Workspace, real ( kind = 8 ) FX(MEQUA+MCON).
!
!    Throughput, integer IOPT(*), parameters to be passed
!    to FUNC.
!
!    Input, integer LDFJ, the leading dimension of FJ, which must
!    be at least MEQUA+MCON.
!
!    Input, integer MCON, the number of constraints.
!
!    Input, integer MEQUA, the number of nonlinear functions.
!
!    Input, integer NVARS, the number of variables.
!
!    Throughput, real ( kind = 8 ) ROPT(*), parameters to be passed to FUNC.
!
!    Input, real ( kind = 8 ) X(NVARS), the point at which the
!    jacobian should be evaluated.
!
  implicit none

  integer ldfj
  integer mcon
  integer mequa
  integer nvars

  real ( kind = 8 ) dxj
  real ( kind = 8 ) eps
  real ( kind = 8 ) fj(ldfj,nvars)
  external func
  real ( kind = 8 ) fx(mequa+mcon)
  integer iopt(*)
  integer j
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) x(nvars)
  real ( kind = 8 ) xsave
!
!  Evaluate F(X) and save it in FX.
!
  call func ( fx, iopt, mcon, mequa, nvars, ropt, x )
!
!  Get the square root of the machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  Consider each component X(J) of the set of variables.
!
  do j = 1, nvars
!
!  Set the appropriate increment DXJ to X(J).
!
    dxj = eps * ( abs ( x(j) ) + 1.0D+00 )
!
!  Make a copy XP of X, with X(J) incremented by DXJ.
!
    xsave = x(j)

    x(j) = xsave + dxj
!
!  Evaluate F(XP) and store it in column J.
!
    call func ( fj(1,j), iopt, mcon, mequa, nvars, ropt, x )
!
!  Estimate the partial derivative d F/d X(J) by (F(XP)-F(X))/DXJ
!
    fj(1:mequa+mcon,j) = ( fj(1:mequa+mcon,j) - fx(1:mequa+mcon) ) / dxj
!
!  Restore the value of X(J).
!
    x(j) = xsave

  end do

  return
end
function idamax ( n, x, incx )

!*******************************************************************************
!
!! IDAMAX finds the index of the vector element of maximum absolute value.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of SX.
!
!    Output, integer IDAMAX, the index of the element of SX of maximum
!    absolute value.
!
  implicit none

  integer i
  integer incx
  integer idamax
  integer ix
  integer n
  real ( kind = 8 ) damax
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

    idamax = 0

  else if ( n == 1 ) then

    idamax = 1

  else if ( incx == 1 ) then

    idamax = 1
    damax = abs ( x(1) )

    do i = 2, n

      if ( abs ( x(i) ) > damax ) then
        idamax = i
        damax = abs ( x(i) )
      end if

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    idamax = 1
    damax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > damax ) then
        idamax = i
        damax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine ivout ( n, ix, title, idigit )

!***********************************************************************
!
!! IVOUT prints integer vectors.
!
!  Author:
!
!    John Wisniewski and Richard Hanson,
!    SANDIA LABS ALBUQUERQUE.
!
!  Parameters:
!
!    Input, integer N, the size of the vector IX.
!
!    Input, integer IX(N), the array to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
!    Input, integer IDIGIT, indicates the number of digits to print.
!    PRINT UP TO IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
!    THE SUBPROGRAM WILL CHOOSE THAT integer 4,6,10 OR 14
!    WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF
!    PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
!    TO WRITE EACH LINE OF OUTPUT OF THE ARRAY IX(*). (THIS
!    CAN BE USED ON MOST TIME-SHARING TERMINALS). IF
!    IDIGIT.GE.0, 133 PRINTING COLUMNS ARE UTILIZED. (THIS CAN
!    BE USED ON MOST LINE PRINTERS).
!
  implicit none

  integer n

  integer idigit
  integer ix(n)
  integer k1
  integer k2
  integer ndigit
  character ( len = * ) title

  write ( *, '(a)' ) trim ( title )

  if ( n <= 0 ) then
    return
  end if

  ndigit = idigit

  if ( idigit == 0 ) then
    ndigit = 4
  end if

  if ( idigit >= 0 ) go to 80

  ndigit = -idigit

  if ( ndigit <= 4 )then

    do k1 = 1, n, 10
      k2 = min ( n, k1+9 )
      write(*,1000) k1, k2, ix(k1:k2)
    end do

    return

  end if

  if ( ndigit > 6) go to 40

  do k1=1,n,7
    k2 = min(n,k1+6)
    write(*,1001) k1,k2, ix(k1:k2)
  end do

  return

   40 continue
  if ( ndigit > 10) go to 60

  do k1=1,n,5
    k2=min(n,k1+4)
    write(*,1002) k1,k2, ix(k1:k2)
  end do

  return

   60 continue

  do k1=1,n,3
    k2 = min(n,k1+2)
    write(*,1003) k1,k2, ix(k1:k2)
  end do

  return

   80 continue

  if ( ndigit > 4) go to 100

  do k1=1,n,20
    k2 = min(n,k1+19)
    write(*,1000) k1,k2, ix(k1:k2)
  end do

  return

  100 continue

  if ( ndigit > 6) go to 120

  do k1=1,n,15
    k2 = min(n,k1+14)
    write(*,1001) k1,k2, ix(k1:k2)
  end do

  return

  120 continue

  if ( ndigit > 10) go to 140

  do k1=1,n,10
    k2 = min(n,k1+9)
    write(*,1002) k1,k2, ix(k1:k2)
  end do

  return

  140 continue

  do k1=1,n,7
    k2 = min(n,k1+6)
    write(*,1003) k1,k2, ix(k1:k2)
  end do

  return
 1000 format(1x,i4,' - ',i4,20(1x,i5))
 1001 format(1x,i4,' - ',i4,15(1x,i7))
 1002 format(1x,i4,' - ',i4,10(1x,i11))
 1003 format(1x,i4,' - ',i4,7(1x,i15))
end
function damax ( n, x, incx )

!*******************************************************************************
!
!! DAMAX returns the maximum absolute value of the entries in a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DAMAX, the maximum absolute value of an 
!    element of X.
!
  implicit none

  integer i
  integer incx
  integer ix
  integer n
  real ( kind = 8 ) damax
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

    damax = 0.0D+00

  else if ( n == 1 ) then

    damax = abs ( x(1) )

  else if ( incx == 1 ) then

    damax = abs ( x(1) )

    do i = 2, n
      if ( abs ( x(i) ) > damax ) then
        damax = abs ( x(i) )
      end if
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    damax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > damax ) then
        damax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function dasum ( n, x, incx )

!*******************************************************************************
!
!! DASUM sums the absolute values of the entries of a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DASUM, the sum of the absolute values of 
!    the vector.
!
  implicit none

  integer i
  integer incx
  integer ix
  integer m
  integer n
  real ( kind = 8 ) dasum
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)

  stemp = 0.0D+00

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 6 )

    do i = 1, m
      stemp = stemp + abs ( x(i) )
    end do

    do i = m+1, n, 6
      stemp = stemp + abs ( x(i)   ) + abs ( x(i+1) ) + abs ( x(i+2) ) &
                    + abs ( x(i+3) ) + abs ( x(i+4) ) + abs ( x(i+5) )
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      stemp = stemp + abs ( x(ix) )
      ix = ix + incx
    end do

  end if

  dasum = stemp

  return
end
subroutine daxpy ( n, sa, x, incx, y, incy )

!*******************************************************************************
!
!! DAXPY adds a constant times one vector to another.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input, real ( kind = 8 ) X(*), the vector to be scaled and added to Y.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), the vector to which a multiple 
!    of X is to be added.
!
!    Input, integer INCY, the increment between successive entries of Y.
!
  implicit none

  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( sa == 0.0D+00 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = y(1:n) + sa * x(1:n)

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = y(iy) + sa * x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine dbocls ( w, mdw, mcon, mrows, ncols, bl, bu, ind, iopt, x, &
  rnormc, rnorm, mode, rw, iw )

!***********************************************************************
!
!! DBOCLS solves a bounded and constrained least squares problem.
!
!  Discussion:
!
!    DBOCLS solves the bounded and constrained least squares
!    problem consisting of solving the equation
!      E*X = F  (in the least squares sense)
!    subject to the linear constraints
!      C*X = Y.
!
!    This subprogram solves the bounded and constrained least squares
!    problem. The problem statement is:
!
!    Solve E*X = F (least squares sense), subject to constraints
!    C*X=Y.
!
!    In this formulation both X and Y are unknowns, and both may
!    have bounds on any of their components.  This formulation
!    of the problem allows the user to have equality and inequality
!    constraints as well as simple bounds on the solution components.
!
!    This constrained linear least squares subprogram solves E*X=F
!    subject to C*X=Y, where E is MROWS by NCOLS, C is MCON by NCOLS.
!
!    The user must have dimension statements of the form
!
!      DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON), BU(NCOLS+MCON),
!     * X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON)
!       integer IND(NCOLS+MCON), IOPT(17+NI), IW(2*(NCOLS+MCON))
!
!    (here NX=number of extra locations required for the options; NX=0
!    if no options are in use. Also NI=number of extra locations
!    for options 1-9.
!
!  Author:
!
!    Richard Hanson, Sandia National Laboratory
!
!  Reference:
!
!    Richard Hanson,
!    Linear Least Squares with Bounds and Linear Constraints,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 7, Number 3, July 1986.
!
!    INPUT
!    -----
!
!    -------------------------
!    W(MDW,*),MCON,MROWS,NCOLS
!    -------------------------
!     The array W contains the (possibly null) matrix [C:*] followed by
!     [E:F].  This must be placed in W as follows:
!          [C  :  *]
!     W  = [       ]
!          [E  :  F]
!     The (*) after C indicates that this data can be undefined. The
!     matrix [E:F] has MROWS rows and NCOLS+1 columns. The matrix C is
!     placed in the first MCON rows of W(*,*) while [E:F]
!     follows in rows MCON+1 through MCON+MROWS of W(*,*). The vector F
!     is placed in rows MCON+1 through MCON+MROWS, column NCOLS+1. The
!     values of MDW and NCOLS must be positive; the value of MCON must
!     be nonnegative. An exception to this occurs when using option 1
!     for accumulation of blocks of equations. In that case MROWS is an
!     OUTPUT variable only, and the matrix data for [E:F] is placed in
!     W(*,*), one block of rows at a time. See IOPT(*) contents, option
!     number 1, for further details. The row dimension, MDW, of the
!     array W(*,*) must satisfy the inequality:
!
!     If using option 1,
!                     MDW >= MCON + max(max. number of
!                     rows accumulated, NCOLS)
!
!     If using option 8, MDW >= MCON + MROWS.
!     Else, MDW >= MCON + max(MROWS, NCOLS).
!
!     Other values are errors, but this is checked only when using
!     option=2.  The value of MROWS is an output parameter when
!     using option number 1 for accumulating large blocks of least
!     squares equations before solving the problem.
!     See IOPT(*) contents for details about option 1.
!
!    ------------------
!    BL(*),BU(*),IND(*)
!    ------------------
!     These arrays contain the information about the bounds that the
!     solution values are to satisfy. The value of IND(J) tells the
!     type of bound and BL(J) and BU(J) give the explicit values for
!     the respective upper and lower bounds on the unknowns X and Y.
!     The first NVARS entries of IND(*), BL(*) and BU(*) specify
!     bounds on X; the next MCON entries specify bounds on Y.
!
!    1.    For IND(J)=1, require X(J) >= BL(J);
!          IF J > NCOLS,        Y(J-NCOLS) >= BL(J).
!          (the value of BU(J) is not used.)
!    2.    For IND(J)=2, require X(J) <= BU(J);
!          IF J > NCOLS,        Y(J-NCOLS) <= BU(J).
!          (the value of BL(J) is not used.)
!    3.    For IND(J)=3, require X(J) >= BL(J) and
!                                X(J) <= BU(J);
!          IF J > NCOLS,        Y(J-NCOLS) >= BL(J) and
!                                Y(J-NCOLS) <= BU(J).
!          (to impose equality constraints have BL(J)=BU(J)=
!          constraining value.)
!    4.    For IND(J)=4, no bounds on X(J) or Y(J-NCOLS) are required.
!          (the values of BL(J) and BU(J) are not used.)
!
!     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
!     IND(J)=3 (upper and lower bounds) the condition BL(J)  >  BU(J)
!     is  an  error.   The values BL(J), BU(J), J  >  NCOLS, will be
!     changed.  Significant changes mean that the constraints are
!     infeasible.  (Users must make this decision themselves.)
!     The new values for BL(J), BU(J), J  >  NCOLS, define a
!     region such that the perturbed problem is feasible.  If users
!     know that their problem is feasible, this step can be skipped
!     by using option number 8 described below.
!
!    -------
!    IOPT(*)
!    -------
!     This is the array where the user can specify nonstandard options
!     for DBOCLS( ). Most of the time this feature can be ignored by
!     setting the input value IOPT(1)=99. Occasionally users may have
!     needs that require use of the following subprogram options. For
!     details about how to use the options see below: IOPT(*) CONTENTS.
!
!     Option Number   Brief Statement of Purpose
!     ------ ------   ----- --------- -- -------
!           1         Return to user for accumulation of blocks
!                     of least squares equations.  The values
!                     of IOPT(*) are changed with this option.
!                     The changes are updates to pointers for
!                     placing the rows of equations into position
!                     for processing.
!           2         Check lengths of all arrays used in the
!                     subprogram.
!           3         Column scaling of the data matrix, [C].
!                                                        [E]
!           4         User provides column scaling for matrix [C].
!                                                             [E]
!           5         Provide option array to the low-level
!                     subprogram DBOLS( ).
!                     {Provide option array to the low-level
!                     subprogram DBOLSM( ) by imbedding an
!                     option array within the option array to
!                     DBOLS(). Option 6 is now disabled.}
!           7         Move the IOPT(*) processing pointer.
!           8         Do not preprocess the constraints to
!                     resolve infeasibilities.
!           9         Do not pretriangularize the least squares matrix.
!          99         No more options to change.
!
!    ----
!    X(*)
!    ----
!     This array is used to pass data associated with options 4,5 and
!     6. Ignore this parameter (on input) if no options are used.
!     Otherwise see below: IOPT(*) CONTENTS.
!
!
!    OUTPUT
!    ------
!
!    -----------------
!    X(*),RNORMC,RNORM
!    -----------------
!     The array X(*) contains a solution (if MODE >=0 or  == -22) for
!     the constrained least squares problem. The value RNORMC is the
!     minimum residual vector length for the constraints C*X - Y = 0.
!     The value RNORM is the minimum residual vector length for the
!     least squares equations. Normally RNORMC=0, but in the case of
!     inconsistent constraints this value will be nonzero.
!     The values of X are returned in the first NVARS entries of X(*).
!     The values of Y are returned in the last MCON entries of X(*).
!
!    ----
!    MODE
!    ----
!     The sign of MODE determines whether the subprogram has completed
!     normally, or encountered an error condition or abnormal status. A
!     value of MODE >= 0 signifies that the subprogram has completed
!     normally. The value of mode (>= 0) is the number of variables
!     in an active status: not at a bound nor at the value zero, for
!     the case of free variables. A negative value of MODE will be one
!     of the cases (-57)-(-41), (-37)-(-22), (-19)-(-2). Values  <  -1
!     correspond to an abnormal completion of the subprogram. These
!     error messages are in groups for the subprograms DBOCLS(),
!     DBOLSM(), and DBOLS().  An approximate solution will be returned
!     to the user only when max. iterations is reached, MODE=-22.
!
!    -----------
!    RW(*),IW(*)
!    -----------
!     These are working arrays.  (normally the user can ignore the
!     contents of these arrays.)
!
!    IOPT(*) CONTENTS
!    ------- --------
!     The option array allows a user to modify some internal variables
!     in the subprogram without recompiling the source code. A central
!     goal of the initial software design was to do a good job for most
!     people. Thus the use of options will be restricted to a select
!     group of users. The processing of the option array proceeds as
!     follows: a pointer, here called LP, is initially set to the value
!     1. At the pointer position the option number is extracted and
!     used for locating other information that allows for options to be
!     changed. The portion of the array IOPT(*) that is used for each
!     option is fixed; the user and the subprogram both know how many
!     locations are needed for each option. The value of LP is updated
!     for each option based on the amount of storage in IOPT(*) that is
!     required. A great deal of error checking is done by the
!     subprogram on the contents of the option array. Nevertheless it
!     is still possible to give the subprogram optional input that is
!     meaningless. For example option 4 uses the locations
!     X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing scaling data.
!     The user must manage the allocation of these locations.
!
!   1
!   -
!     This option allows the user to solve problems with a large number
!     of rows compared to the number of variables. The idea is that the
!     subprogram returns to the user (perhaps many times) and receives
!     new least squares equations from the calling program unit.
!     Eventually the user signals "that's all" and a solution is then
!     computed. The value of MROWS is an output variable when this
!     option is used. Its value is always in the range 0 <= MROWS
!     <= NCOLS+1. It is the number of rows after the
!     triangularization of the entire set of equations. If LP is the
!     processing pointer for IOPT(*), the usage for the sequential
!     processing of blocks of equations is
!
!
!        IOPT(LP)=1
!         Move block of equations to W(*,*) starting at
!         the first row of W(*,*).
!        IOPT(LP+3)=# of rows in the block; user defined
!
!     The user now calls DBOCLS( ) in a loop. The value of IOPT(LP+1)
!     directs the user's action. The value of IOPT(LP+2) points to
!     where the subsequent rows are to be placed in W(*,*). Both of
!     these values are first defined in the subprogram. The user
!     changes the value of IOPT(LP+1) (to 2) as a signal that all of
!     the rows have been processed.
!
!
!      .<LOOP
!      . CALL DBOCLS( )
!      . IF(IOPT(LP+1) .EQ. 1) THEN
!      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
!      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
!      .    W(*,*) STARTING AT ROW MCON + IOPT(LP+2).
!      .
!      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN
!      .       IOPT(LP+1)=2
!      .<------CYCLE LOOP
!      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN
!      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0
!      . ELSE
!      . ERROR CONDITION; SHOULD NOT HAPPEN.
!      .<END LOOP
!
!     Use of this option adds 4 to the required length of IOPT(*).
!
!   2
!   -
!     This option is useful for checking the lengths of all arrays used
!     by DBOCLS( ) against their actual requirements for this problem.
!     The idea is simple: the user's program unit passes the declared
!     dimension information of the arrays. These values are compared
!     against the problem-dependent needs within the subprogram. If any
!     of the dimensions are too small an error message is printed and a
!     negative value of MODE is returned, -41 to -47. The printed error
!     message tells how long the dimension should be. If LP is the
!     processing pointer for IOPT(*),
!
!        IOPT(LP)=2
!        IOPT(LP+1)=Row dimension of W(*,*)
!        IOPT(LP+2)=Col. dimension of W(*,*)
!        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*)
!        IOPT(LP+4)=Dimension of X(*)
!        IOPT(LP+5)=Dimension of RW(*)
!        IOPT(LP+6)=Dimension of IW(*)
!        IOPT(LP+7)=Dimension of IOPT(*)
!         .
!        CALL DBOCLS( )
!
!     Use of this option adds 8 to the required length of IOPT(*).
!
!   3
!   -
!     This option can change the type of scaling for the data matrix.
!     Nominally each nonzero column of the matrix is scaled so that the
!     magnitude of its largest entry is equal to the value ONE. If LP
!     is the processing pointer for IOPT(*),
!
!        IOPT(LP)=3
!        IOPT(LP+1)=1,2 or 3
!            1= Nominal scaling as noted;
!            2= Each nonzero column scaled to have length ONE;
!            3= Identity scaling; scaling effectively suppressed.
!         .
!        CALL DBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   4
!   -
!     This options allows the user to provide arbitrary (positive)
!     column scaling for the matrix. If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=4
!        IOPT(LP+1)=IOFF
!        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
!        = Positive scale factors for cols. of E.
!         .
!        CALL DBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*)
!     and NCOLS to the required length of X(*).
!
!   5
!   -
!     This option allows the user to provide an option array to the
!     low-level subprogram DBOLS( ). If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=5
!        IOPT(LP+1)= Position in IOPT(*) where option array
!                    data for DBOLS( ) begins.
!         .
!        CALL DBOCLS( )
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   6
!   -
!     This option is no longer operative.  To pass an option array
!     to the low-level subprogram DBOLSM( ), imbed it within an option
!     array passed to DBOLS() using option 5.
!
!   7
!   -
!     Move the processing pointer (either forward or backward) to the
!     location IOPT(LP+1). The processing pointer moves to locations
!     LP+2 if option number 7 is used with the value -7.  For
!     example to skip over locations 3,...,NCOLS+2,
!
!       IOPT(1)=7
!       IOPT(2)=NCOLS+3
!       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
!       IOPT(NCOLS+3)=99
!       CALL DBOCLS( )
!
!     CAUTION: Misuse of this option can yield some very hard-to-find
!     bugs. Use it with care. It is intended to be used for passing
!     option arrays to other subprograms.
!
!   8
!   -
!     This option allows the user to suppress the algorithmic feature
!     of DBOCLS( ) that processes the constraint equations C*X = Y and
!     resolves infeasibilities. The steps normally done are to solve
!     C*X - Y = 0 in a least squares sense using the stated bounds on
!     both X and Y. Then the "reachable" vector Y = C*X is computed
!     using the solution X obtained. Finally the stated bounds for Y are
!     enlarged to include C*X. To suppress the feature:
!
!
!       IOPT(LP)=8
!         .
!       CALL DBOCLS( )
!
!     Use of this option adds 1 to the required length of IOPT(*).
!
!   9
!   -
!     This option allows the user to suppress the pretriangularizing
!     step of the least squares matrix that is done within DBOCLS( ).
!     This is primarily a means of enhancing the subprogram efficiency
!     and has little effect on accuracy. To suppress the step, set:
!
!       IOPT(LP)=9
!         .
!       CALL DBOCLS( )
!
!     Use of this option adds 1 to the required length of IOPT(*).
!
!   99
!   --
!     There are no more options to change.
!
!     Only option numbers -99, -9,-8,...,-1, 1,2,...,9, and 99 are
!     permitted. Other values are errors. Options -99,-1,...,-9 mean
!     that the respective options 99,1,...,9 are left at their default
!     values. An example is the option to suppress the preprocessing of
!     contraints:
!
!       IOPT(1)=-8 Option is recognized but not changed
!       IOPT(2)=99
!       CALL DBOCLS( )
!
!***END PROLOGUE  DBOCLS
!     REVISED 880722-1100
!     REVISED YYMMDD-HHMM
!
!    PURPOSE
!    -------
!     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE LEAST SQUARES
!     PROBLEM CONSISTING OF LINEAR CONSTRAINTS
!
!              C*X = Y
!
!     AND LEAST SQUARES EQUATIONS
!
!              E*X = F
!
!     IN THIS FORMULATION THE VECTORS X AND Y ARE BOTH UNKNOWNS.
!     FURTHER, X AND Y MAY BOTH HAVE USER-SPECIFIED BOUNDS ON EACH
!     COMPONENT.  THE USER MUST HAVE DIMENSION STATEMENTS OF THE
!     FORM
!
!     DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON),BU(NCOLS+MCON),
!               X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON)
!
!     integer IND(NCOLS+MCON), IOPT(16+NI), IW(2*(NCOLS+MCON))
!
  implicit none

  integer mdw

  logical :: accum
  real ( kind = 8 ) anorm
  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) bu(*)
  !logical, save :: checkl
  logical :: checkl
  real ( kind = 8 ) cnorm
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) drelpr
  logical filter
  integer i
  integer icase
  integer idope
  integer idum
  integer :: igo = 0
  integer iiw
  integer ind(*)
  integer inrows
  integer iopt(*)
  integer ip
  integer irw
  integer iscale
  integer iw(*)
  integer j
  integer jopt(5)
  integer jp
  integer lds
  integer lenx
  integer level
  integer liopt
  integer liw
  integer llb
  integer lliw
  integer llrw
  integer llx
  integer lmdw
  integer lndw
  integer locacc
  integer locdim
  integer lopt
  integer lp
  integer lrw
  integer mcon
  integer mdwl
  integer mnew
  integer mode
  integer modec
  integer mout
  integer mrows
  integer ncols
  integer nerr
  real ( kind = 8 ), parameter :: one = 1.0D+00
  logical pretri
  real ( kind = 8 ) rdum
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) rnormc
  real ( kind = 8 ) rw(*)
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) w(mdw,*)
  real ( kind = 8 ) wt
  real ( kind = 8 ) x(*)

  save

  common /dbocom/ idope(5)

  idum = 0
  rdum = 0.0D+00
  nerr = 0
  mode = 0
  level = 1

  if ( igo == 0) then
!
!  Check validity of input data.
!
!  SEE THAT MDW IS .GT.0. GROSS CHECK ONLY.
!
      if ( mdw<=0) then
          nerr = 53
          call xerrwv('dbocls(). mdw=(i1) must be positive.', &
                      nerr,level,1,mdw,idum,0,rdum,rdum)
          go to 260
      end if
!
!  SEE THAT NUMBER OF CONSTRAINTS IS NONNEGATIVE.
!
      if ( mcon < 0) then
          nerr = 54
          call xerrwv('dbocls(). mcon=(i1) must be nonnegative.', &
                      nerr,level,1,mcon,idum,0,rdum,rdum)
          go to 260
      end if
!
!  SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
!
      if ( ncols<=0) then
          nerr = 55
          call xerrwv( &
       'dbocls(). ncols=(i1) the no. of variables must be positive.' &
                      ,nerr,level,1,ncols,idum,0,rdum,rdum)
          go to 260
      end if
!
!  SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
!
      do j = 1,ncols + mcon
          if ( ind(j) < 1 .or. ind(j) > 4) then
              nerr = 56
              call xerrwv( &
                    'dbocls(). for j=(i1), ind(j)=(i2) must be 1-4.' &
                          ,nerr,level,2,j,ind(j),0,rdum,rdum)
              go to 260
          end if
      end do
!
!  SEE THAT BOUNDS ARE CONSISTENT.
!
      do j = 1,ncols + mcon
          if ( ind(j) == 3) then
              if ( bl(j) > bu(j)) then
                  nerr = 57
                  call xerrwv( &
        'dbocls(). for j=(i1), bound bl(j)=(r1) is  >  bu(j)=(r2).' &
         ,nerr,level,1,j,idum,2,bl(j),bu(j))
                  go to 260
              end if
          end if
      end do
!
!  PROCESS OPTION ARRAY
!
      drelpr = epsilon ( drelpr )
      checkl = .false.
      filter = .true.
      lenx = 2* (ncols+mcon) + 2
      iscale = 1
      igo = 1
      accum = .false.
      pretri = .true.
      lopt = 0
      lp = 0
      lds = 0

   30     continue
      lp = lp + lds
      ip = iopt(lp+1)
      jp = abs(ip)
!
!  TEST FOR NO MORE OPTIONS TO CHANGE.
!
      if ( ip == 99) then
          if ( lopt == 0) lopt = lp+1
!
!  SEND COL. SCALING TO DBOLS().
!
          idope(4)=1
!
!  NOTE THAT DBOLS() WAS CALLED BY DBOCLS()
!
          idope(5)=1
!
!  CHANGE PRETRIANGULARIZATION FACTOR IN DBOLSM().
!
          idope(1) = ncols + mcon + 1
!
!  PASS WEIGHT TO DBOLSM() FOR RANK TEST.
!
          idope(2) = ncols + mcon + 2
          idope(3) = mcon
          go to 50
      else if ( jp == 99) then
          lds = 1
          go to 50
      else if ( jp == 1) then
          if ( ip > 0) then
!
!  SET UP DIRECTION FLAG LOCATION, ROW STACKING POINTER
!  LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
!
              locacc = lp + 2
!
!                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
!     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
!                  IOPT(LOCACC+1)=ROW STACKING POINTER.
!                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
!
!     USER ACTION WITH THIS OPTION..
!
!      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).)
!      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
!       ROW OF W(*,*) BELOW THE ROWS FOR THE CONSTRAINT MATRIX C.
!       SET IOPT(LOCACC+2)=NO. OF LEAST SQUARES EQUATIONS IN BLOCK.
!
!              LOOP
!              CALL DBOCLS()
!
!                  IF(IOPT(LOCACC) .EQ. 1) THEN
!                      STACK EQUAS. INTO W(*,*), STARTING AT
!                      ROW IOPT(LOCACC+1).
!                       INTO W(*,*).
!                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
!                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
!                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN
!                      (PROCESS IS OVER. EXIT LOOP.)
!                  ELSE
!                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
!                  END IF
!              END LOOP
!
              iopt(locacc+1) = mcon + 1
              accum = .true.
              iopt(locacc) = igo
          end if
          lds = 4
          go to 30
      else if ( jp == 2) then
          if ( ip > 0) then
!
!  GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
!
              locdim = lp + 2
!
!  LMDW.GE.MCON+MAX(MOUT,NCOLS), IF MCON.GT.0 .AND FILTER
!  LMDW.GE.MCON+MOUT, OTHERWISE
!
!  LNDW.GE.NCOLS+MCON+1
!  LLB .GE.NCOLS+MCON
!  LLX .GE.2*(NCOLS+MCON)+2+EXTRA REQD. IN OPTIONS.
!  LLRW.GE.6*NCOLS+5*MCON
!  LLIW.GE.2*(NCOLS+MCON)
!  LIOP.GE. AMOUNT REQD. FOR OPTION ARRAY.
!
              lmdw = iopt(locdim)
              lndw = iopt(locdim+1)
              llb = iopt(locdim+2)
              llx = iopt(locdim+3)
              llrw = iopt(locdim+4)
              lliw = iopt(locdim+5)
              liopt = iopt(locdim+6)
              checkl = .true.
          end if
          lds = 8
          go to 30
!
!  OPTION TO MODIFY THE COLUMN SCALING.
!
      else if ( jp == 3) then
          if ( ip > 0) then
              iscale = iopt(lp+2)
!
!     SEE THAT ISCALE IS 1 THRU 3.
!
              if ( iscale < 1 .or. iscale > 3) then
                  nerr = 48
                  call xerrwv('dbocls(). iscale option=(i1) must be 1-3.' &
                    ,nerr,level,1,iscale,idum,0,rdum,rdum)
                  go to 260
              end if
          end if
          lds = 2
          go to 30
!
!  IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE
!  SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)).
!
      else if ( jp == 4) then
          if ( ip > 0) then
              iscale = 4
              if ( iopt(lp+2)<=0) then
                  nerr = 49
                  call xerrwv('dbocls(). offset past x(ncols) (i1) for' // &
                    'user-provided column scaling must be positive.', &
                    nerr,level,1,iopt(lp+2),idum,0,rdum,rdum)
                  go to 260
              end if
              call dcopy(ncols,x(ncols+iopt(lp+2)),1,rw,1)
              lenx = lenx + ncols
              do j = 1,ncols
                  if ( rw(j)<=0.0D+00 ) then
                      nerr = 50
                      call xerrwv('dbocls(). each provided col. scale ' // &
                        'factor must be positive. comp. (i1)   now = (r1).', &
                        nerr,level,1,j,idum,1,rw(j),rdum)
                      go to 260
                  end if
              end do
          end if
          lds = 2
          go to 30
!
!  IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLS().
!
      else if ( jp == 5) then
          if ( ip > 0) then
              lopt = iopt(lp+2)
          end if
          lds = 2
          go to 30
!
!  IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM().
!  (NO LONGER USED.) OPTION NOW MUST BE PASSED IMBEDDED IN
!  OPTION ARRAY FOR DBOLS().
!
      else if ( jp == 6) then
          lds = 2
          go to 30
!
!  THIS OPTION USES THE NEXT LOC OF IOPT(*) AS A
!  POINTER VALUE TO SKIP TO NEXT.
!
      else if ( jp == 7) then
          if ( ip > 0) then
              lp = iopt(lp+2)-1
              lds = 0
          else
              lds = 2
          end if
          go to 30
!
!  THIS OPTION AVOIDS THE CONSTRAINT RESOLVING PHASE FOR
!  THE LINEAR CONSTRAINTS C*X=Y.
!
      else if ( jp == 8) then
          filter = .not. (ip > 0)
          lds = 1
          go to 30
!
!  THIS OPTION SUPPRESSES PRETRIANGULARIZATION OF THE LEAST
!  SQUARES EQATIONS.
!
      else if ( jp == 9) then
          pretri = .not. (ip > 0)
          lds = 1
          go to 30
!
!  NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
!
      else
          nerr = 51
          call xerrwv('dbocls(). the option number=(i1) is not defined.' &
            ,nerr,level,1,jp,idum,0,rdum,rdum)
          go to 260
      end if

   50     continue

      if ( checkl) then
!
!  CHECK LENGTHS OF ARRAYS
!
!  THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
!  ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
!
       if ( filter .and. .not.accum) then
         mdwl=mcon+max(mrows,ncols)
       else if ( accum) then
         mdwl=mcon+ncols+1
       else
         mdwl=mcon+ncols
       end if

          if ( lmdw < mdwl) then
              nerr = 41
              call xerrwv('dbocls(). the row dimension of w(,)=(i1) ' // &
                'must be >= the number of effective rows=(i2).', &
                nerr,level,2,lmdw,mdwl,0,rdum,rdum)
              go to 260
          end if
          if ( lndw < ncols+mcon+1) then
              nerr = 42
              call xerrwv('dbocls(). the column dimension of w(,)=(i1) ' // &
                'must be >= ncols+mcon+1=(i2).',nerr,level,2,lndw, &
                ncols+mcon+1,0,rdum,rdum)
              go to 260
          end if
          if ( llb < ncols+mcon) then
              nerr = 43
              call xerrwv('dbocls(). the dimensions of the arrays bl(),' // &
                'bu(), and ind()=(i1) must be >= ncols+mcon=(i2).', &
                nerr,level,2,llb,ncols+mcon,0,rdum,rdum)
              go to 260
          end if

          if ( llx < lenx) then
              nerr = 44
              call xerrwv( 'dbocls(). the dimension of x()=(i1) must be ' // &
                '>= the reqd. length=(i2).',nerr,level,2,llx,lenx, &
                0,rdum,rdum)
              go to 260
          end if

          if ( llrw < 6*ncols+5*mcon) then
              nerr = 45
              call xerrwv('dbocls(). the dimension of rw()=(i1) must be ' // &
                '>= 6*ncols+5*mcon=(i2).',nerr,level,2, &
                llrw,6*ncols+5*mcon,0,rdum,rdum)
              go to 260
          end if

          if ( lliw < 2*ncols+2*mcon) then
              nerr = 46
              call xerrwv('dbocls() the dimension of iw()=(i1) must be ' // &
                '>= 2*ncols+2*mcon=(i2).',nerr,level,2,lliw, &
                2*ncols+2*mcon,0,rdum,rdum)
              go to 260
          end if

          if ( liopt < lp+17) then
              nerr = 47
              call xerrwv('dbocls(). the dimension of iopt()=(i1) must be ' // &
                '>= the reqd. len.=(i2).',nerr,level,2,liopt,lp+17, 0,rdum,rdum)
              go to 260
          end if
      end if
  end if
!
!  OPTIONALLY GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
!  EQUATIONS AND DIRECTIONS FOR PROCESSING THESE EQUATIONS.
!
!  ACCUMULATE LEAST SQUARES EQUATIONS
!
  if ( accum) then
      mrows = iopt(locacc+1) - 1 - mcon
      inrows = iopt(locacc+2)
      mnew = mrows + inrows
      if ( mnew < 0 .or. mnew+mcon > mdw) then
          nerr = 52
          call xerrwv('dbocls(). no. of rows=(i1) must be >= 0 ' // &
            '.and. <=mdw-mcon=(i2)',nerr,level,2,mnew,mdw-mcon,0,rdum,rdum)
          go to 260
      end if
  end if
!
!  USE THE SOFTWARE OF DBOLS( ) FOR THE TRIANGULARIZATION OF THE
!  LEAST SQUARES MATRIX.  THIS MAY INVOLVE A SYSTALTIC INTERCHANGE
!  OF PROCESSING POINTERS BETWEEN THE CALLING AND CALLED (DBOLS())
!  PROGRAM UNITS.
!
  jopt(01) = 1
  jopt(02) = 2
  jopt(04) = mrows
  jopt(05) = 99
  irw = ncols + 1
  iiw = 1
  if ( accum .or. pretri) then
!
!  NOTE THAT DBOLS() WAS CALLED BY DBOCLS()
!
          idope(5)=0
      call dbols(w(mcon+1,1),mdw,mout,ncols,bl,bu,ind,jopt,x,rnorm, &
             mode,rw(irw),iw(iiw))
  else
      mout = mrows
  end if

  if ( accum) then
    accum = iopt(locacc)  ==  1
    iopt(locacc+1) = jopt(03) + mcon
    mrows = min(ncols+1,mnew)
  end if

  if ( accum) return
!
!  SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM
!
!  MOVE RIGHT HAND SIDE OF LEAST SQUARES EQUATIONS.
!
  call dcopy(mout,w(mcon+1,ncols+1),1,w(mcon+1,ncols+mcon+1),1)
  if ( mcon > 0 .and. filter) then
!
!  PROJECT THE LINEAR CONSTRAINTS INTO A REACHABLE SET.
!
      do i = 1,mcon
        call dcopy(ncols,w(i,1),mdw,w(mcon+1,ncols+i),1)
      end do
!
!  PLACE (-)IDENTITY MATRIX AFTER CONSTRAINT DATA.
!
      do j = ncols + 1, ncols + mcon + 1
        w(1,j) = 0.0D+00
        call dcopy ( mcon, w(1,j), 0, w(1,j), 1 )
      end do

      w(1,ncols+1) = -1.0D+00
      call dcopy ( mcon, w(1,ncols+1), 0, w(1,ncols+1), mdw+1 )
!
!  OBTAIN A 'FEASIBLE POINT' FOR THE LINEAR CONSTRAINTS.
!
      jopt(01) = 99
      irw = ncols + 1
      iiw = 1
!
!  NOTE THAT DBOLS() WAS CALLED BY DBOCLS()
!
          idope(5)=0
      call dbols(w,mdw,mcon,ncols+mcon,bl,bu,ind,jopt,x,rnormc, &
        modec,rw(irw),iw(iiw))
!
!  ENLARGE THE BOUNDS SET, IF REQUIRED, TO INCLUDE POINTS THAT
!  CAN BE REACHED.
!
      do j = ncols + 1,ncols + mcon
          icase = ind(j)
          if ( icase < 4) then
              t = ddot ( ncols, w(mcon+1,j), 1, x, 1 )
          end if
          go to (80,90,100,110),icase
          go to 120
   80         bl(j) = min(t,bl(j))
          go to 120
   90         bu(j) = max(t,bu(j))
          go to 120
  100         bl(j) = min(t,bl(j))
          bu(j) = max(t,bu(j))
          go to 120
  110         continue
  120         continue
      end do
!
!  MOVE CONSTRAINT DATA BACK TO THE ORIGINAL AREA.
!
      do j = ncols + 1,ncols + mcon
          call dcopy(ncols,w(mcon+1,j),1,w(j-ncols,1),mdw)
      end do

  end if

  if ( mcon > 0) then
      do j = ncols + 1,ncols + mcon
        w(mcon+1:mcon+mout,j) = 0.0D+00
      end do
!
!  PUT IN (-)IDENTITY MATRIX (POSSIBLY) ONCE AGAIN.
!
      do j = ncols + 1,ncols + mcon + 1
        w(1,j) = 0.0D+00
        call dcopy ( mcon, w(1,j), 0, w(1,j), 1 )
      end do

      w(1,ncols+1) = -1.0D+00
      call dcopy ( mcon, w(1,ncols+1), 0, w(1,ncols+1), mdw+1 )

  end if
!
!  COMPUTE NOMINAL COLUMN SCALING FOR THE UNWEIGHTED MATRIX.
!
  cnorm = 0.0D+00
  anorm = 0.0D+00
  do j = 1,ncols
      t1 = dasum(mcon,w(1,j),1)
      t2 = dasum(mout,w(mcon+1,1),1)
      t = t1 + t2
      if ( t == 0.0D+00 ) t = one
      cnorm = max(cnorm,t1)
      anorm = max(anorm,t2)
      x(ncols+mcon+j) = one/t
  end do

  go to (180,190,210,220),iscale
  go to 230
  180 continue
  go to 230
!
!  SCALE COLS. (BEFORE WEIGHTING) TO HAVE LENGTH ONE.
!
  190 continue

  do j = 1,ncols
    t = dnrm2(mcon+mout,w(1,j),1)
    if ( t == 0.0D+00 ) t = one
    x(ncols+mcon+j) = one/t
  end do

  go to 230
!
!  SUPPRESS SCALING (USE UNIT MATRIX).
!
  210 continue

  x(ncols+mcon+1) = one
  call dcopy(ncols,x(ncols+mcon+1),0,x(ncols+mcon+1),1)
  go to 230
!
!  THE USER HAS PROVIDED SCALING.
!
  220 call dcopy(ncols,rw,1,x(ncols+mcon+1),1)
  230 continue

  do j = ncols + 1,ncols + mcon
    x(ncols+mcon+j) = one
  end do
!
!  WEIGHT THE LEAST SQUARES EQUATIONS.
!
  wt = sqrt(drelpr)
  if ( anorm > 0.0D+00 ) wt = wt/anorm
  if ( cnorm > 0.0D+00 ) wt = wt*cnorm

  do i = 1,mout
      call dscal(ncols,wt,w(i+mcon,1),mdw)
  end do
  call dscal(mout,wt,w(mcon+1,mcon+ncols+1),1)
  lrw = 1
  liw = 1
!
!  SET THE NEW TRIANGULARIZATION FACTOR.
!
  x(ncols+mcon+idope(1))= 0.0D+00
!
!  SET THE WEIGHT TO USE IN COMPONENTS .GT. MCON,
!  WHEN MAKING LINEAR INDEPENDENCE TEST.
!
  x(ncols+mcon+idope(2))= one/wt
  idope(5)=1
  call dbols(w,mdw,mout+mcon,ncols+mcon,bl,bu,ind,iopt(lopt),x, &
    rnorm,mode,rw(lrw),iw(liw))

  260 continue

  if ( mode >= 0 ) then
    mode=-nerr
  end if

  igo = 0

  return
end
subroutine dbols ( w, mdw, mrows, ncols, bl, bu, ind, iopt, x, rnorm, &
  mode, rw, iw )

!***********************************************************************
!
!! DBOLS solves the linear system E*X = F in the least squares sense.
!
!  Discussion:
!
!    This routine solves the problem
!
!      E*X = F 
!
!    in the least squares sense with bounds on selected X values.
!
!  Author:
!
!    Richard Hanson, Sandia National Laboratory
!
!  Reference:
!
!    Richard Hanson,
!    Linear Least Squares with Bounds and Linear Constraints,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 7, Number 3, July 1986.
!
!***DESCRIPTION
!
!     The user must have dimension statements of the form:
!
!       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
!      * X(NCOLS+NX), RW(5*NCOLS)
!       integer IND(NCOLS), IOPT(1+NI), IW(2*NCOLS)
!
!     (here NX=number of extra locations required for option 4; NX=0
!     for no options; NX=NCOLS if this option is in use. Here NI=number
!     of extra locations required for options 1-6; NI=0 for no
!     options.)
!
!   INPUT
!   -----
!
!    --------------------
!    W(MDW,*),MROWS,NCOLS
!    --------------------
!     The array W(*,*) contains the matrix [E:F] on entry. The matrix
!     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in
!     the array W(*,*) with E occupying the first NCOLS columns and the
!     right side vector F in column NCOLS+1. The row dimension, MDW, of
!     the array W(*,*) must satisfy the inequality MDW >= MROWS.
!     Other values of MDW are errrors. The values of MROWS and NCOLS
!     must be positive. Other values are errors. There is an exception
!     to this when using option 1 for accumulation of blocks of
!     equations. In that case MROWS is an OUTPUT variable ONLY, and the
!     matrix data for [E:F] is placed in W(*,*), one block of rows at a
!     time.  MROWS contains the number of rows in the matrix after
!     triangularizing several blocks of equations. This is an OUTPUT
!     parameter ONLY when option 1 is used. See IOPT(*) CONTENTS
!     for details about option 1.
!
!    ------------------
!    BL(*),BU(*),IND(*)
!    ------------------
!     These arrays contain the information about the bounds that the
!     solution values are to satisfy. The value of IND(J) tells the
!     type of bound and BL(J) and BU(J) give the explicit values for
!     the respective upper and lower bounds.
!
!    1.    For IND(J)=1, require X(J) >= BL(J).
!          (the value of BU(J) is not used.)
!    2.    For IND(J)=2, require X(J) <= BU(J).
!          (the value of BL(J) is not used.)
!    3.    For IND(J)=3, require X(J) >= BL(J) and
!                                X(J) <= BU(J).
!    4.    For IND(J)=4, no bounds on X(J) are required.
!          (the values of BL(J) and BU(J) are not used.)
!
!     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
!     IND(J)=3 (upper and lower bounds) the condition BL(J)  >  BU(J)
!     is an error.
!
!    -------
!    IOPT(*)
!    -------
!     This is the array where the user can specify nonstandard options
!     for DBOLSM( ). Most of the time this feature can be ignored by
!     setting the input value IOPT(1)=99. Occasionally users may have
!     needs that require use of the following subprogram options. For
!     details about how to use the options see below: IOPT(*) CONTENTS.
!
!     Option Number   Brief Statement of Purpose
!     ------ ------   ----- --------- -- -------
!           1         Return to user for accumulation of blocks
!                     of least squares equations.
!           2         Check lengths of all arrays used in the
!                     subprogram.
!           3         Standard scaling of the data matrix, E.
!           4         User provides column scaling for matrix E.
!           5         Provide option array to the low-level
!                     subprogram DBOLSM( ).
!           6         Move the IOPT(*) processing pointer.
!           7         User has called DBOLS() directly.
!          99         No more options to change.
!
!    ----
!    X(*)
!    ----
!     This array is used to pass data associated with option 4. Ignore
!     this parameter if this option is not used. Otherwise see below:
!     IOPT(*) CONTENTS.
!
!    OUTPUT
!    ------
!
!    ----------
!    X(*),RNORM
!    ----------
!     The array X(*) contains a solution (if MODE >=0 or  == -22) for
!     the constrained least squares problem. The value RNORM is the
!     minimum residual vector length.
!
!    ----
!    MODE
!    ----
!     The sign of MODE determines whether the subprogram has completed
!     normally, or encountered an error condition or abnormal status. A
!     value of MODE >= 0 signifies that the subprogram has completed
!     normally. The value of MODE (.GE. 0) is the number of variables
!     in an active status: not at a bound nor at the value ZERO, for
!     the case of free variables. A negative value of MODE will be one
!     of the cases -37,-36,...,-22, or -17,...,-2. Values  <  -1
!     correspond to an abnormal completion of the subprogram. To
!     understand the abnormal completion codes see below: ERROR
!     MESSAGES for DBOLS( ). AN approximate solution will be returned
!     to the user only when max. iterations is reached, MODE=-22.
!     Values for MODE=-37,...,-22 come from the low-level subprogram
!     DBOLSM(). See the section ERROR MESSAGES for DBOLSM() in the
!     documentation for DBOLSM().
!
!    -----------
!    RW(*),IW(*)
!    -----------
!     These are working arrays with 5*NCOLS and 2*NCOLS entries.
!     (normally the user can ignore the contents of these arrays,
!     but they must be dimensioned properly.)
!
!    IOPT(*) CONTENTS
!    ------- --------
!     The option array allows a user to modify internal variables in
!     the subprogram without recompiling the source code. A central
!     goal of the initial software design was to do a good job for most
!     people. Thus the use of options will be restricted to a select
!     group of users. The processing of the option array proceeds as
!     follows: a pointer, here called LP, is initially set to the value
!     1. This value is updated as each option is processed. At the
!     pointer position the option number is extracted and used for
!     locating other information that allows for options to be changed.
!     The portion of the array IOPT(*) that is used for each option is
!     fixed; the user and the subprogram both know how many locations
!     are needed for each option. A great deal of error checking is
!     done by the subprogram on the contents of the option array.
!     Nevertheless it is still possible to give the subprogram optional
!     input that is meaningless. For example option 4 uses the
!     locations X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing
!     scaling data. The user must manage the allocation of these
!     locations.
!
!   1
!   -
!     This option allows the user to solve problems with a large number
!     of rows compared to the number of variables. The idea is that the
!     subprogram returns to the user (perhaps many times) and receives
!     new least squares equations from the calling program unit.
!     Eventually the user signals "that's all" and then computes the
!     solution with one final call to subprogram DBOLS( ). The value of
!     MROWS is an OUTPUT variable when this option is used. Its value
!     is always in the range 0 <= MROWS .le. NCOLS+1. It is equal to
!     the number of rows after the triangularization of the entire set
!     of equations. If LP is the processing pointer for IOPT(*), the
!     usage for the sequential processing of blocks of equations is
!
!        IOPT(LP)=1
!        Move block of equations to W(*,*) starting at
!        the first row of W(*,*).
!        IOPT(LP+3)=# of rows in the block; user defined
!
!     The user now calls DBOLS( ) in a loop. The value of IOPT(LP+1)
!     directs the user's action. The value of IOPT(LP+2) points to
!     where the subsequent rows are to be placed in W(*,*).
!
!      .<LOOP
!      . CALL DBOLS()
!      . IF(IOPT(LP+1) .EQ. 1) THEN
!      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
!      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
!      .    W(*,*) STARTING AT ROW IOPT(LP+2).
!      .
!      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN
!      .       IOPT(LP+1)=2
!      .<------CYCLE LOOP
!      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN
!      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0
!      . ELSE
!      . ERROR CONDITION; SHOULD NOT HAPPEN.
!      .<END LOOP
!
!     Use of this option adds 4 to the required length of IOPT(*).
!
!
!   2
!   -
!     This option is useful for checking the lengths of all arrays used
!     by DBOLS() against their actual requirements for this problem.
!     The idea is simple: the user's program unit passes the declared
!     dimension information of the arrays. These values are compared
!     against the problem-dependent needs within the subprogram. If any
!     of the dimensions are too small an error message is printed and a
!     negative value of MODE is returned, -11 to -17. The printed error
!     message tells how long the dimension should be. If LP is the
!     processing pointer for IOPT(*),
!
!        IOPT(LP)=2
!        IOPT(LP+1)=Row dimension of W(*,*)
!        IOPT(LP+2)=Col. dimension of W(*,*)
!        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*)
!        IOPT(LP+4)=Dimension of X(*)
!        IOPT(LP+5)=Dimension of RW(*)
!        IOPT(LP+6)=Dimension of IW(*)
!        IOPT(LP+7)=Dimension of IOPT(*)
!         .
!        CALL DBOLS()
!
!     Use of this option adds 8 to the required length of IOPT(*).
!
!   3
!   -
!     This option changes the type of scaling for the data matrix E.
!     Nominally each nonzero column of E is scaled so that the
!     magnitude of its largest entry is equal to the value ONE. If LP
!     is the processing pointer for IOPT(*),
!
!        IOPT(LP)=3
!        IOPT(LP+1)=1,2 or 3
!            1= Nominal scaling as noted;
!            2= Each nonzero column scaled to have length ONE;
!            3= Identity scaling; scaling effectively suppressed.
!         .
!        CALL DBOLS()
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   4
!   -
!     This option allows the user to provide arbitrary (positive)
!     column scaling for the matrix E. If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=4
!        IOPT(LP+1)=IOFF
!        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
!        = Positive scale factors for cols. of E.
!         .
!        CALL DBOLS()
!
!     Use of this option adds 2 to the required length of IOPT(*) and
!     NCOLS to the required length of X(*).
!
!   5
!   -
!     This option allows the user to provide an option array to the
!     low-level subprogram DBOLSM(). If LP is the processing pointer
!     for IOPT(*),
!
!        IOPT(LP)=5
!        IOPT(LP+1)= Position in IOPT(*) where option array
!                    data for DBOLSM() begins.
!         .
!        CALL DBOLS()
!
!     Use of this option adds 2 to the required length of IOPT(*).
!
!   6
!   -
!     Move the processing pointer (either forward or backward) to the
!     location IOPT(LP+1). The processing point is moved to entry
!     LP+2 of IOPT(*) if the option is left with -6 in IOPT(LP).  For
!     example to skip over locations 3,...,NCOLS+2 of IOPT(*),
!
!       IOPT(1)=6
!       IOPT(2)=NCOLS+3
!       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
!       IOPT(NCOLS+3)=99
!       CALL DBOLS()
!
!     CAUTION: Misuse of this option can yield some very hard
!     -to-find bugs.  Use it with care.
!
!   7
!   -
!     If the user is calling DBOLS() directly, use this option.
!     (This is necessary because DBOCLS() uses DBOLS() as a
!     low-level subprogram.  Due to weighting required within
!     DBOCLS(), the two cases must be known.) For example,
!
!       IOPT(1)=7
!       IOPT(1)=99
!
!   99
!   --
!     There are no more options to change.
!
!     Only option numbers -99, -7,-6,-5,...,-1, 1,2,...,7, and 99 are
!     permitted. Other values are errors. Options -99,-1,...,-7 mean
!     that the repective options 99,1,...,7 are left at their default
!     values. An example is the option to modify the (rank) tolerance:
!
!       IOPT(1)=-3 Option is recognized but not changed
!       IOPT(2)=2  Scale nonzero cols. to have length ONE
!       IOPT(3)=99
!
!***ROUTINES CALLED  IDAMAX,DBOLSM,DCOPY,DNRM2,DROT,DROTG,XERRWV
!***COMMON BLOCKS    DBOCOM
!***END PROLOGUE  DBOLS
!
!     SOLVE LINEAR LEAST SQUARES SYSTEM WITH BOUNDS ON
!     SELECTED VARIABLES.
!
  implicit none

  integer mdw

  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) bu(*)
  logical checkl
  real ( kind = 8 ) dnrm2
  integer i
  integer ibig
  integer idamax
  integer idope
  integer idum
  integer :: igo = 0
  integer ind(*)
  integer inrows
  integer iopt(*)
  integer ip
  integer :: iscale
  integer iw(*)
  integer j
  integer jp
  integer lds
  integer lenx
  integer level
  integer liopt
  integer llb
  integer lliw
  integer llrw
  integer llx
  integer lmdw
  integer lndw
  integer :: locacc
  integer locdim
  integer :: lopt
  integer lp
  integer mnew
  integer mode
  integer mrows
  integer ncols
  integer nerr
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) rdum
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) rw(*)
  real ( kind = 8 ) sc
  real ( kind = 8 ) ss
  real ( kind = 8 ) w(mdw,*)
  real ( kind = 8 ) x(*)

  save

  common /dbocom/ idope(5)

  idum = 0
  rdum = 0.0D+00
  nerr = 0
  mode = 0
  level = 1
  if ( igo == 0 ) then
!
!  CHECK VALIDITY OF INPUT DATA
!
!  SEE THAT MDW IS .GT.0. GROSS CHECK ONLY.
!
      if ( mdw<=0) then
          nerr = 2
          call xerrwv('dbols(). mdw=(i1) must be positive.', &
            nerr,level,1,mdw,idum,0,rdum,rdum)
          go to 190
      end if
!
!  SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
!
      if ( ncols<=0) then
          nerr = 3
          call xerrwv('dbols(). ncols=(i1) the no. of variables must ' // &
            'be positive.',nerr,level,1,ncols,idum,0,rdum,rdum)
          go to 190
      end if
!
!  SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
!
      do j = 1,ncols
          if ( ind(j) < 1 .or. ind(j) > 4) then
              nerr = 4
              call xerrwv( &
                     'dbols(). for j=(i1), ind(j)=(i2) must be 1-4.' &
                          ,nerr,level,2,j,ind(j),0,rdum,rdum)
              go to 190
          end if
      end do
!
!  SEE THAT BOUNDS ARE CONSISTENT.
!
      do j = 1,ncols
          if ( ind(j) == 3) then
              if ( bl(j) > bu(j)) then
                  nerr = 5
                  call xerrwv( &
         'dbols(). for j=(i1), bound bl(j)=(r1) is  >  bu(j)=(r2).' &
                   ,nerr,level,1,j,idum,2,bl(j),bu(j))
                  go to 190
              end if
          end if
      end do
!
!  PROCESS OPTION ARRAY
!
      checkl = .false.
      lenx = ncols
      iscale = idope(4)
      igo = 2
      lopt = 0
      lp = 0
      lds = 0

30    continue

      lp = lp + lds
      ip = iopt(lp+1)
      jp = abs(ip)
!
!  TEST FOR NO MORE OPTIONS.
!
      if ( ip == 99) then
          if ( lopt == 0) lopt = lp + 1
          go to 50
      else if ( jp == 99) then
          lds = 1
          go to 30
      else if ( jp == 1) then
          if ( ip > 0) then
!
!  SET UP DIRECTION FLAG, ROW STACKING POINTER
!  LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
!
              locacc = lp + 2
!
!                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
!     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
!                  IOPT(LOCACC+1)=ROW STACKING POINTER.
!                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
!     USER ACTION WITH THIS OPTION..
!      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).
!      MUST ALSO START PROCESS WITH IOPT(LOCACC)=1.)
!      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
!       ROW OF W(*,*).  SET IOPT(LOCACC+2)=NO. OF ROWS IN BLOCK.)
!              LOOP
!              CALL DBOLS()
!
!                  IF(IOPT(LOCACC) .EQ. 1) THEN
!                      STACK EQUAS., STARTING AT ROW IOPT(LOCACC+1),
!                       INTO W(*,*).
!                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
!                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
!                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN
!                      (PROCESS IS OVER. EXIT LOOP.)
!                  ELSE
!                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
!                  END IF
!              END LOOP
!              SET IOPT(LOCACC-1)=-OPTION NUMBER FOR SEQ. ACCUMULATION.
!              CALL DBOLS( )
              iopt(locacc+1) = 1
              igo = 1
          end if
          lds = 4
          go to 30
      else if ( jp == 2) then
          if ( ip > 0) then
!
!  GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
!
              locdim = lp + 2
!
!  LMDW.GE.MROWS
!  LNDW.GE.NCOLS+1
!  LLB .GE.NCOLS
!  LLX .GE.NCOLS+EXTRA REQD. IN OPTIONS.
!  LLRW.GE.5*NCOLS
!  LLIW.GE.2*NCOLS
!  LIOP.GE. AMOUNT REQD. FOR IOPTION ARRAY.
!
              lmdw = iopt(locdim)
              lndw = iopt(locdim+1)
              llb = iopt(locdim+2)
              llx = iopt(locdim+3)
              llrw = iopt(locdim+4)
              lliw = iopt(locdim+5)
              liopt = iopt(locdim+6)
              checkl = .true.
          end if
          lds = 8
          go to 30
!
!  OPTION TO MODIFY THE COLUMN SCALING.
!
      else if ( jp == 3) then
          if ( ip > 0) then
              iscale = iopt(lp+2)
!
!  SEE THAT ISCALE IS 1 THRU 3.
!
              if ( iscale < 1 .or. iscale > 3) then
                  nerr = 7
                  call xerrwv( 'dbols(). iscale option=(i1) must be 1-3.' &
                    ,nerr,level,1,iscale,idum,0,rdum,rdum)
                  go to 190
              end if
          end if
          lds = 2
          go to 30
!
!  IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE
!  SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)).
!
      else if ( jp == 4) then
          if ( ip > 0) then
              iscale = 4
              if ( iopt(lp+2)<=0) then
                  nerr = 8
                  call xerrwv('dbols(). offset past x(ncols) (i1) ' // &
                    'for user-provided column scaling must be positive.', &
                    nerr,level,1,iopt(lp+2),idum,0,rdum,rdum)
                  go to 190
              end if
              call dcopy(ncols,x(ncols+iopt(lp+2)),1,rw,1)
              lenx = lenx + ncols
              do j = 1,ncols
                  if ( rw(j)<= 0.0D+00 ) then
                      nerr = 9
                      call xerrwv('dbols(). each provided col. scale factor ' // &
                        'must be positive. component (i1) now = (r1).', &
                        nerr,level,1,j,idum,1,rw(j),rdum)
                      go to 190
                  end if
              end do
          end if
          lds = 2
          go to 30
!
!  IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM().
!
      else if ( jp == 5) then
          if ( ip > 0) then
              lopt = iopt(lp+2)
          end if
          lds = 2
          go to 30
!
!  THIS OPTION USES THE NEXT LOC OF IOPT(*) AS THE PLACE TO
!  MOVE AT THE NEXT STEP OF PROCESSESING.
!
      else if ( jp == 6) then
          if ( ip > 0) then
              lp = iopt(lp+2)-1
              lds = 0
          else
              lds = 2
          end if
          go to 30
!
!  THIS OPTION PROVIDES INFORMATION ABOUT WHO CALLED DBOLS.
!  IT WAS EITHER DBOCLS() OR THE USER.

      else if ( jp == 7) then
          lds=1
          if ( ip > 0) then
            idope(5)=0
            iscale=1
          end if
          go to 30
!
!  NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
!
      else
          nerr = 6
          call xerrwv('dbols(). the option number=(i1) is not defined.' &
            ,nerr,level,1,jp,idum,0,rdum,rdum)
          go to 190
      end if
   50     continue

      if ( checkl) then
!
!  CHECK LENGTHS OF ARRAYS
!
!  THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
!  ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
!
          if ( lmdw < mrows) then
              nerr = 11
              call xerrwv('dbols(). the row dimension of w(,)=(i1) ' // &
                'must be >=the number of rows=(i2).',nerr,level,2, &
                lmdw,mrows,0,rdum,rdum)
              go to 190
          end if
          if ( lndw < ncols+1) then
              nerr = 12
              call xerrwv('dbols(). the column dimension of w(,)=(i1) ' // &
                'must be >= ncols+1=(i2).',nerr,level,2,lndw,ncols+1,0,rdum,rdum)
              go to 190
          end if
          if ( llb < ncols) then
              nerr = 13
              call xerrwv('dbols(). the dimensions of the arrays bl(),bu(), ' // &
                'and ind()=(i1) must be >= ncols=(i2).', &
                nerr,level,2,llb,ncols,0,rdum,rdum)
              go to 190
          end if
          if ( llx < lenx) then
              nerr = 14
              call xerrwv('dbols(). the dimension of x()=(i1) must be ' // &
                '>= the reqd. length=(i2).',nerr,level,2,llx,lenx,0,rdum,rdum)
              go to 190
          end if
          if ( llrw < 5*ncols) then
              nerr = 15
              call xerrwv('dbols(). the dimension of rw()=(i1) must be ' // &
                '>= 5*ncols=(i2).',nerr,level,2,llrw,5*ncols,0,rdum,rdum)
              go to 190
          end if
          if ( lliw < 2*ncols) then
              nerr = 16
              call xerrwv('dbols() the dimension of iw()=(i1) must be ' // &
                '>= 2*ncols=(i2).',nerr,level,2,lliw,2*ncols,0,rdum,rdum)
              go to 190
          end if
          if ( liopt < lp+1) then
              nerr = 17
              call xerrwv('dbols(). the dimension of iopt()=(i1) must be ' // &
                '>= the reqd. len.=(i2).',nerr,level,2,liopt,lp+1,0,rdum,rdum)
              go to 190
          end if

      end if
  end if
  go to (60,90),igo
  go to 180
!
!  GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
!  EQUATIONS AND DIRECTIONS TO QUIT PROCESSING.
!
   60 continue
!
!  ACCUMULATE LEAST SQUARES EQUATIONS
!
  mrows = iopt(locacc+1) - 1
  inrows = iopt(locacc+2)
  mnew = mrows + inrows

  if ( mnew < 0 .or. mnew > mdw) then
      nerr = 10
      call xerrwv('dbols(). no. of rows=(i1) must be >= 0 .and. <= mdw=(i2).' &
                  ,nerr,level,2,mnew,mdw,0,rdum,rdum)
      go to 190
  end if

  do j = 1,min(ncols+1,mnew)
      do i = mnew,max(mrows,j) + 1,-1
          ibig = idamax(i-j,w(j,j),1) + j - 1
!
!  PIVOT FOR INCREASED STABILITY.
!
          call drotg(w(ibig,j),w(i,j),sc,ss)
          call drot(ncols+1-j,w(ibig,j+1),mdw,w(i,j+1),mdw,sc,ss)
          w(i,j) = 0.0D+00
    end do
  end do

  mrows = min(ncols+1,mnew)
  iopt(locacc+1) = mrows + 1
  igo = iopt(locacc)

  if ( igo == 2) then
      igo = 0
  end if
  go to 180

   90 continue
!
!  INITIALIZE VARIABLES AND DATA VALUES
!
  do 150 j = 1,ncols
      go to (100,110,120,130),iscale
      go to 140
  100     continue
!
!  THIS IS THE NOMINAL SCALING. EACH NONZERO
!  COL. HAS MAX. NORM EQUAL TO ONE.
!
      ibig = idamax(mrows,w(1,j),1)
      rw(j) = abs(w(ibig,j))
      if ( rw(j) == 0.0D+00 ) then
          rw(j) = one
      else
          rw(j) = one/rw(j)
      end if
      go to 140
  110     continue
!
!  THIS CHOICE OF SCALING MAKES EACH NONZERO COLUMN
!  HAVE EUCLIDEAN LENGTH EQUAL TO ONE.
!
      rw(j) = dnrm2(mrows,w(1,j),1)
      if ( rw(j) == 0.0D+00 ) then
          rw(j) = one
      else
          rw(j) = one/rw(j)
      end if
      go to 140
  120     continue
!
!  THIS CASE EFFECTIVELY SUPPRESSES SCALING BY SETTING
!  THE SCALING MATRIX TO THE IDENTITY MATRIX.
!
      rw(1) = one
      call dcopy(ncols,rw,0,rw,1)
      go to 160
  130     continue
      go to 160
  140     continue
  150 continue
  160 continue
!
!  SOLVE BOUNDED LEAST SQUARES PROBLEM
!
!  INITIALIZE IBASIS(*), J=1,NCOLS, AND IBB(*), J=1,NCOLS,
!  TO =J,AND =1, FOR USE IN DBOLSM( ).
!
  do j = 1,ncols
      iw(j) = j
      iw(j+ncols) = 1
      rw(3*ncols+j) = bl(j)
      rw(4*ncols+j) = bu(j)
  end do

  call dbolsm(w,mdw,mrows,ncols,rw(3*ncols+1),rw(4*ncols+1),ind, &
    iopt(lopt),x,rnorm,mode,rw(ncols+1),rw(2*ncols+1),rw,iw,iw(ncols+1))

  igo = 0
  180 continue
  return

  190 if ( mode>=0)mode = -nerr
  igo = 0
  return
end
subroutine dbolsm ( w, mdw, minput, ncols, bl, bu, ind, iopt, x, &
  rnorm, mode, rw, ww, scl, ibasis, ibb )

!***********************************************************************
!
!! DBOLSM solves E*X = F in the least squares sense with bounds on some X values.
!
!***BEGIN PROLOGUE  DBOLSM
!***DESCRIPTION
!
!          Solve E*X = F (least squares sense) with bounds on
!            selected X values.
!     The user must have dimension statements of the form:
!
!       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
!      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS)
!       integer IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS)
!
!     (here NX=number of extra locations required for options 1,...,7;
!     NX=0 for no options; here NI=number of extra locations possibly
!     required for options 1-7; NI=0 for no options; NI=14 if all the
!     options are simultaneously in use.)
!
!    INPUT
!    -----
!
!    --------------------
!    W(MDW,*),MROWS,NCOLS
!    --------------------
!     The array w(*,*) contains the matrix [E:F] on entry. The matrix
!     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in
!     the array W(*,*) with E occupying the first NCOLS columns and the
!     right side vector F in column NCOLS+1. The row dimension, MDW, of
!     the array W(*,*) must satisfy the inequality MDW >= MROWS.
!     Other values of MDW are errors. The values of MROWS and NCOLS
!     must be positive. Other values are errors.
!
!    ------------------
!    BL(*),BU(*),IND(*)
!    ------------------
!     These arrays contain the information about the bounds that the
!     solution values are to satisfy. The value of IND(J) tells the
!     type of bound and BL(J) and BU(J) give the explicit values for
!     the respective upper and lower bounds.
!
!    1.    For IND(J)=1, require X(J) >= BL(J).
!    2.    For IND(J)=2, require X(J) <= BU(J).
!    3.    For IND(J)=3, require X(J) >= BL(J) and
!                                X(J) <= BU(J).
!    4.    For IND(J)=4, no bounds on X(J) are required.
!     The values of BL(*),BL(*) are modified by the subprogram. Values
!     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3
!     (upper and lower bounds) the condition BL(J)  >  BU(J) is an
!     error.
!
!    -------
!    IOPT(*)
!    -------
!     This is the array where the user can specify nonstandard options
!     for DBOLSM( ). Most of the time this feature can be ignored by
!     setting the input value IOPT(1)=99. Occasionally users may have
!     needs that require use of the following subprogram options. For
!     details about how to use the options see below: IOPT(*) CONTENTS.
!
!     Option Number   Brief Statement of Purpose
!     ----- ------   ----- --------- -- -------
!           1         Move the IOPT(*) processing pointer.
!           2         Change rank determination tolerance.
!           3         Change blow-up factor that determines the
!                     size of variables being dropped from active
!                     status.
!           4         Reset the maximum number of iterations to use
!                     in solving the problem.
!           5         The data matrix is triangularized before the
!                     problem is solved whenever (NCOLS/MROWS)  <
!                     FAC. Change the value of FAC.
!           6         Redefine the weighting matrix used for
!                     linear independence checking.
!           7         Debug output is desired.
!          99         No more options to change.
!
!    ----
!    X(*)
!    ----
!     This array is used to pass data associated with options 1,2,3 and
!     5. Ignore this input parameter if none of these options are used.
!     Otherwise see below: IOPT(*) CONTENTS.
!
!    ----------------
!    IBASIS(*),IBB(*)
!    ----------------
!     These arrays must be initialized by the user. The values
!         IBASIS(J)=J, J=1,...,NCOLS
!         IBB(J)   =1, J=1,...,NCOLS
!     are appropriate except when using nonstandard features.
!
!    ------
!    SCL(*)
!    ------
!     This is the array of scaling factors to use on the columns of the
!     matrix E. These values must be defined by the user. To suppress
!     any column scaling set SCL(J)=1.0, J=1,...,NCOLS.
!
!    OUTPUT
!    ------
!
!    ----------
!    X(*),RNORM
!    ----------
!     The array X(*) contains a solution (if MODE >=0 or  == -22) for
!     the constrained least squares problem. The value RNORM is the
!     minimum residual vector length.
!
!    ----
!    MODE
!    ----
!     The sign of mode determines whether the subprogram has completed
!     normally, or encountered an error condition or abnormal status.
!     A value of MODE >= 0 signifies that the subprogram has completed
!     normally. The value of MODE (>= 0) is the number of variables
!     in an active status: not at a bound nor at the value ZERO, for
!     the case of free variables. A negative value of MODE will be one
!     of the 20 cases -40,-39,...,-22, or -1. Values  <  -1 correspond
!     to an abnormal completion of the subprogram. To understand the
!     abnormal completion codes see below: ERROR MESSAGES for DBOLSM( )
!     An approximate solution will be returned to the user only when
!     max. iterations is reached, MODE=-22.
!
!    -----------
!    RW(*),WW(*)
!    -----------
!     These are working arrays each with NCOLS entries. The array RW(*)
!     contains the working (scaled, nonactive) solution values. The
!     array WW(*) contains the working (scaled, active) gradient vector
!     values.
!
!    ----------------
!    IBASIS(*),IBB(*)
!    ----------------
!     These arrays contain information about the status of the solution
!     when MODE >= 0. The indices IBASIS(K), K=1,...,MODE, show the
!     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are
!     the active variables. The value (IBB(J)-1) is the number of times
!     variable J was reflected from its upper bound. (normally the user
!     can ignore these parameters.)
!
!    IOPT(*) CONTENTS
!    ------- --------
!     The option array allows a user to modify internal variables in
!     the subprogram without recompiling the source code. A central
!     goal of the initial software design was to do a good job for most
!     people. Thus the use of options will be restricted to a select
!     group of users. The processing of the option array proceeds as
!     follows: a pointer, here called LP, is initially set to the value
!     1. The value is updated as the options are processed.  At the
!     pointer position the option number is extracted and used for
!     locating other information that allows for options to be changed.
!     The portion of the array IOPT(*) that is used for each option is
!     fixed; the user and the subprogram both know how many locations
!     are needed for each option. A great deal of error checking is
!     done by the subprogram on the contents of the option array.
!     Nevertheless it is still possible to give the subprogram optional
!     input that is meaningless. For example some of the options use
!     the location X(NCOLS+IOFF) for passing data. The user must manage
!     the allocation of these locations when more than one piece of
!     option data is being passed to the subprogram.
!
!   1
!   -
!     Move the processing pointer (either forward or backward) to the
!     location IOPT(LP+1). The processing pointer is moved to location
!     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over
!     locations 3,...,NCOLS+2 of IOPT(*),
!
!       IOPT(1)=1
!       IOPT(2)=NCOLS+3
!       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
!       IOPT(NCOLS+3)=99
!       CALL DBOLSM( )
!
!     CAUTION: Misuse of this option can yield some very hard
!     -to-find bugs.  Use it with care.
!
!   2
!   -
!     The algorithm that solves the bounded least squares problem
!     iteratively drops columns from the active set. This has the
!     effect of joining a new column vector to the QR factorization of
!     the rectangular matrix consisting of the partially triangularized
!     nonactive columns. After triangularizing this matrix a test is
!     made on the size of the pivot element. The column vector is
!     rejected as dependent if the magnitude of the pivot element is
!     <= TOL* magnitude of the column in components strictly above
!     the pivot element. Nominally the value of this (rank) tolerance
!     is TOL = DRELPR, where DRELPR is relative machine
!     precision. To change only the value of TOL, for example,
!
!       X(NCOLS+1)=TOL
!       IOPT(1)=2
!       IOPT(2)=1
!       IOPT(3)=99
!       CALL DBOLSM()
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=TOL
!       IOPT(LP)=2
!       IOPT(LP+1)=IOFF
!        .
!       CALL DBOLSM()
!
!     The required length of IOPT(*) is increased by 2 if option 2 is
!     used; The required length of X(*) is increased by 1. A value of
!     IOFF <= 0 is an error. A value of TOL .le. DRELPR gives a
!     warning message; it is not considered an error.
!     Here DRELPR is the relative machine precision.
!
!   3
!   -
!     A solution component is left active (not used) if, roughly
!     speaking, it seems too large. Mathematically the new component is
!     left active if the magnitude is >=((vector norm of F)/(matrix
!     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(DRELPR)
!     where DRELPR is the relative machine precision. To change only
!     the value of BLOWUP, for example,
!
!       X(NCOLS+2)=BLOWUP
!       IOPT(1)=3
!       IOPT(2)=2
!       IOPT(3)=99
!       CALL DBOLSM()
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=BLOWUP
!       IOPT(LP)=3
!       IOPT(LP+1)=IOFF
!        .
!       CALL DBOLSM()
!
!     The required length of IOPT(*) is increased by 2 if option 3 is
!     used; the required length of X(*) is increased by 1. A value of
!     IOFF <= 0 is an error. A value of BLOWUP .le. 0.0 is an error.
!
!   4
!   -
!     Normally the algorithm for solving the bounded least squares
!     problem requires between NCOLS/3 and NCOLS drop-add steps to
!     converge. (this remark is based on examining a small number of
!     test cases.) The amount of arithmetic for such problems is
!     typically about twice that required for linear least squares if
!     there are no bounds and if plane rotations are used in the
!     solution method. Convergence of the algorithm, while
!     mathematically certain, can be much slower than indicated. To
!     avoid this potential but unlikely event ITMAX drop-add steps are
!     permitted. Nominally ITMAX=5*(MAX(MROWS,NCOLS)). To change the
!     value of ITMAX, for example,
!
!       IOPT(1)=4
!       IOPT(2)=ITMAX
!       IOPT(3)=99
!       CALL DBOLSM()
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       IOPT(LP)=4
!       IOPT(LP+1)=ITMAX
!        .
!       CALL DBOLSM()
!
!     The value of ITMAX must be  >  0. Other values are errors. Use
!     of this option increases the required length of IOPT(*) by 2.
!
!   5
!   -
!     For purposes of increased efficiency the MROWS by NCOLS+1 data
!     matrix [E:F] is triangularized as a first step whenever MROWS
!     satisfies FAC*MROWS  >  NCOLS. Nominally FAC=0.  To change the
!     value of FAC,
!
!       X(NCOLS+3)=FAC
!       IOPT(1)=5
!       IOPT(2)=3
!       IOPT(3)=99
!       CALL DBOLSM()
!
!     Generally, if LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=FAC
!       IOPT(LP)=5
!       IOPT(LP+1)=IOFF
!        .
!       CALL DBOLSM()
!
!     The value of FAC must be nonnegative. Other values are errors.
!     Using the value FAC=0.0 suppresses the initial triangularization step.
!     Use of this option increases the required length of IOPT(*) by 2;
!     The required length of of X(*) is increased by 1.
!
!   6
!   -
!     The norm used in testing the magnitudes of the pivot element
!     compared to the mass of the column above the pivot line can be
!     changed. The type of change that this option allows is to weight
!     the components with an index larger than MVAL by the parameter
!     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and
!     WT, where LP is the processing pointer for IOPT(*),
!
!       X(NCOLS+IOFF)=WT
!       IOPT(LP)=6
!       IOPT(LP+1)=IOFF
!       IOPT(LP+2)=MVAL
!
!     Use of this option increases the required length of IOPT(*) by 3.
!     The length of X(*) is increased by 1. Values of MVAL must be
!     nonnegative and not greater than MROWS. Other values are errors.
!     The value of WT must be positive. Any other value is an error. If
!     either error condition is present a message will be printed.
!
!   7
!   -
!     Debug output, showing the detailed add-drop steps for the
!     constrained least squares problem, is desired. This option is
!     intended to be used to locate suspected bugs.  To print,
!
!       IOPT(LP)=7
!
!   99
!   --
!     There are no more options to change.
!
!     The values for options are 1,...,7,99, and are the only ones
!     permitted. Other values are errors. Options -99,-1,...,-7 mean
!     that the repective options 99,1,...,7 are left at their default
!     values. An example is the option to modify the (rank) tolerance:
!
!       X(NCOLS+1)=TOL
!       IOPT(1)=-2
!       IOPT(2)=1
!       IOPT(3)=99
!
!***END PROLOGUE  DBOLSM
!
!     PURPOSE
!     -------
!     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED
!     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS:
!
!     SOLVE E*X =  F  (LEAST SQUARES SENSE)
!     WITH BOUNDS ON SELECTED X VALUES.
!
  implicit none

  integer mdw

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) big
  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) bou
  real ( kind = 8 ) bu(*)
  real ( kind = 8 ) cl1
  real ( kind = 8 ) cl2
  real ( kind = 8 ) cl3
  logical cnz
  real ( kind = 8 ) colabv
  real ( kind = 8 ) colblo
  logical constr
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) fac
  logical found
  integer i
  integer ibasis(*)
  integer ibb(*)
  integer icase
  integer idope
  integer idum
  integer igopr
  integer ind(*)
  integer inext
  integer ioff
  integer iopt(*)
  integer ip
  integer iprint
  integer itemp
  integer iter
  integer itmax
  integer j
  integer jbig
  integer jcol
  integer jdrop
  integer jdrop1
  integer jdrop2
  integer jp
  integer lds
  integer level
  integer lgopr
  integer lp
  integer minput
  integer mode
  integer mrows
  integer mval
  integer ncols
  integer nerr
  integer nlevel
  integer nsetb
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) rdum
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) rw(*)
  real ( kind = 8 ) sc
  real ( kind = 8 ) scl(*)
  real ( kind = 8 ) ss
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tolind
  real ( kind = 8 ) tolsze
  real ( kind = 8 ) w(mdw,*)
  real ( kind = 8 ) wlarge
  real ( kind = 8 ) wla
  real ( kind = 8 ) wlb
  real ( kind = 8 ) wt
  real ( kind = 8 ) ww(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xnew

  save

  common /dbocom/ idope(5)

  inext(idum) = min(idum+1,mrows)

  level = 1
!
!  VERIFY THAT THE PROBLEM DIMENSIONS ARE DEFINED PROPERLY.
!
  if ( minput<=0) then
      nerr = 31
      call xerrwv('dbolsm(). the number of rows=(i1) must be positive.' &
        ,nerr,level,1,minput,idum,0,rdum,rdum)
      go to 600
  end if

  if ( ncols<=0) then
      nerr = 32
      call xerrwv('dbolsm(). the number of cols.=(i1) must be positive.' &
        ,nerr,level,1,ncols,idum,0,rdum,rdum)
      go to 600
  end if

  if ( mdw < minput) then
      nerr = 33
      call xerrwv('dbolsm(). the row dimension of w(,)=(i1) must be >= the ' &
        //'number of rows=(i2).',nerr,level,2,mdw,mrows,0,rdum,rdum)
      go to 600
  end if
!
!  VERIFY THAT BOUND INFORMATION IS CORRECT.
!
  do j = 1,ncols
     if ( ind(j) < 1 .or. ind(j) > 4) then
         nerr = 34
         call xerrwv('dbolsm(). for j=(i1) the constraint indicator must be1-4.'&
                   ,nerr,level,2,j,ind(j),0,rdum,rdum)
         go to 600
     end if
  end do

  do j = 1,ncols
     if ( ind(j) == 3) then
         if ( bu(j) < bl(j)) then
             nerr = 35
             call xerrwv('dbolsm(). for j=(i1) the lower bound=(r1) is  ' // &
               '>  the upper bound=(r2).',nerr,level,1,j,idum,2,bl(j),bu(j))
             go to 600
         end if
     end if
  end do
!
!  CHECK THAT PERMUTATION AND POLARITY ARRAYS HAVE BEEN SET.
!
  do j = 1,ncols
     if ( ibasis(j) < 1 .or. ibasis(j) > ncols) then
         nerr = 36
         call xerrwv('dbolsm(). the input order of columns=(i1) is ' // &
           'not between 1 and ncols=(i2).',nerr,level,2,ibasis(j), &
           ncols,0,rdum,rdum)
         go to 600
     end if

     if ( ibb(j)<=0) then
         nerr = 37
         call xerrwv('dbolsm(). the bound polarity flag in component ' // &
           'j=(i1) must be positive. now=(i2).',nerr,level,2,j,ibb(j), &
           0,rdum,rdum)
         go to 600

     end if
  end do
!
!  PROCESS OPTION ARRAY
!
  go to 570

   40 continue
!
!  INITIALIZE VARIABLES AND DATA VALUES
!
  go to 460

   50 continue
  if ( iprint > 0) then
      call dmout(mrows,ncols+1,mdw,w,'('' pretri. input matrix'')',-4)
      call dvout(ncols,bl,'('' lower bounds'')',-4)
      call dvout(ncols,bu,'('' upper bounds'')',-4)
  end if

   60 continue
  iter = iter + 1
  if ( iter<=itmax) go to 80
  nerr = 22
  call xerrwv('dbolsm(). more than (i1)=itmax iterations solving bounded ' // &
    'least squares problem.', &
    nerr,level,1,itmax,idum,0,rdum,rdum)
!
!  RESCALE AND TRANSLATE VARIABLES
!
  igopr = 1
  go to 130

   70 continue

  go to 600

   80 continue
!
!  FIND A VARIABLE TO BECOME NON-ACTIVE
!
  go to 180

   90 continue
  if ( found) go to 110
!
!  RESCALE AND TRANSLATE VARIABLES
!
  igopr = 2
  go to 130

  100 continue
  mode = nsetb
  return

  110 continue
!
!  MAKE MOVE AND UPDATE FACTORIZATION
!
  go to 280

  120 continue
  go to 60
!
!  RESCALE AND TRANSLATE VARIABLES
!
  130 continue
  call dcopy(nsetb,x,1,rw,1)
  x(1) = 0.0D+00
  call dcopy(ncols,x,0,x,1)

  do j = 1,nsetb
     jcol = abs(ibasis(j))
     x(jcol) = rw(j)*abs(scl(jcol))
  end do

  do j = 1,ncols
     if ( mod(ibb(j),2) == 0) x(j) = bu(j) - x(j)
  end do

  do j = 1,ncols
     jcol = ibasis(j)
     if ( jcol < 0) x(-jcol) = bl(-jcol) + x(-jcol)
  end do

  do j = 1,ncols
     if ( scl(j) < 0.0D+00 ) x(j) = -x(j)
  end do

  call dscal(mrows-mval,wt,w(inext(mval),ncols+1),1)
  rnorm = dnrm2(mrows-max(nsetb,mval),w(inext(max(nsetb,mval)),ncols+1),1)

  go to (70,100),igopr
!
!  FIND A VARIABLE TO BECOME NON-ACTIVE
!
  180 continue
!
!  COMPUTE (NEGATIVE) OF GRADIENT VECTOR, W=(TRANSPOSE OF E)*(F-E*X).
!
  ww(1) = 0.0D+00
  call dcopy(ncols,ww,0,ww,1)

  do j = nsetb + 1,ncols
     jcol = abs(ibasis(j))
     ww(j) = ddot(mrows-nsetb,w(inext(nsetb),j),1, &
             w(inext(nsetb),ncols+1),1)*abs(scl(jcol))
  end do

  if ( iprint > 0) then
      call dvout(ncols,ww,'('' gradient values'')',-4)
      call ivout(ncols,ibasis,'('' internal variable order'')',-4)
      call ivout(ncols,ibb,'('' bound polarity'')',-4)
  end if

  200 continue
!
!  IF ACTIVE SET = NUMBER OF TOTAL ROWS, QUIT.
!
  if ( nsetb == mrows) then
      found = .false.
      go to 90
  end if
!
!  CHOOSE AN EXTREMAL COMPONENT OF GRADIENT VECTOR
!  FOR A CANDIDATE TO BECOME NON-ACTIVE.
!
  wlarge = -big
  jbig = 0
  cnz = .false.

  do 210 j = nsetb + 1,ncols

     t = ww(j)
!
!  SKIP LOOKING AT COMPONENTS FLAGGED AS NON-CANDIDATES.
!
     if ( t == big) go to 210
     itemp = ibasis(j)
     jcol = abs(itemp)
     if ( nsetb < mval) then
         cl1 = dnrm2(nsetb,w(1,j),1)
         cl2 = dnrm2(mval-nsetb,w(inext(nsetb),j),1)
         colabv = cl1
         colblo = cl2
     else
         cl1 = dnrm2(mval,w(1,j),1)
         cl2 = abs(wt)*dnrm2(nsetb-mval,w(inext(mval),j),1)
         cl3 = abs(wt)*dnrm2(mrows-nsetb,w(inext(nsetb),j),1)
         call drotg(cl1,cl2,sc,ss)
         colabv = abs(cl1)
         colblo = cl3
     end if

     if ( itemp < 0) then
         if ( mod(ibb(jcol),2) == 0) t = -t
!
!  SKIP LOOKING AT COMPONENTS THAT WOULD NOT DECREASE OBJECTIVE.
!
         if ( t < 0.0D+00 ) go to 210
     end if
!
!  THIS IS A COLUMN PIVOTING STEP THAT MAXIMIZES THE RATIO OF
!  COLUMN MASS ON AND BELOW THE PIVOT LINE RELATIVE TO THAT
!  STRICTLY ABOVE THE PIVOT LINE.
!
     if ( colabv == 0.0D+00 .and. .not. cnz) then
         t = colblo*abs(scl(jcol))
         if ( wlarge < t) then
             wlarge = t
             jbig = j
         end if
     else
         if ( .not. cnz) then
             wla = 0.0D+00
             wlb = 0.0D+00
             cnz = .true.
         end if

       if ( sqrt(colblo)*sqrt(wla) >= sqrt(colabv)*sqrt(wlb)) then
            wlb=colblo
            wla=colabv
            jbig=j
       end if
     end if

  210 continue
  if ( jbig == 0) then
      found = .false.
      if ( iprint > 0) then
        write ( *, '(a)' ) '  Found no variable to enter.'
      end if

      go to 90

  end if
!
!  SEE IF THE INCOMING COL. IS SUFFICIENTLY INDEPENDENT.
!  THIS TEST IS MADE BEFORE AN ELIMINATION IS PERFORMED.
!
  if ( iprint > 0) then
    write ( *, '(a,i6)' ) '  Try to bring in column ', jbig
  end if

  if ( cnz) then
  if ( wlb<=wla*tolind) then
      found = .false.
      if ( iprint > 0) then
        write ( *, '(a)' ) '  Variable is dependent, not used.'
      end if

      ww(jbig) = big
      go to 200

  end if
  end if
!
!  SWAP MATRIX COLS. NSETB+1 AND JBIG, PLUS POINTER INFO., AND
!  GRADIENT VALUES.
!
  nsetb = nsetb + 1
  if ( nsetb/=jbig) then
      call dswap(mrows,w(1,nsetb),1,w(1,jbig),1)
      call dswap(1,ww(nsetb),1,ww(jbig),1)
      itemp = ibasis(nsetb)
      ibasis(nsetb) = ibasis(jbig)
      ibasis(jbig) = itemp
  end if
!
!  ELIMINATE ENTRIES BELOW THE PIVOT LINE IN COL. NSETB.
!
  if ( mrows > nsetb) then

      do i = mrows,nsetb + 1,-1
         if ( i /= mval+1 ) then
           call drotg(w(i-1,nsetb),w(i,nsetb),sc,ss)
           w(i,nsetb) = 0.0D+00
           call drot(ncols-nsetb+1,w(i-1,nsetb+1),mdw,w(i,nsetb+1),mdw,sc,ss)
         end if
      end do

      if ( (mval>=nsetb) .and. (mval < mrows)) then
          t = w(nsetb,nsetb)
          if ( t/= 0.0D+00 ) then
              t = wt*w(mval+1,nsetb)/t
          else
              t = big
          end if

          if ( tolind*abs(t) > one) then
              call dswap(ncols-nsetb+2,w(nsetb,nsetb),mdw,w(mval+1,nsetb),mdw)
              call dscal(ncols-nsetb+2,wt,w(nsetb,nsetb),mdw)
              call dscal(ncols-nsetb+2,one/wt,w(mval+1,nsetb),mdw)
          end if

          call drotg(w(nsetb,nsetb),w(mval+1,nsetb),sc,ss)
          w(mval+1,nsetb) = 0.0D+00
          call drot(ncols-nsetb+1,w(nsetb,nsetb+1),mdw, &
                   w(mval+1,nsetb+1),mdw,sc,ss)
      end if

  end if

  if ( w(nsetb,nsetb) == 0.0D+00 ) then
      ww(nsetb) = big
      nsetb = nsetb - 1
      if ( iprint > 0) then
        write ( *, '(a)' ) '  Pivot is zero, not used.'
      end if

      go to 200
  end if
!
!  CHECK THAT NEW VARIABLE IS MOVING IN THE RIGHT DIRECTION.
!
  itemp = ibasis(nsetb)
  jcol = abs(itemp)
  xnew = (w(nsetb,ncols+1)/w(nsetb,nsetb))/abs(scl(jcol))

  if ( itemp < 0) then
      if ( ww(nsetb)>= 0.0D+00 .and. xnew<= 0.0D+00 ) go to 230
      if ( ww(nsetb)<= 0.0D+00 .and. xnew>= 0.0D+00 ) go to 230
  end if

  go to 240

  230 continue
  ww(nsetb) = big
  nsetb = nsetb - 1
  if ( iprint > 0) then
    write ( *, '(a)' ) '  Variable has bad direction, not used.'
  end if

  go to 200

  240 continue
  found = .true.
  go to 250

  250 continue

  go to 90
!
!  SOLVE THE TRIANGULAR SYSTEM
!
  260 continue
  call dcopy(nsetb,w(1,ncols+1),1,rw,1)

  do j = nsetb,1,-1
     rw(j) = rw(j)/w(j,j)
     jcol = abs(ibasis(j))
     t = rw(j)
     if ( mod(ibb(jcol),2) == 0) rw(j) = -rw(j)
     call daxpy(j-1,-t,w(1,j),1,rw,1)
     rw(j) = rw(j)/abs(scl(jcol))
  end do

  if ( iprint > 0) then
      call dvout(nsetb,rw,'('' soln. values'')',-4)
      call ivout(nsetb,ibasis,'('' cols. used'')',-4)
  end if

  go to (290,430),lgopr
!
!  MAKE MOVE AND UPDATE FACTORIZATION
!
  280 continue
!
!  SOLVE THE TRIANGULAR SYSTEM
!
  lgopr = 1
  go to 260

  290 continue
!
!  SEE IF THE UNCONSTRAINED SOL. (OBTAINED BY SOLVING THE
!  TRIANGULAR SYSTEM) SATISFIES THE PROBLEM BOUNDS.
!
  alpha = 2.0D+00
  beta = 2.0D+00
  x(nsetb) = 0.0D+00

  do j = 1,nsetb

     itemp = ibasis(j)
     jcol = abs(itemp)
     t1 = 2.0D+00
     t2 = 2.0D+00

     if ( itemp < 0) then
         bou = 0.0D+00
     else
         bou = bl(jcol)
     end if

     if ( (-bou)/=big) bou = bou/abs(scl(jcol))
     if ( rw(j)<=bou) t1 = (x(j)-bou)/ (x(j)-rw(j))
     bou = bu(jcol)
     if ( bou/=big) bou = bou/abs(scl(jcol))
     if ( rw(j)>=bou) t2 = (bou-x(j))/ (rw(j)-x(j))
!
!  IF NOT, THEN COMPUTE A STEP LENGTH SO THAT THE
!  VARIABLES REMAIN FEASIBLE.
!
     if ( t1 < alpha) then
         alpha = t1
         jdrop1 = j
     end if

     if ( t2 < beta) then
         beta = t2
         jdrop2 = j
     end if

  end do

  constr = alpha < 2.0D+00 .or. beta < 2.0D+00
  if ( constr) go to 310
!
!  ACCEPT THE CANDIDATE BECAUSE IT SATISFIES THE STATED BOUNDS
!  ON THE VARIABLES.
!
  call dcopy(nsetb,rw,1,x,1)
  go to 120

  310 continue
!
!  TAKE A STEP THAT IS AS LARGE AS POSSIBLE WITH ALL
!  VARIABLES REMAINING FEASIBLE.
!
  do j = 1,nsetb
     x(j) = x(j) + min(alpha,beta)* (rw(j)-x(j))
  end do

  if ( alpha<=beta) then
      jdrop2 = 0

  else
      jdrop1 = 0
  end if

  330 if ( jdrop1+jdrop2 > 0 .and. nsetb > 0) go to 340
  go to 450

  340 continue

  jdrop = jdrop1 + jdrop2
  itemp = ibasis(jdrop)
  jcol = abs(itemp)
  if ( jdrop2 > 0) then
!
!  VARIABLE IS AT AN UPPER BOUND.  SUBTRACT MULTIPLE OF THIS COL.
!  FROM RIGHT HAND SIDE.
!
      t = bu(jcol)
      if ( itemp > 0) then
          bu(jcol) = t - bl(jcol)
          bl(jcol) = -t
          itemp = -itemp
          scl(jcol) = -scl(jcol)

          do i = 1,jdrop
             w(i,jdrop) = -w(i,jdrop)
          end do

      else
          ibb(jcol) = ibb(jcol) + 1
          if ( mod(ibb(jcol),2) == 0) t = -t
      end if
!
!  VARIABLE IS AT A LOWER BOUND.
!
  else
      if ( itemp < 0.0D+00 ) then
          t = 0.0D+00
      else
          t = -bl(jcol)
          bu(jcol) = bu(jcol) + t
          itemp = -itemp
      end if

  end if

  call daxpy(jdrop,t,w(1,jdrop),1,w(1,ncols+1),1)
!
!  MOVE CERTAIN COLS. LEFT TO ACHIEVE UPPER HESSENBERG FORM.
!
  call dcopy(jdrop,w(1,jdrop),1,rw,1)

  do j = jdrop + 1,nsetb
     ibasis(j-1) = ibasis(j)
     x(j-1) = x(j)
     call dcopy(j,w(1,j),1,w(1,j-1),1)
  end do

  ibasis(nsetb) = itemp
  w(1,nsetb) = 0.0D+00
  call dcopy(mrows-jdrop,w(1,nsetb),0,w(jdrop+1,nsetb),1)
  call dcopy(jdrop,rw,1,w(1,nsetb),1)
!
!  TRANSFORM THE MATRIX FROM UPPER HESSENBERG FORM TO
!  UPPER TRIANGULAR FORM.
!
  nsetb = nsetb - 1

  do i = jdrop,nsetb
!
!  LOOK FOR SMALL PIVOTS AND AVOID MIXING WEIGHTED AND NONWEIGHTED ROWS.
!
     if ( i == mval) then
         t = 0.0D+00
         do j = i,nsetb
            jcol = abs(ibasis(j))
            t1 = abs(w(i,j)*scl(jcol))
            if ( t1 > t) then
                jbig = j
                t = t1
            end if
         end do
         go to 390
     end if

     call drotg(w(i,i),w(i+1,i),sc,ss)
     w(i+1,i) = 0.0D+00
     call drot(ncols-i+1,w(i,i+1),mdw,w(i+1,i+1),mdw,sc,ss)
  end do

  go to 420

  390 continue
!
!  THE TRIANGULARIZATION IS COMPLETED BY GIVING UP
!  THE HESSENBERG FORM AND TRIANGULARIZING A RECTANGULAR MATRIX.
!
  call dswap(mrows,w(1,i),1,w(1,jbig),1)
  call dswap(1,ww(i),1,ww(jbig),1)
  call dswap(1,x(i),1,x(jbig),1)
  itemp = ibasis(i)
  ibasis(i) = ibasis(jbig)
  ibasis(jbig) = itemp
  jbig = i
  do j = jbig,nsetb
     do i = j + 1,mrows
        call drotg(w(j,j),w(i,j),sc,ss)
        w(i,j) = 0.0D+00
        call drot(ncols-j+1,w(j,j+1),mdw,w(i,j+1),mdw,sc,ss)
     end do
  end do

  420 continue
!
!  SEE IF THE REMAINING COEFFICIENTS ARE FEASIBLE.  THEY SHOULD
!  BE BECAUSE OF THE WAY MIN(ALPHA,BETA) WAS CHOSEN.  ANY THAT ARE
!  NOT FEASIBLE WILL BE SET TO THEIR BOUNDS AND
!  APPROPRIATELY TRANSLATED.
!
  jdrop1 = 0
  jdrop2 = 0
!
!  SOLVE THE TRIANGULAR SYSTEM
!
  lgopr = 2
  go to 260

  430 continue
  call dcopy(nsetb,rw,1,x,1)

  do j = 1,nsetb

     itemp = ibasis(j)
     jcol = abs(itemp)

     if ( itemp < 0) then
         bou = 0.0D+00
     else
         bou = bl(jcol)
     end if

     if ( (-bou)/=big) bou = bou/abs(scl(jcol))

     if ( x(j)<=bou) then
         jdrop1 = j
         go to 330
     end if

     bou = bu(jcol)

     if ( bou/=big) bou = bou/abs(scl(jcol))

     if ( x(j)>=bou) then
         jdrop2 = j
         go to 330
     end if

  end do

  go to 330

  450 continue

  go to 120
!
!  INITIALIZE VARIABLES AND DATA VALUES
!
  460 continue
!
!  PRETRIANGULARIZE RECTANGULAR ARRAYS OF CERTAIN SIZES
!  FOR INCREASED EFFICIENCY.
!
  if ( fac*minput > ncols) then
      do j = 1,ncols + 1
         do i = minput,j + mval + 1,-1
            call drotg(w(i-1,j),w(i,j),sc,ss)
            w(i,j) = 0.0D+00
            call drot(ncols-j+1,w(i-1,j+1),mdw,w(i,j+1),mdw,sc,ss)
         end do
      end do
      mrows = ncols + mval + 1
  else
      mrows = minput
  end if
!
!  SET THE X(*) ARRAY TO ZERO SO ALL COMPONENTS ARE DEFINED.
!
  x(1) = 0.0D+00
  call dcopy(ncols,x,0,x,1)
!
!  THE ARRAYS IBASIS(*), IBB(*) ARE INITIALIZED BY THE CALLING
!  PROGRAM UNIT.
!  THE COL. SCALING IS DEFINED IN THE CALLING PROGRAM UNIT.
!  'BIG' IS PLUS INFINITY ON THIS MACHINE.
!
  big = huge ( big )

  do j = 1,ncols

     icase = ind(j)

     go to (490,500,510,520),icase

     go to 530

  490    bu(j) = big
     go to 530

  500    bl(j) = -big
     go to 530

  510    go to 530

  520    bl(j) = -big
     bu(j) = big

  530    continue

  end do

  do j = 1,ncols

     if ( (bl(j)<= 0.0D+00 .and. 0.0D+00 <= bu(j) .and. &
        abs(bu(j)) < abs(bl(j))) .or. bu(j) < 0.0D+00 ) then
         t = bu(j)
         bu(j) = -bl(j)
         bl(j) = -t
         scl(j) = -scl(j)
         do i = 1,mrows
            w(i,j) = -w(i,j)
         end do
     end if
!
!  INDICES IN SET T(=TIGHT) ARE DENOTED BY NEGATIVE VALUES OF IBASIS(*).
!
     if ( bl(j)>=0.0D+00 ) then
         ibasis(j) = -ibasis(j)
         t = -bl(j)
         bu(j) = bu(j) + t
         call daxpy(mrows,t,w(1,j),1,w(1,ncols+1),1)
     end if

  end do

  nsetb = 0
  iter = 0

  go to 50
!
!  PROCESS OPTION ARRAY
!
  570 continue
  if ( idope(5) == 1) then
      fac = x(ncols+idope(1))
      wt = x(ncols+idope(2))
      mval = idope(3)
  else
      fac = 0.0D+00
      wt = 1.0D+00
      mval = 0
  end if

  tolind = sqrt( epsilon ( tolind ) )
  tolsze = sqrt( epsilon ( tolsze ) )
  itmax = 5 * max ( minput, ncols )
  iprint = 0
!
!  CHANGES TO SOME PARAMETERS CAN OCCUR THROUGH THE OPTION
!  ARRAY, IOPT(*).  PROCESS THIS ARRAY LOOKING CAREFULLY
!  FOR INPUT DATA ERRORS.
!
  lp = 0
  lds = 0

  580 continue

  lp = lp + lds
!
!  TEST FOR NO MORE OPTIONS.
!
  ip = iopt(lp+1)
  jp = abs(ip)
  if ( ip == 99) then
      go to 590
  else if ( jp == 99) then
      lds = 1
      go to 580
  else if ( jp == 1) then
!
!  MOVE THE IOPT(*) PROCESSING POINTER.
!
      if ( ip > 0) then
          lp = iopt(lp+2) - 1
          lds = 0
      else
          lds = 2
      end if

      go to 580

  else if ( jp == 2) then
!
!  CHANGE TOLERANCE FOR RANK DETERMINATION.
!
      if ( ip > 0) then
          ioff = iopt(lp+2)
          if ( ioff<=0) then
              nerr = 24
              call xerrwv('dbolsm(). the offset=(i1) beyond postion ' // &
                'ncols=(i2) must be positive for option number 2.', &
                nerr,level,2,ioff,ncols,0,rdum,rdum)
              go to 600
          end if

          tolind = x(ncols+ioff)
          if ( tolind < epsilon ( tolind ) ) then
              nerr = 25
              nlevel = 0
              call xerrwv( 'dbolsm(). the tolerance for rank ' // &
                'determination=(r1) is less than machine precision=(r2).', &
                nerr,nlevel,0,idum,idum,2,tolind, epsilon ( tolind ) )
          end if
      end if

      lds = 2
      go to 580

  else if ( jp == 3) then
!
!  CHANGE BLOWUP FACTOR FOR ALLOWING VARIABLES TO BECOME INACTIVE.
!
      if ( ip > 0) then
          ioff = iopt(lp+2)
          if ( ioff<=0) then
              nerr = 26
              call xerrwv( 'dbolsm(). the offset=(i1) beyond position ' // &
                'ncols=(i2) must be postive for option number 3.', &
                nerr,level,2,ioff,ncols,0,rdum,rdum)
              go to 600
          end if

          tolsze = x(ncols+ioff)
          if ( tolsze<= 0.0D+00 ) then
              nerr = 27
              call xerrwv('dbolsm(). the reciprocal of the blow-up factor ' // &
                'for rejecting variables must be positive. now=(r1).', &
                nerr,level,0,idum,idum,1,tolsze,rdum)
              go to 600
          end if
      end if

      lds = 2
      go to 580

  else if ( jp == 4) then
!
!  CHANGE THE MAXimum Number OF ITERATIONS ALLOWED.
!
      if ( ip > 0) then
          itmax = iopt(lp+2)
          if ( itmax<=0) then
              nerr = 28
              call xerrwv('dbolsm(). the maximum number of iterations=(i1) ' // &
                'must be positive.',nerr,level,1,itmax,idum,0,rdum,rdum)
              go to 600
          end if
      end if

      lds = 2
      go to 580

  else if ( jp == 5) then
!
!  CHANGE THE FACTOR FOR PRETRIANGULARIZING THE DATA MATRIX.
!
      if ( ip > 0) then
          ioff = iopt(lp+2)
          if ( ioff<=0) then
              nerr = 29
              call xerrwv('dbolsm(). the offset=(i1) beyond position ' // &
                'ncols=(i2) must be postive for option number 5.', &
                nerr,level,2,ioff,ncols,0,rdum,rdum)
              go to 600
          end if

          fac = x(ncols+ioff)
          if ( fac < 0.0D+00 ) then
              nerr = 30
              nlevel = 0
              call xerrwv('dbolsm(). the factor (ncols/mrows) where ' // &
                'pre-triangularizing is performed must be nonnegative. ' // &
                'now=(r1).',nerr,nlevel,0,idum,idum,1,fac,rdum)
              go to 600
          end if
      end if

      lds = 2
      go to 580
  else if ( jp == 6) then
!
!  CHANGE THE WEIGHTING FACTOR (FROM ONE) TO APPLY TO COMPONENTS
!  NUMBERED .GT. MVAL (INITIALLY SET TO 1.)  THIS TRICK IS NEEDED
!  FOR APPLICATIONS OF THIS SUBPROGRAM TO THE HEAVILY WEIGHTED
!  LEAST SQUARES PROBLEM THAT COME FROM EQUALITY CONSTRAINTS.
!
      if ( ip > 0) then
        ioff = iopt(lp+2)
        mval = iopt(lp+3)
        wt = x(ncols+ioff)
      end if

      if ( mval < 0 .or. mval > minput .or. wt<= 0.0D+00 ) then
          nerr = 38
          nlevel = 0
          call xerrwv('dbolsm(). the row separator to apply weighting (i1) ' // &
            'must lie between 0 and mrows (i2). weight (r1) must be positive.', &
            nerr,nlevel,2,mval,minput,1,wt,rdum)
          go to 600
      end if

      lds = 3
      go to 580
!
!  TURN ON DEBUG OUTPUT.
!
  else if ( jp == 7) then
      if ( ip > 0) iprint = 1
      lds = 1
      go to 580
  else
      nerr = 23
      call xerrwv('dbolsm. the option number=(i1) is not defined.', &
           nerr,level,1,ip,idum,0,rdum,rdum)
      go to 600
  end if

  590 continue

  go to 40

  600 continue
  mode = -nerr
  return
end
subroutine dcopy ( n, x, incx, y, incy )

!*******************************************************************************
!
!! DCOPY copies one double precision vector into another.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be copied into Y.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) Y(*), the copy of X.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
  implicit none

  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = x(1:n)

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
function ddot ( n, x, incx, y, incy )

!*******************************************************************************
!
!! DDOT forms the dot product of two vectors.
!
!  Modified:
!
!    02 June 2000
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) X(*), one of the vectors to be multiplied.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input, real ( kind = 8 ) Y(*), one of the vectors to be multiplied.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
!    Output, real ( kind = 8 ) DDOT, the dot product of X and Y.
!
  implicit none

  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real ( kind = 8 ) ddot
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

    ddot = 0.0D+00

  else if ( incx == 1 .and. incy == 1 ) then

    ddot = dot_product ( x(1:n), y(1:n) )

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    stemp = 0.0D+00
    do i = 1, n
      stemp = stemp + x(ix) * y(iy)
      ix = ix + incx
      iy = iy + incy
    end do

    ddot = stemp

  end if

  return
end
subroutine dgeco ( a, lda, n, ipvt, rcond, z )

!*******************************************************************************
!
!! DGECO factors a double precision matrix and estimates its condition.
!
!  Discussion:
!
!     if  rcond  is not needed, dgefa is slightly faster.
!     to solve  a*x = b , follow dgeco by dgesl.
!     to compute  inverse(a)*c , follow dgeco by dgesl.
!     to compute  determinant(a) , follow dgeco by sgedi.
!     to compute  inverse(a) , follow dgeco by sgedi.
!
!  Parameters:
!
!     on entry
!
!        a       real(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   real
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond == 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       real(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ek
  integer info
  integer ipvt(n)
  integer j
  integer k
  integer kp1
  integer l
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Compute 1-norm of A.
!
  anorm = 0.0D+00
  do j = 1, n
     anorm = max ( anorm, dasum ( n, a(1,j), 1 ) )
  end do
!
!  Factor.
!
  call dgefa ( a, lda, n, ipvt, info )
!
!  RCOND = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!
!  estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!
!  trans(a)  is the transpose of a .  the components of  e  are
!  chosen to cause maximum local growth in the elements of w  where
!  trans(u)*w = e .  the vectors are frequently rescaled to avoid overflow.
!
!  Solve trans(u)*w = e.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do k = 1, n

     if ( z(k) /= 0.0D+00 ) then
       ek = sign ( ek, -z(k) )
     end if

     if ( abs ( ek - z(k) ) > abs ( a(k,k) ) ) then
        s = abs ( a(k,k) ) / abs ( ek - z(k) )
        z(1:n) = s * z(1:n)
        ek = s * ek
     end if

     wk = ek - z(k)
     wkm = -ek - z(k)
     s = abs ( wk )
     sm = abs ( wkm )

     if ( a(k,k) /= 0.0D+00 ) then
        wk = wk / a(k,k)
        wkm = wkm / a(k,k)
     else
        wk = 1.0D+00
        wkm = 1.0D+00
     end if

     kp1 = k + 1

     if ( k+1 <= n ) then

        do j = k+1, n
           sm = sm + abs ( z(j) + wkm * a(k,j) )
           z(j) = z(j) + wk * a(k,j)
           s = s + abs ( z(j) )
        end do

        if ( s < sm ) then
           t = wkm - wk
           wk = wkm
           do j = k+1, n
              z(j) = z(j) + t * a(k,j)
           end do

        end if

     end if

     z(k) = wk

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
!
!  Solve trans(l)*y = w.
!
  do k = n, 1, -1

     if ( k < n ) then
       z(k) = z(k) + dot_product ( a(k:n,k), z(k:n) )
     end if

     if ( abs ( z(k) ) > 1.0D+00 ) then
        s = 1.0D+00 / abs ( z(k) )
        z(1:n) = s * z(1:n)
     end if

     l = ipvt(k)
     t = z(l)
     z(l) = z(k)
     z(k) = t
  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = 1.0D+00
!
!  Solve l*v = y.
!
  do k = 1, n

     l = ipvt(k)
     t = z(l)
     z(l) = z(k)
     z(k) = t

     if ( k < n ) then
       call daxpy ( n-k, t, a(k+1,k), 1, z(k+1), 1 )
     end if

     if ( abs ( z(k) ) > 1.0D+00 ) then
        s = 1.0D+00 / abs ( z(k) )
        z(1:n) = s * z(1:n)
        ynorm = s * ynorm
     end if

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve u*z = v.
!
  do k = n, 1, -1

     if ( abs ( z(k) ) > abs ( a(k,k) ) ) then
        s = abs ( a(k,k) ) / abs ( z(k) )
        z(1:n) = s * z(1:n)
        ynorm = s * ynorm
     end if

     if ( a(k,k) /= 0.0D+00 ) then
       z(k) = z(k) / a(k,k)
     else
       z(k) = 1.0D+00
     end if

     t = -z(k)
     call daxpy ( k-1, t, a(1,k), 1, z(1), 1 )

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dgefa ( a, lda, n, ipvt, info )

!*******************************************************************************
!
!! DGEFA factors a double precision matrix.
!
!  Discussion:
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!  Parameters:
!
!     on entry
!
!        a       real(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) == 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or sgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  integer info
  integer ipvt(n)
  integer idamax
  integer j
  integer k
  integer kp1
  integer l
  real ( kind = 8 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1

     kp1 = k + 1
!
!  Find L = pivot index.
!
     l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
     ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
     if ( a(l,k) /= 0.0D+00 ) then
!
!  Interchange if necessary.
!
        if ( l /= k ) then
           t = a(l,k)
           a(l,k) = a(k,k)
           a(k,k) = t
        end if
!
!  Compute multipliers.
!
        t = -1.0D+00 / a(k,k)
        call dscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
        do j = k+1, n
           t = a(l,j)
           if ( l /= k ) then
              a(l,j) = a(k,j)
              a(k,j) = t
           end if
           call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
         end do

     else
        info = k
     end if

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgesl ( a, lda, n, ipvt, b, job )

!*******************************************************************************
!
!! DGESL solves a system factored by DGECO or DGEFA.
!
!  Discussion:
!
!    DGESL can solve either of the systems a * x = b  or  trans(a) * x = b
!    using the factors computed by dgeco or dgefa.
!
!  Parameters:
!
!     on entry
!
!        a       real(lda, n)
!                the output from dgeco or dgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!
!        b       real(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond > 0.0
!        or dgefa has set info == 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco ( a, lda, n, ipvt, rcond, z )
!           if (rcond is too small) go to ...
!           do j = 1, p
!              call dgesl ( a, lda, n, ipvt, c(1,j), 0 )
!           end do
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ipvt(n)
  integer job
  integer k
  integer l
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
!
!  JOB = 0, solve  a * x = b.
!
!  First solve  l*y = b.
!
  if ( job == 0 ) then

     do k = 1, n-1
        l = ipvt(k)
        t = b(l)
        if ( l /= k ) then
           b(l) = b(k)
           b(k) = t
        end if
        call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )
     end do
!
!  Now solve  u*x = y.
!
     do k = n, 1, -1
        b(k) = b(k) / a(k,k)
        t = -b(k)
        call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
     end do
!
!  JOB = nonzero, solve  trans(a) * x = b.
!  First solve  trans(u)*y = b
!
  else

     do k = 1, n
        t = dot_product ( a(1:k-1,k), b(1:k-1) )
        b(k) = (b(k) - t)/a(k,k)
     end do
!
!  Now solve trans(l)*x = y.
!
     do k = n-1, 1, -1
        b(k) = b(k) + ddot ( n-k, a(k+1,k), 1, b(k+1), 1 )
        l = ipvt(k)
        if ( l /= k ) then
           t = b(l)
           b(l) = b(k)
           b(k) = t
        end if
     end do

  end if

  return
end
subroutine dmout ( m, n, lda, a, ifmt, idigit )

!***********************************************************************
!
!! DMOUT prints double precision matrices.
!
!  EXAMPLE..
!
!  PRINT AN ARRAY CALLED (SIMPLEX TABLEAU   ) OF SIZE 10 BY 20 SHOWING
!  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
!  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
!
!     double precision TABLEU(20,20)
!     M = 10
!     N = 20
!     LDTABL = 20
!     IDIGIT = -6
!     CALL DMOUT(M,N,LDTABL,TABLEU,'(''1SIMPLEX TABLEAU'')',IDIGIT)
!
!  Author:
!
!    JOHN A. WISNIEWSKI and RICHARD J. HANSON,
!    SANDIA LABS ALBUQUERQUE.
!
!  INPUT..
!
!  M,N,LDA,A(*,*) PRINT THE double precision ARRAY A(I,J),I=1,...,M,
!                 J=1,...,N, ON OUTPUT UNIT *. LDA IS THE DECLARED
!                 FIRST DIMENSION OF A(*,*) AS SPECIFIED IN THE CALLING
!                 PROGRAM. THE HEADING IN THE FORTRAN FORMAT STATEMENT
!                 IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST STEP.
!                 THE COMPONENTS A(I,J) ARE INDEXED, ON OUTPUT, IN A
!                 PLEASANT FORMAT.
!  IFMT(*)        A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON
!                 OUTPUT UNIT * WITH THE VARIABLE FORMAT FORTRAN
!                 STATEMENT
!                       write(*,IFMT).
!  IDIGIT         PRINT AT LEAST IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
!                 THE SUBPROGRAM WILL CHOOSE THAT integer 6,14,20 OR 28
!                 WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF
!                 PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE
!                 UTILIZED TO WRITE EACH LINE OF OUTPUT OF THE ARRAY
!                 A(*,*). (THIS CAN BE USED ON MOST TIME-SHARING
!                 TERMINALS).  IF IDIGIT.GE.0, 133 PRINTING COLUMNS ARE
!                 UTILIZED. (THIS CAN BE USED ON MOST LINE PRINTERS).
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  integer i
  integer idigit
  character ( len = * ) ifmt
  character icol(3)
  integer j
  integer k1
  integer k2
  integer m
  integer ndigit

  data icol(1),icol(2),icol(3)/'c','o','l'/

  write ( *, ifmt )

  if ( m <= 0 .or. n <= 0 .or. lda <= 0 ) then
    return
  end if

  ndigit = idigit
  if ( idigit == 0) ndigit = 6
  if ( idigit>=0) go to 80

  ndigit = -idigit
  if ( ndigit > 6) go to 20

  do k1=1,n,4
    k2 = min0(n,k1+3)
    write(*,1000) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1004) i,(a(i,j),j=k1,k2)
    end do
  end do

  return

   20 continue
  if ( ndigit > 14) go to 40

  do k1=1,n,2
    k2 = min ( n, k1+1 )
    write(*,1001) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1005) i,(a(i,j),j=k1,k2)
    end do
  end do

  return

   40 continue
  if ( ndigit > 20) go to 60

  do k1=1,n,2
    k2=min(n,k1+1)
    write(*,1002) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1006) i,(a(i,j),j=k1,k2)
    end do
  end do
  return

   60 continue
  do 70 k1=1,n
  k2 = k1
  write(*,1003) (icol,i,i=k1,k2)
  do 70 i=1,m
  write(*,1007) i,(a(i,j),j=k1,k2)
   70 continue
  return

   80 continue
  if ( ndigit > 6) go to 100

  do k1=1,n,8
    k2 = min0(n,k1+7)
    write(*,1000) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1004) i,(a(i,j),j=k1,k2)
    end do
  end do
  return

  100 continue
  if ( ndigit > 14) go to 120

  do k1=1,n,5
    k2 = min(n,k1+4)
    write(*,1001) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1005) i,(a(i,j),j=k1,k2)
    end do
  end do

  return

  120 continue
  if ( ndigit > 20) go to 140

  do k1=1,n,4
    k2 = min(n,k1+3)
    write(*,1002) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1006) i,(a(i,j),j=k1,k2)
    end do
  end do

  return

  140 continue

  do k1=1,n,3
    k2 = min0(n,k1+2)
    write(*,1003) (icol,i,i=k1,k2)
    do i=1,m
      write(*,1007) i,(a(i,j),j=k1,k2)
    end do
  end do

  return
 1000 format(10x,8(5x,3a1,i4,2x))
 1001 format(10x,5(9x,3a1,i4,6x))
 1002 format(10x,4(12x,3a1,i4,9x))
 1003 format(10x,3(16x,3a1,i4,13x))
 1004 format(1x,'row',i4,2x,1p8d14.5)
 1005 format(1x,'row',i4,2x,1p5d22.13)
 1006 format(1x,'row',i4,2x,1p4d28.19)
 1007 format(1x,'row',i4,2x,1p3d36.27)
end
function dnrm2 ( n, x, incx )

!*******************************************************************************
!
!! DNRM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    The original DNRM2 algorithm is accurate but written in a bizarre,
!    unreadable and obsolete format.  This version goes for clarity.
!
!  Modified:
!
!    01 June 2000
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!
  implicit none

  integer i
  integer incx
  integer ix
  integer n
  real ( kind = 8 ) damax
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xmax

  if ( n <= 0 ) then

    dnrm2 = 0.0D+00

  else

    xmax = damax ( n, x, incx )

    if ( xmax == 0.0D+00 ) then

      dnrm2 = 0.0D+00

    else

      if ( incx >= 0 ) then
        ix = 1
      else
        ix = ( - n + 1 ) * incx + 1
      end if

      stemp = 0.0D+00
      do i = 1, n
        stemp = stemp + ( x(ix) / xmax )**2
        ix = ix + incx
      end do

      dnrm2 = xmax * sqrt ( stemp )

    end if

  end if

  return
end
subroutine dpchek ( df, dqedev, fj, iopt, ldfj, nvars, ropt, x, y )

!***********************************************************************
!
!! DPCHEK checks the user's jacobian routine.
!
!  Modified:
!
!    11 September 2002
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) DF(NVARS).
!
!    Input, external DQEDEV, the name of the user written jacobian
!    and function evaluation routine.
!
!    Workspace, real ( kind = 8 ) FJ(LDFJ,NVARS+1), space to store
!    the jacobian and the function, as required by DQEDEV.
!
!    Throughput, integer IOPT(*), parameters to be passed to DQEDEV.
!
!    Input, integer LDFJ, the leading dimension of FJ, which must
!    be at least NVARS.
!
!    Input, integer NVARS, the number of variables.
!
!    Throughput, real ( kind = 8 ) ROPT(*), parameters to be passed to DQEDEV.
!
!    Input, real ( kind = 8 ) X(NVARS), the point at which the
!    jacobian should be evaluated.
!
!    Workspace, real ( kind = 8 ) Y(NVARS).
!
  implicit none

  integer ldfj
  integer nvars

  real ( kind = 8 ) df(nvars)
  external dqedev
  real ( kind = 8 ) eps
  real ( kind = 8 ) fj(ldfj,nvars+1)
  integer igo
  integer iopt(*)
  integer j
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) t
  real ( kind = 8 ) x(nvars)
  real ( kind = 8 ) y(nvars)
!
!  Get the square root of the machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  Consider each component X(J) of the set of variables.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DPCHEK:'
  write ( *, '(a)' ) '  Compare user jacobian and function for'
  write ( *, '(a)' ) '  consistency, using finite differences.'
  write ( *, '(a)' ) ' '

  call dvout ( nvars, x, '(///,'' Evaluation point X'')', -4 )

  do j = 1, nvars
!
!  Set the appropriate increment T to X(J).
!
    t = eps * ( abs ( x(j) ) + 1.0D+00 )
!
!  Make a copy YP of X, with Y(J) incremented by T.
!
    y(1:nvars) = x(1:nvars)
    y(j) = x(j) + t
!
!  Evaluate F(YP).
!
    igo = 0
    call dqedev ( y, fj, ldfj, igo, iopt, ropt )
!
!  Save F(YP).
!
    df(1:nvars) = fj(1:nvars,nvars+1)
!
!  Make a copy YM of X, with Y(J) decremented by T.
!
    y(j) = x(j) - t
!
!  Evaluate F(YM).
!
    igo = 0
    call dqedev ( y, fj, ldfj, igo, iopt, ropt )
!
!  Estimate the partial derivative d F/d X(J) by (F(YP)-F(YM))/2*T
!
    df(1:nvars) = ( df(1:nvars) - fj(1:nvars,nvars+1) ) / ( 2.0D+00 * t )
!
!  Evaluate the user's formula for the partial derivatives.
!
    igo = 1
    call dqedev ( x, fj, ldfj, igo, iopt, ropt )

    call dvout(nvars,df,'(/,'' Numerical derivative'')',-4)

    write ( *, '(a,i6)' ) '  Variable number ', j

    call dvout(nvars,fj(1,j),'('' Analytic partial'')',-4)

  end do

  return
end
subroutine dqed ( dqedev, mequa, nvars, mcon, ind, bl, bu, x, fjac, &
  ldfjac, fnorm, igo, iopt, ropt, iwa, wa )

!***********************************************************************
!
!! DQED solves bounded and constrained least squares and nonlinear equations.
!
!  Modified:
!
!    11 September 2002
!
!  Author:
!
!    Richard Hanson
!    Fred Krogh
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide, 
!    SIAM (Society for Industrial and Applied Mathematics),
!    Philadelphia, 1979.
!
!    Richard Hanson,
!    Least Squares with Bounds and Linear Constraints,
!    SIAM Journal of Scientific and Statistical Computing, 
!    Volume 7, Number 3, July 1986, pages 826-834.
!
!    Robert Schnabel and P Frank,
!    Tensor Methods for Nonlinear Equations,
!    SIAM Journal on Numerical Analysis, 
!    Volume 21, Number 5, October 1984, pages 815-843.
! 
!***PURPOSE  SOLVE NONLINEAR LEAST SQUARES AND NONLINEAR
!            EQUATIONS.  USER PROVIDES SIMPLE BOUNDS, LINEAR
!            CONSTRAINTS AND EVALUATION CODE FOR THE FUNCTIONS.
!***LONG DESCRIPTION
!        SUBROUTINE DQED (DQEDEV, MEQUA, NVARS, MCON, IND, BL, BU, X,
!       *            FJ, LDFJ, RNORM, IGO, IOPT, ROPT,
!       *            IWORK, WORK)
!
!
!  Table of Sections
!  -----------------
!  1. Introduction
!     ------------
!  2. Calling Sequence Explained
!     ------- -------- ---------
!  3. Remarks on the Usage Examples
!     ------- -- --- ----- --------
!  5. References
!     ----------
!
!  1. Introduction
!     ------------
!  This software package is available in both single and double
!  precision.  The double precision version  (type  REAL*8)  is
!  described  below.   For the REAL  version  of the
!  documentation  substitute  'REAL' for 'DOUBLE  PRECISION' in the
!  type  statements.  Change the names of the subprograms: 'DQED()'
!  to 'SQED()', 'DQEDEV()' to 'SQEDEV()', and 'D1MACH' to 'R1MACH.'
!
!  The Fortran subprogram, DQED(), solves the constrained nonlinear
!  least squares problem:
!
!    Minimize  the  sum  of  squares  of  MEQUA (generally nonlinear)
!    equations,
!
!       f (x) = 0, I=1,...,MEQUA                    Eq. (1)
!        I
!
!  where x is  a set  of  NVARS unknowns.  (The vector
!  function with  these  MEQUA  components  is  called f(x) in the
!  discussion that  follows.)   The components of x may have upper
!  and lower bounds  given  by  the  user.   (In  fact all of the
!  possible cases, no bounds, bounds at one end only, or upper and
!  lower bounds  can  be  specified.)   Linear  constraints on the
!  unknowns, more  general than simple bounds,  can also be given.
!  These constraints can be of the equality or inequality type:
!
!       a  x + ... + a       x      =  y , L = 1,...,MCON,
!        L1 1         L,NVARS NVARS     L
!                                                   Eq. (2)
!
!  with bounds specified on the y , again given by the user.  The
!                                L
!  constraints can actually be slightly nonlinear.  In this case
!  the constraints can be described as:
!
!       g (x) =  y , L = 1,...,MCON,                Eq. (2')
!        L        L
!  where bounds are specified on each y .  The functions g (x) must
!                                      L                  L
!  be  defined for all x in the set described by the simple bounds.
!  Experienced users may wish to turn directly to Examples 1 and 2,
!  listed  below,  before  reading  the  subprogram  documentation.
!  There  is  no  size relation required for the problem dimensions
!  MEQUA,  NVARS,  and  MCON  except  that MEQUA and NVARS are both
!  positive,  and   MCON is nonnegative.
!
!  This code package will do a decent job of solving most nonlinear
!  least squares problems that can be expressed as Eqs. (1) and (2)
!  above, provided  that  continuous  derivatives of the functions
!  with respect  to the parameters can be computed.  This can also
!  include problems  where  the derivatives must be computed using
!  some  form  of numerical    differentiation.    Numerical
!  differentiation is not provided with this software for solving
!  nonlinear least squares problems.  Refer to the subprogram
!  JACG for numerical differentiation.  (Note: D. Salane has this
!  submitted to TOMS.  It is not included here.)
!
!  The authors also  plan  to develop methods that will do a much
!  better job of  coping  with  constraints more general than the
!  essentially linear ones indicated above in Eqs. (2)-(2').  There
!  are nonlinear  least squares problems with innocent looking but
!  highly nonlinear  constraints  where  this package will fail to
!  work.   The authors also hope to reduce the overhead required by
!  the software.  This high overhead is due primarily to the method
!  used  to  solve  the  inner-loop  quadratic  model problem.  The
!  authors  recommend  that  users consider using the option number
!  14, described below, to suppress use of the quadratic model. The
!  user  may  find  that  the software works quite well without the
!  quadratic  model.  This  may  be important when the function and
!  derivatives  evaluations  are  not expensive but many individual
!  problems are being solved.
!
!  There are two fundamental  ways to use the subprogram DQED().
!  The most  staightforward way is to make one Fortran CALL to the
!  subprogram and  obtain  values  for  the unknowns, x.  The user
!  provides a subprogram DQEDEV(), described below, that gives the
!  subprogram DQED() values of the functions f(x) and g(x), and the
!  derivative or  Jacobian  matrices  for  f(x)  and  g(x) at each
!  desired point x.  This usage is called 'forward communication.'
!  An alternate  way to use the subprogram is to provide an option
!  that allows  the  user  to communicate these values by 'reverse
!  communication.'   The  subprogram returns to the calling program
!  unit and  requests  values  for f(x) and g(x), and the Jacobian
!  matrices  for  f(x)  and  g(x)  for  a  given value of x.  (This
!  framework   is   often   required   in  applications  that  have
!  complicated  algorithmic  requirements  for  evaluation  of  the
!  functions.)   An  example  using  both  'forward'  and 'reverse'
!  communication  is  provided  below  (see  Remarks  on  the Usage
!  Examples) for least squares fitting of two exponential functions
!  to five data points.
!
!  2. Calling Sequence Explained
!     ------- -------- ---------
!  There   are  arrays  used  by  the  subprogram  that  must  have
!  dimensions equivalent to the following declarations.
!
!        integer MEQUA, NVARS, MCON, LDFJ, IGO
!        integer IND(NVARS+MCON), IOPT(LIOPT), IWORK(LIWORK)
!
!        double precision BL(NVARS+MCON), BU(NVARS+MCON), X(NVARS), RNORM,
!       *ROPT(LROPT), FJ(LDFJ,NVARS+1), WORK(LWORK)
!
!        EXTERNAL DQEDEV
!
!  The array dimensions must satisfy the bounds:
!
!        LIOPT >=  Number required for options in use.
!        LROPT >= Number required for options in use.
!         LDFJ >= MEQUA+MCON,
!
!  The  array  dimensions  for  the arrays IWORK(*) and WORK(*) can
!  change  if either option 14 or option 15 are in use.  For use in
!  the formulas, define:
!
!       MC=MCON
!       ME=MEQUA
!       NV=NVARS
!       MX=MAX(MEQUA,NVARS)
!
!  If the user is not using option 15, then
!
!       NT=5.
!
!  If the user is using option 15, then
!
!       NT=new number, must be >= 2.
!
!  If the user is not using option 14, then
!
!       NA=MC+2*NV+NT.
!
!  If the user is using option 14, then
!
!       NA=MC+NV+1.
!
!
!  In terms of these values defined above,
!        LIWORK >= 3*MC+9*NV+4*NT+NA+10
!         LWORK >= NA*(NA+4)+NV*(NT+33)+(ME+MX+14)*NT+9*MC+26
!
!  The  subprogram  DQEDEV  must  be declared in a Fortran EXTERNAL
!  statement:
!
!        EXTERNAL DQEDEV
!
!  Initialize the values of the parameters:
!
!        MEQUA, NVARS, MCON, IND(*), BL(*), BU(*), X(*), LDFJ,
!        IOPT(*), IWORK(1), IWORK(2),
!
!        CALL DQED  (DQEDEV, MEQUA, NVARS, MCON, IND, BL, BU, X,
!       *            FJ, LDFJ, RNORM, IGO, IOPT, ROPT,
!       *            IWORK, WORK)
!
!  Subprogram parameters:
!
!  DQEDEV (Input)
!  -----
!  This is  the  name  of  a subprogram that the user will usually
!  supply for  evaluation  of  the  values  of the constraints and
!  model, and  the  derivatives  of these functions. The user must
!  provide this subprogram unless 'reverse communication' is used.
!  A  model for  writing the subprogram DQEDEV() is provided under
!  the heading Example 1 Using Forward Communication, listed below.
!  Users  may  find  it  convenient to modify this model subprogram
!  when  writing  a  subprogram  for  their  own  application.   If
!  'reverse communication' is used, the user does not need to write
!  a stub or dummy subroutine named DQEDEV().  All that is required
!  is  to  declare exactly this name in an EXTERNAL statement.  The
!  code  package  has a dummy subroutine DQEDEV() that will be used
!  in   the   linking  or  load  step.   Example  2  Using  Reverse
!  Communication, listed below, illustrates this detail.
!
!  MEQUA, NVARS, MCON (Input)
!  ------------------
!  Respectively  they  are:  The number of least squares equations,
!  the  number  of unknowns or variables, and the number of general
!  constraints  for the solution, not including simple bounds.  The
!  values  of  MEQUA  and NVARS must be positive; the value of MCON
!  must  be  nonnegative.   Other  values  for these parameters are
!  errors.
!
!  IND(*),BL(*),BU(*) (Input)
!  ------------------
!  These  arrays  describe  the  form of the simple bounds that the
!  components of x are to satisfy.  Components numbered 1,...,NVARS
!  are  used  to  describe  the  form of the simple bounds that the
!  unknown     are     to     satisfy.      Components     numbered
!  NVARS+1,...,NVARS+MCON  are  used  to  describe  the form of the
!  general  MCON linear constraints.  The first NVARS components of
!  IND(*)  indicate  the type of simple bounds that the solution is
!  to  satisfy.   The  corresponding entries of BL(*) and BU(*) are
!  the bounding value.  The only values of IND(*) allowed are 1,2,3
!  or 4.  Other values are errors.  Specifically:
!
!  IND(J)=1, if x >= BL(J) is required; BU(J) is not used.
!                J
!        =2, if x <= BU(J) is required; BL(J) is not used.
!                J
!        =3, if x >= BL(J) and x <= BU(J) is required.
!                J                J
!        =4, if no bounds on x  are required;
!                             J
!                BL(*),BU(*) are not used.
!  General  linear constraints of the form shown in Eq. (2) require
!  that bounds be given for the linear functions y .  Specifically:
!                                                 L
!
!  IND(NVARS+L)=1,  if y >= BL(NVARS+L) is required; BU(*) is not
!                       L
!                 needed.
!
!              =2, if y <= BU(NVARS+L) is required; BL(*) is not
!                      L
!                  needed.
!              =3, if y >= BL(NVARS+L) and y <= BU(NVARS+L)
!                      L                      L
!
!              =4, if no bounds on y  are required;
!                                   L
!                  BL(*),BU(*) are not used.
!
!  The  values of the bounds for the unknowns may be changed by the
!  user  during  the  evaluation of the functions f(x) and g(x) and
!  their Jacobian matrices.
!
!  X(*),FJ(*,*),LDFJ (Input and Output, except LDFJ which is Input)
!  -----------------
!  The  array  X(*)  contains  the  NVARS  values,  x,   where  the
!  functions  f(x)  and  g(x)  and  their Jacobian matrices will be
!  evaluated  by  the subprogram DQED().  After the computation has
!  successfully  completed, the array X(*) will contain a solution,
!  namely  the  unknowns of the problem, x.  (Success is determined
!  by  an  appropriate  value for IGO.  This parameter is described
!  below.)  Initially  the array X(*) must contain a starting guess
!  for  the  unknowns of the problem, x.  The initial values do not
!  need  to  satisfy the constraints or the bounds.  If they do not
!  satisfy the bounds, then the point will be simply projected onto
!  the  bounds  as  a  first  step.  The first linear change to the
!  values  of  x must satisfy the general constraints.  (It is here
!  that the assumption of their linearity is utilized.)
!
!  The  Fortran  two-dimensional array FJ(*,*) is used to store the
!  linear  constraints  of  Eq. (2) (or more generally the Jacobian
!  matrix  of  the functions g(x) with respect to the variables x),
!  and  the  Jacobian  matrix  of  the function f(x).  The Jacobian
!  matrix of the (linear) constraints is placed in rows 1,...,MCON.
!  The  Jacobian  matrix  of  f(x)  is  placed in rows MCON+1, ...,
!  MCON+MEQUA.   The parameter LDFJ is the leading or row dimension
!  of  the  array  FJ(*,*).  Normally the array FJ(*,*) is assigned
!  values   by   the   user  when  the  nonlinear  solver  requests
!  evaluations  of  the  constraints  g(x)  and  thefunction f(x)
!  together  with  the Jacobian matrices G(x) and J(x).  The values
!  of  the  constraintfunctions  g (x)  are  placed  in the array
!                                   L
!  FJ(L,NVARS+1),  L=1,...,MCON.  The values of the model functions
!  f (x)  are  placed  in  the array at entries FJ(MCON+I,NVARS+1),
!   I
!  I=1,...,MEQUA.   Note  that the second dimension of FJ(*,*) must
!  be at least NVARS+1 in value.
!
!  RNORM (Output)
!  -----
!  This is the value of the Euclidean length or square root of sums
!  of  squares  of  components  of  thefunction  f(x)  after  the
!  approximate solution, x, has been found.  During the computation
!  it  is  updated  and equals the best value of the length of f(x)
!  that has been found.
!
!  IGO (Output; it can be an Input if interrupting the code)
!  ---
!  This  flag  directs  user  action and informs the user about the
!  type  of  results obtained by the subprogram.  The user may find
!  it  convenient  to  treat  the cases abs(IGO) <= 1 the same as
!  IGO=1.  This has no effect on the solution process.
!
!  The  user  can  interrupt the computation at any time, obtaining
!  the  best  values  of the vector x up to this point,  by setting
!  IGO  to any value  >  1 and then return control to DQED().  For
!  example  if  a  calculation  must be done in a certain length of
!  time,  the  user  can,  as  the end of the time draws near, set
!  IGO=20 and return control to DQED().  It is important that this
!  method be  used  to  stop further computing, rather than simply
!  proceeding.  The reason for this is that certain flags in DQED()
!  must be reset before any further solving on subsequent problems
!  can take  place.  The value of IGO  >  1 used to interrupt the
!  computation  is  arbitrary and is the value of IGO returned.  If
!  values of IGO =2,...,18 are used to flag this interrupt, they do
!  not mean the same thing as indicated below.  For this reason the
!  value  IGO=20 is recommended for signaling interrupts in DQED().
!  Another situation  that  may  occur  is  the  request  for  an
!  evaluation of  the functions and derivatives at a point x where
!  these can't be evaluated.  If this occurs, set IGO=99 and return
!  control to  DQED().   This will have the effect of defining the
!  derivatives to  be  all  zero  and the functions to be 'large.'
!  Thus a  reduction  in  the trust region around the current best
!  estimate  will occur.  Assigning the value IGO=99 will not cause
!  DQED() to stop computing.
!
!      =0   Place  the  value  of  f(x)  in FJ(MCON+*,NVARS+1).  If
!  'reverse  communication'  is  being used, CALL DQED() again.  If
!  'forward communication' is being used, do a RETURN.
!
!      =1  or  (-1)   Evaluate the Jacobians for the functions g(x)
!  and  f(x) as well as evaluating g(x) and f(x).  Use the vector x
!  that is now in the array X(*) as the values  where this
!  evaluation  will  be performed.  Place the Jacobian matrix
!  for  g(x) in the first MCON rows of FJ(*,*).  Place the Jacobian
!  matrix for f(x) in rows MCON+1,...,MCON+MEQUA in FJ(*,*).  Place
!  the  value of g(x) in FJ(*,NVARS+1).  Place the value of f(x) in
!  FJ(MCON+*,NVARS+1).
!
!  (Users who have complicated functions whose derivatives cannot be
!  computed analytically may want to use the numerical differentiation
!  subroutine JAGC.  This is available on the SLATEC library.)
!
!  If  'reverse communication' is being used, CALL DQED() again.
!  If 'forward communication' is being used, do a RETURN.
!
!  A  value  IGO=(-1)  flags  that  that the number of terms in the
!  quadratic  model  is  being  restricted by the amount of storage
!  given  for  that  purpose.   It  is  suggested,  but  it  is not
!  required,  that  additional  storage  be given for the quadratic
!  model  parameters.   See the following description of The Option
!  Array,  option  number 15, for the way to designate more storage
!  for this purpose.
!
!      =2   The function f(x) has a length less than TOLF.  This is
!  the  value  for  IGO to be expected when an actual zero value of
!  f(x)  is  anticipated.   See the description of The Option Array
!  for the value.
!
!      =3   Thefunction  f(x)  has  reached a value that may be a
!  local   minimum.   However,  the  bounds  on  the  trust  region
!  defining  the size of the step are being hit at each step.  Thus
!  the  situation  is  suspect.  (Situations of this type can occur
!  when  the  solution  is at infinity in some of the components of
!  the   unknowns,  x.  See the description of The Option Array for
!  ways to avoid this value of output value of IGO.
!
!  =4   The function f(x) has reached a local minimum.  This is
!  the value of IGO that is expected when a nonzero value of f(x)
!  is anticipated.  See the description of The Option Array for the
!  conditions that have been satisfied.
!
!      =5   The  model  problem  solver  has  noted a value for the
!  linear or quadratic model problem residual vector length that is
!  >=  the  current  value  of  thefunction, i.e. the Euclidean
!  length   of  f(x).   This  situation  probably  means  that  the
!  evaluation  of  f(x)  has  more  uncertainty  or  noise  than is
!  possible  to  account for in the tolerances used to note a local
!  minimum.  The value for x is suspect, but a minimum has probably
!  been found.
!
!      =6  A small change (absolute) was noted for the vector x.  A
!  full  model problem step was taken.  The condition for IGO=4 may
!  also  be  satisfied, so that a minimum has been found.  However,
!  this test is made before the test for IGO=4.
!
!      =7   A  small change (relative to the length of x) was noted
!  for  the  vector  x.   A full model problem step was taken.  The
!  condition for IGO=4 may also be satisfied, so that a minimum has
!  been  found.   However,  this  test  is made before the test for
!  IGO=4.
!
!      =8   More  than  ITMAX  iterations  were taken to obtain the
!  solution.   The  value obtained for x is suspect, although it is
!  the   best   set  of  x  values  that  occurred  in  the  entire
!  computation.   See  the  description  of  The  Option  Array for
!  directions  on  how  to  increase  this  value.   (Note that the
!  nominal  value  for ITMAX, 75, is sufficient to solve all of the
!  nonlinear test problems described in Ref. (2).)
!
!      =9-18     Errors  in the usage of the subprogram were noted.
!  The  exact condition will be noted using an error processor that
!  prints  an  informative  message  unless  this printing has been
!  suppressed.   A  minimum  value  has  not been found for x.  The
!  relation  between  IGO  and  the  error number are IGO=NERR + 8.
!  Here  NERR is the identifying number.  See below, Error Messages
!  for DQED().
!  The Option Array
!  --- ------ -----
!  Glossary of Items Modified by Options.  Those items with Nominal
!  Values listed can be changed.
!
!       Names    Nominal         Definition
!                Values
!       -----    -------         ----------
!       FC                       Current value of length of f(x).
!       FB                       Best value of length of f(x).
!       FL                       Value of length of f(x) at the
!                                previous step.
!       PV                       Predicted value of length of f(x),
!                                after the step is taken, using the
!                                approximating model.
!  The  quantity  'eps',  used  below,  is  the  machine  precision
!  parameter.   Its  value  is obtained by a call to the Bell Labs.
!  Port subprogram D1MACH(4).  It is machine dependent.
!                MIN(1.0D-05,
!       TOLF     sqrt(eps))      Tolerance for stopping when
!                                FC <= TOLF.
!                MIN(1.0D-05,
!       TOLD     sqrt(eps))      Tolerance for stopping when
!                                change to x values has length
!                                <= TOLD.
!                MIN(1.0D-05,
!       TOLX     sqrt(eps))      Tolerance for stopping when
!                                change to x values has length
!                                <= TOLX*length of x values.
!       TOLSNR    1.0D-05          Tolerance used in stopping
!                                condition IGO=4.  Explained below.
!       TOLP      1.0D-05          Tolerance used in stopping
!                                condition IGO=4.  Explained below.
!
!  (The  conditions  (abs(FB-PV)<=TOLSNR*FB  and  abs(FC-PV) .le.
!  TOLP*FB)   and  (ABS(FC-FL)<=TOLSNR*FB) together with  taking
!  a  full  model  step, must be satisfied before the condition  IGO=4
!  is returned.  Decreasing any of the values for  TOLF,  TOLD,  TOLX,
!  TOLSNR,  or  TOLP  will likely increase the number of iterations
!  required for convergence.)
!
!       COND       30.           Largest condition number to allow
!                                when solving for the quadratic
!                                model coefficients.  Increasing
!                                this value may result in more
!                                terms being used in the quadratic
!                                model.
!       TOLUSE   sqrt(eps)       A tolerance that is used to avoid
!                                values of x in the quadratic
!                                model's interpolation of previous
!                                points.  Decreasing this value may
!                                result in more terms being used in
!                                the quadratic model.
!        ITMAX     75            The number of iterations to take
!                                with the algorithm before giving
!                                up and noting it with the value
!                                IGO=8.
!        IPRINT     0            Control the level of printed
!                                output in the solver.  A value
!                                of IPRINT  >  0 will result in
!                                output of information about each
!                                iteration.
!        LEVEL      1            Error processor error level.  See
!                                the SLATEC library documentation
!                                for XERROR() for an explanation.
!        NTERMS     5            One more than the maximum number
!                                of terms used in the quadratic
!                                model.
!
!  IOPT(*) (Input)
!  -------
!  In  order  to  use the option array technique to change selected
!  data within a subprogram, it is necessary to understand how this
!  array  is  processed  within the software.  Let LP designate the
!  processing pointer that moves to positions of the IOPT(*) array.
!  Initially  LP=1,  and  as  each option is noted and changed, the
!  value  of  LP is updated.  The values of IOPT(LP) determine what
!  options get changed.  The amount that LP changes is known by the
!  software  to  be  equal to the value two except for two options.
!  These exceptional cases are the last option (=99) and the 'leap'
!  option  (=13)  which  advances LP by the value in IOPT(LP+1).  A
!  negative  value for IOPT(LP) means that this option is not to be
!  changed.   This aids the programmer in using options;  often the
!  code  for  using  an  option can be in the calling program but a
!  negative value of the option number avoids rewriting code.
!
!  Option Usage Example
!  ------ ----- -------
!  In  the  Fortran code fragment that follows, an example is given
!  where  we  change  the  value  of  TOLF and decrease the maximum
!  number  of  iterations  allowed  from  75  to 30.
!  In this example the dimensions of IOPT(*) and ROPT(*) must
!  satisfy:
!
!        double precision ROPT(01)
!        integer IOPT(005)
!        .
!        .
!        .
!  C     SET THE OPTION TO CHANGE THE VALUE OF TOLF.
!
!        IOPT(01)=4
!
!  C     THE NEXT ENTRY POINTS TO THE PLACE IN ROPT(*) WHERE
!  C     THE NEW VALUE OF TOLF IS LOCATED.
!
!        IOPT(02)=1
!  C     THIS IS THE NEW VALUE OF TOLF.  THE SPECIFIC VALUE
!  C     1.0D-09 IS USED HERE ONLY FOR ILLUSTRATION.
!
!        ROPT(01)=1.0D-09
!
!  C     CHANGE THE NUMBER OF ITERATIONS.
!
!        IOPT(03)=2
!
!  C     THIS NEXT ENTRY IS THE NEW VALUE FOR THE MAXIMUM NUMBER OF
!  C     ITERATIONS.
!
!        IOPT(04)=30
!
!  C     THIS NEXT OPTION IS A SIGNAL THAT THERE ARE NO MORE
!  C     OPTIONS.
!
!        IOPT(05)=99
!        .
!        .
!        .
!        CALL DQED()
!        .
!        .
!        .
!  Option Values   Explanation
!  ------ ------   -----------
!     =99          There are no more options to change.
!                  Normally this is the first and only
!                  option that a user needs to specify,
!                  and it can be simply IOPT(01)=99.  The
!                  total dimension of IOPT(*) must be at
!                  least 17, however.  This can lead to a
!                  hard-to-find program bug if the dimension
!                  is too small.
!
!     = 1          Change the amount of printed output.
!                  The next value of IOPT(*) is the print
!                  level desired, IPRINT.  Any value of
!                  IPRINT  >  0 gives all the available
!                  output.
!
!     = 2          Change the value of ITMAX.  The next value
!                  of IOPT(*) is the value of ITMAX desired.
!
!     = 3          Pass prior determined bounds for the box
!                  containing the initial point.  This box is the
!                  trust region for the first move from the initial
!                  point.  The next entry in IOPT(*) points to
!                  the place in ROPT(*) where the NVARS values for
!                  the edges of the box are found.
!
!     = 4          Change the value of TOLF.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of TOLF is found.
!
!     = 5          Change the value of TOLX.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of TOLX is found.
!
!     = 6          Change the value of TOLD.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of TOLD is found.
!
!     = 7          Change the value of TOLSRN.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of TOLSNR is found.
!
!     = 8          Change the value of TOLP.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of TOLP is found.
!
!     = 9          Change the value of TOLUSE.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of TOLUSE is found.
!
!     =10          Change the value of COND.  The next entry of
!                  IOPT(*) points to the place in ROPT(*) where the
!                  new value of COND is found.
!
!     =11          Change the value of LEVEL.  The next entry of
!                  IOPT(*) is the new value of LEVEL.
!
!     =12          Pass an option array to the subprogram DQEDGN()
!                  used as the inner loop solver for the
!                  model problem.  The next entry of IOPT(*) is the
!                  starting location for the option array for
!                  DQEDGN() within the array IOPT(*).  Thus the
!                  option array for DQEDGN() must be a part of
!                  the array IOPT(*).
!
!     =13          Move (or leap) the processing pointer LP for the
!                  option array by the next value in IOPT(*).
!
!     =14          Change a logical flag that suppresses the
!                  use of the quadratic model in the inner
!                  loop.  Use the next value in IOPT(*) for
!                  this flag.  If this value = 1, then never
!                  use the quadratic model.  (Just use the
!                  linear model).  Otherwise, use the quadratic
!                  model when appropriate.  This option decreases
!                  the amount of scratch storage as well as the
!                  computing overhead required by the code package.
!                  A user may want to determine if the application
!                  really requires the use of the quadratic model.
!                  If it does not, then use this option to save
!                  both storage and computing time.
!
!     =15          Change, NTERMS,  the maximum number of array
!                  columns that can be used for saving quadratic
!                  model data.  (The value of NTERMS is one more
!                  than the maximum number of terms used.)  Each
!                  unit increase for NTERMS increases the required
!                  dimension of the array WORK(*) by 2*MEQUA+NVARS.
!                  Use the value in IOPT(LP+1) for the new value
!                  of NTERMS.  Decreasing this value to 2 (its
!                  minimum) decreases the amount of storage
!                  required by the code package.
!
!     =16          Change a logical flag so that 'reverse
!                  communication' is used instead of 'forward
!                  communication.'  Example EX01, listed below,
!                  uses 'forward communication.'  Example EX02,
!                  also listed below, uses 'reverse communication.'
!                  Use the next value in IOPT(*) for
!                  this flag.  If this value = 1, then
!                  use 'reverse communication.'  Otherwise,
!                  use 'forward communication.'  WARNING:  This
!                  usage may not work unless the operating system
!                  saves variables between subroutine calls to DQED.
!
!     =17          Do not allow the flag IGO to return with the
!                  value IGO=3.  This means that convergence will
!                  not be claimed unless a full model step is taken.
!                  Normal output values will then be IGO = 2,4,6 or 7.
!                  Use the next value in IOPT(*) for this flag.  If
!                  this value = 1, then force a full model step.
!                  Otherwise,  do not force a full model step if small
!                  steps are noted.
!
!  IWORK(*), WORK(*) (Input and Output)
!  ----------------
!  These  are  scratch arrays that the software uses for storage of
!  intermediate  results.   It  is  important  not  to  modify  the
!  contents of this storage during the computation.
!
!  The  array  locations  IWORK(1)  and  IWORK(2)  must contain the
!  actual  lengths  of  the  arrays WORK(*) and IWORK(*) before the
!  call to the subprogram.  These array entries are replaced by the
!  actual amount of storage required for each array.  If the amount
!  of  storage  for either array is too small, an informative error
!  message will be printed, and the value IGO=13 or 14 returned.
!
!  The  user may find it useful to let the subprogram DQED() return
!  the  amounts  of storage required for these arrays.  For example
!  set  IWORK(1)=1,  IWORK(2)=1.   The  subprogram will return with
!  IGO=13,     IWORK(1)=required    length    of    WORK(*),    and
!  IWORK(2)=required   length   of  IWORK(*).   (Appropriate  error
!  messages will normally be printed.)
!
!  3. Remarks on the Usage Examples
!     ------- -- --- ----- --------
!  The  following  complete  program units, EX01 and EX02, show how
!  one  can  use  the  nonlinear  solver  for  fitting  exponential
!  functions  to  given data.  These examples are calculations that
!  match  two  terms  of  an  exponential series to five given data
!  points.   There are some subtle points about exponential fitting
!  that   are  important  to  note.     First,  the  signs  of  the
!  exponential arguments are restricted to be nonpositive.
!  The size of the arguments should not be much larger than the start
!  of the time data (reciprocated).  This is the reason the lower
!  bounds are set a bit less than the reciprocal of the time value.
!  In many applications that require exponential modeling this is a
!  natural assumption.  The nonlinear solver allows these bounds
!  on  the arguments explicitly.  In addition, the coefficients are
!  constrained  to  be  nonnegative.   These  bounds  are harder to
!  justify.  The idea is to avoid the situation where a coefficient
!  is  very  large  and negative, and the corresponding exponential
!  argument is also large and negative.  The resulting contribution
!  to  the  series may be very small, but its presence is spurious.
!  Finally,  the  single  general  linear constraint that keeps the
!  arguments  separated  (by  0.05 in this example) is used for two
!  purposes.   First,  it naturally orders these values so that the
!  first  one  is  algebraically  largest.  Second, this constraint
!  moves the parameters from the local minimum corresponding to the
!  initial  values  used  in  the  examples.   This constraint also
!  retains  the  validity of the model function h(t) = w*exp(x*t) +
!  y*exp(z*t).  Namely, if the arguments are allowed to coalesce to
!  the  same value, then the model itself must change.  The form of
!  the model can become h(t)=(a+b*t)*exp(c*t) or h(t) = d*exp(e*t).
!  Either one could occur, and the choice is problem dependent.
!
!
!  Example 1  Using Forward Communication
!  ---------  ----- ------- -------------
!
!      PROGRAM EX01
!
!C     Illustrate the use of the Hanson-Krogh nonlinear least
!C     squares solver for fitting two exponentials to data.
!C
!C     The problem is to find the four variables x(1),...,x(4)
!C     that are in the model function
!C
!C          h(t) = x(1)*exp(x(2)*t) + x(3)*exp(x(4)*t)
!C     There are values of h(t) given at five values of t,
!C     t=0.05, 0.1, 0.4, 0.5, and 1.0.
!C     We also have problem constraints that x(2), x(4) <= 0, x(1),
!C     x(3) >= 0, and a minimal separation of 0.05 between x(2) and
!C     x(4).  Nothing more about the values of the parameters is known
!C     except that x(2),x(4) are approximately >= 1/min t.
!C     Thus we have no further knowledge of their values.
!C     For that reason all of the initial values are set to zero.
!C
!C     Dimension for the nonlinear solver.
!      double precision FJ(6,5),BL(5),BU(5),X(4),ROPT(001),WA(640)
!C  EDIT on 950228-1300:
!      double precision RNORM
!      integer IND(5),IOPT(24),IWA(084)
!
!      EXTERNAL DQEDEX
!
!      DATA LDFJ,LWA,LIWA/6,640,084/
!
!      MCON = 1
!      MEQUA = 5
!      NVARS = 4
!C     Define the constraints for variables.
!      BL(1) = 0.
!      BL(2) = -25.
!      BU(2) = 0.
!      BL(3) = 0.
!      BL(4) = -25.
!      BU(4) = 0.
!C     Define the constraining value (separation) for the arguments.
!      BL(5) = 0.05
!C     Define all of the constraint indicators.
!      IND(1) = 1
!      IND(2) = 3
!      IND(3) = 1
!      IND(4) = 3
!      IND(5) = 1
!C     Define the initial values of the variables.
!C     We don't know anything more, so all variables are set zero.
!      x(1:nvars) = 0.0D+00
!   10 CONTINUE
!C     Tell how much storage we gave the solver.
!      IWA(1) = LWA
!      IWA(2) = LIWA
!C     No additional options are in use.
!      IOPT(01) = 99
!      CALL DQED(DQEDEX,MEQUA,NVARS,MCON,IND,BL,BU,X,FJ,LDFJ,RNORM,IGO,
!     .          IOPT,ROPT,IWA,WA)
!      WRITE (*,9001) (X(J),J=1,NVARS)
!      WRITE (*,9011) RNORM
!      WRITE (*,9021) IGO
!
!      STOP
!
! 9001 FORMAT (' MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))',/,
!     .  ' X(1),X(2),X(3),X(4) = ',/,4F12.6)
! 9011 FORMAT (' RESIDUAL AFTER THE FIT = ',1PD12.4)
! 9021 FORMAT (' OUTPUT FLAG FROM SOLVER =',17X,I6)
!      END
!      SUBROUTINE DQEDEX(X,FJ,LDFJ,IGO,IOPT,ROPT)
!C     This is the subprogram for evaluating the functions
!C     and derivatives for the nonlinear solver, DQED.
!C
!C     The user problem has MCON constraint functions,
!C     MEQUA least squares equations, and involves NVARS
!C     unknown variables.
!C
!C     When this subprogram is entered, the general (near)
!C     linear constraint partial derivatives, the derivatives
!C     for the least squares equations, and the associated
!C     function values are placed into the array FJ(*,*).
!C     All partials and functions are evaluated at the point
!C     in X(*).  Then the subprogram returns to the calling
!C     program unit. Typically one could do the following
!C     steps:
!C
!C     step 1. Place the partials of the i-th constraint
!C           function with respect to variable j in the
!C             array FJ(i,j), i=1,...,MCON, j=1,...,NVARS.
!C     step 2. Place the values of the i-th constraint
!C             equation into FJ(i,NVARS+1).
!C     step 3. Place the partials of the i-th least squares
!C             equation with respect to variable j in the
!C             array FJ(MCON+i,j), i=1,...,MEQUA,
!C             j=1,...,NVARS.
!C     step 4. Place the value of the i-th least squares
!C             equation into FJ(MCON+i,NVARS+1).
!C     step 5. Return to the calling program unit.
!      double precision FJ(LDFJ,*),X(*),ROPT(*)
!      double precision T(5),F(5)
!      integer IOPT(*)
!
!      DATA T/0.05,0.10,0.40,0.50,1.00/
!      DATA F/2.206D+00,1.994D+00,1.350D+00,1.216D+00,.7358D0/
!
!      DATA MCON,MEQUA,NVARS/1,5,4/
!
!C     Define the derivatives of the constraint with respect to the x(j).
!      FJ(1,1) = 0.0D+00
!      FJ(1,2) = 1.0D+00
!      FJ(1,3) = 0.0D+00
!      FJ(1,4) = -1.0D+00
!C     Define the value of this constraint.
!      FJ(1,5) = X(2) - X(4)
!C     Define the derivatives and residuals for the data model.
!      DO I = 1,MEQUA
!         E1 = EXP(X(2)*T(I))
!         E2 = EXP(X(4)*T(I))
!         FJ(MCON+I,1) = E1
!         FJ(MCON+I,2) = X(1)*T(I)*E1
!         FJ(MCON+I,3) = E2
!         FJ(MCON+I,4) = X(3)*T(I)*E2
!         FJ(MCON+I,5) = X(1)*E1 + X(3)*E2 - F(I)
!      end do
!      RETURN
!      END
!  Output from Example 1 Program
!  ------ ---- --------- -------
!
!   MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))
!  X(1),X(2),X(3),X(4) =
!      1.999475    -.999801     .500057   -9.953988
!   RESIDUAL AFTER THE FIT =   4.2408D-04
!   OUTPUT FLAG FROM SOLVER =                      4
!
!
!  Example 2  Using Reverse Communication
!  ---------  ----- ------- -------------
!      PROGRAM EX02
!
!C     Illustrate the use of the Hanson-Krogh nonlinear least
!C     squares solver for fitting two exponentials to data.
!C
!C     The problem is to find the four variables x(1),...,x(4)
!C     that are in the model function
!C
!C          h(t) = x(1)*exp(x(2)*t) + x(3)*exp(x(4)*t)
!C     There are values of h(t) given at five values of t,
!C     t=0.05, 0.1, 0.4, 0.5, and 1.0.
!C     We also have problem constraints that x(2), x(4) <= 0, x(1),
!C     x(3) >= 0, and a minimal separation of 0.05 between x(2) and
!C     x(4).  Nothing more about the values of the parameters is known
!C     except that x(2),x(4) are approximately >= 1/min t.
!C     Thus we have no further knowledge of their values.
!C     For that reason all of the initial values are set to zero.
!C
!C     Dimension for the nonlinear solver.
!      double precision FJ(6,5),BL(5),BU(5),X(4),ROPT(001),WA(640)
!C  EDIT on 950228-1300:
!      double precision RNORM
!      integer IND(5),IOPT(24),IWA(084)
!      double precision T(5),F(5)
!
!      EXTERNAL DQEDEV
!
!      DATA LDFJ,LWA,LIWA/6,640,084/
!
!      DATA T/0.05,0.10,0.40,0.50,1.00/
!      DATA F/2.206D+00,1.994D+00,1.350D+00,1.216D+00,.7358D0/
!
!      MCON = 1
!      MEQUA = 5
!      NVARS = 4
!C     Define the constraints for variables.
!      BL(1) = 0.
!      BL(2) = -25.
!      BU(2) = 0.
!      BL(3) = 0.
!      BL(4) = -25.
!      BU(4) = 0.
!C     Define the constraining value (separation) for the arguments.
!      BL(5) = 0.05
!C     Define all of the constraint indicators.
!      IND(1) = 1
!      IND(2) = 3
!      IND(3) = 1
!      IND(4) = 3
!      IND(5) = 1
!C     Define the initial values of the variables.
!C     We don't know anything at all, so all variables are set zero.
!      x(1:nvars) = 0.0D+00
!C     Tell how much storage we gave the solver.
!      IWA(1) = LWA
!      IWA(2) = LIWA
!      NITERS = 0
!C     TELL HOW MUCH STORAGE WE GAVE THE SOLVER.
!      IWA(1) = LWA
!      IWA(2) = LIWA
!C     USE REVERSE COMMUMICATION TO EVALUATE THE DERIVATIVES.
!      IOPT(01)=16
!      IOPT(02)=1
!C     NO MORE OPTIONS.
!      IOPT(03) = 99
!   20 CONTINUE
!      CALL DQED(DQEDEV,MEQUA,NVARS,MCON,IND,BL,BU,X,FJ,LDFJ,RNORM,
!     .IGO,IOPT, ROPT,IWA,WA)
!      IF (IGO.GT.1) GO TO 40
!C     COUNT FUNCTION EVALUATIONS.
!      NITERS = NITERS + 1
!C     DEFINE THE DERIVATIVES OF THE CONSTRAINT WITH RESPECT TO THE X(J).
!      FJ(1,1) = 0.0D+00
!      FJ(1,2) = 1.0D+00
!      FJ(1,3) = 0.0D+00
!      FJ(1,4) = -1.0D+00
!C     DEFINE THE VALUE OF THIS CONSTRAINT.
!      FJ(1,5) = X(2) - X(4)
!C     DEFINE THE DERIVATIVES AND RESIDUALS FOR THE DATA MODEL.
!      DO I = 1,MEQUA
!          E1 = EXP(X(2)*T(I))
!          E2 = EXP(X(4)*T(I))
!          FJ(MCON+I,1) = E1
!          FJ(MCON+I,2) = X(1)*T(I)*E1
!          FJ(MCON+I,3) = E2
!          FJ(MCON+I,4) = X(3)*T(I)*E2
!          FJ(MCON+I,5) = X(1)*E1 + X(3)*E2 - F(I)
!      end do
!      GO TO 20
!
!   40 CONTINUE
!      WRITE (*,9001) (X(J),J=1,NVARS)
!      WRITE (*,9011) RNORM
!      WRITE (*,9021) NITERS, IGO
!
! 9001 FORMAT (' MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))',/,
!     . ' X(1),X(2),X(3),X(4) = ',/,4F12.6)
! 9011 FORMAT (' RESIDUAL AFTER THE FIT = ',1PD12.4)
! 9021 FORMAT (' NUMBER OF EVALUATIONS OF PARAMETER MODEL =',I6,/,
!     .          ' OUTPUT FLAG FROM SOLVER =',17X,I6)
!      STOP
!      END
!  Output from Example 2 Program
!  ------ ---- --------- -------
!
!  MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))
!  X(1),X(2),X(3),X(4) =
!      1.999475    -.999801     .500057   -9.953988
!   RESIDUAL AFTER THE FIT =   4.2408D-04
!   NUMBER OF EVALUATIONS OF PARAMETER MODEL =    14
!   OUTPUT FLAG FROM SOLVER =                      4
!
!  5. References
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide, 
!    SIAM (Society for Industrial and Applied Mathematics),
!    Philadelphia, 1979.
!
!    Richard Hanson,
!    Least Squares with Bounds and Linear Constraints,
!    SIAM Journal of Scientific and Statistical Computing, 
!    Volume 7, Number 3, July 1986, pages 826-834.
!
!    Robert Schnabel and P Frank,
!    Tensor Methods for Nonlinear Equations,
!    SIAM Journal on Numerical Analysis, 
!    Volume 21, Number 5, October 1984, pages 815-843.
! 
!***END PROLOGUE  DQED
!     REVISED 870204-1100
!     REVISED YYMMDD-HHMM
!
  implicit none

  integer ldfjac

  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) bu(*)
  external :: dqedev
  real ( kind = 8 ) fjac(ldfjac,*)
  real ( kind = 8 ) fnorm
  integer i
  integer idum
  integer, save :: iflag = 0
  integer igo
  integer, save :: iiwaav
  integer ind(*)
  integer iopt(*)
  integer iwa(*)
  integer iwaav
  integer j
  integer jp
  integer kp
  integer, save :: level
  integer lp
  integer lpdiff
  integer mb
  integer mbb
  integer mblb
  integer mbub
  integer mcon
  integer mdx
  integer mdxl
  integer mequa
  integer mgr
  integer, save :: milast
  integer, save :: mind
  integer, save :: mindb
  integer, save :: mpj
  integer, save :: mqc
  integer, save :: mwa
  integer, save :: mwj
  integer, save :: mwlast
  integer, save :: mxb
  integer, save :: mxp
  integer, save :: mzp
  integer, save :: nall
  logical, save :: noquad
  integer, save :: npmax
  integer nvars
  real ( kind = 8 ), save :: rdum
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) wa(*)
  real ( kind = 8 ) x(*)
!
! Name      Memory Status  Type     Argument   Uses and comments.
!                                    Status
! ----      -------------  ----     --------   ------------------
! BL         DUMMY-ARG     REAL      ADJ-ARY Lower bounds
! BU         DUMMY-ARG     REAL      ADJ-ARY Upper bounds
! FJAC       DUMMY-ARG     REAL      ADJ-ARY Jacobian array
! FNORM      DUMMY-ARG     REAL              Norm at solution
! I          /S$A$V$E/ SAV integer           Dummy loop variable
! IDUM       /S$A$V$E/ SAV integer           Dummy for error pack
! IFLAG      /S$A$V$E/ SAV integer           Gate for reverse comm
! IGO        DUMMY-ARG     integer           Directs user action
! IGOOK      /S$A$V$E/ SAV integer           Internal gate for errors
! IIWAAV     /S$A$V$E/ SAV integer           Length claimed for IWA
! IND        DUMMY-ARG     integer   ADJ-ARY Bound indicators
! IOPT       DUMMY-ARG     integer   ADJ-ARY Option array
! IWA        DUMMY-ARG     integer   ADJ-ARY Work array
! IWAAV      /S$A$V$E/ SAV integer           Length claime for WA
! J          /S$A$V$E/ SAV integer           Dummy loop variable
! MILAST     /S$A$V$E/ SAV integer           Last integer in IWA
! MIND       /S$A$V$E/ SAV integer           Point to start IND
! MINDB      /S$A$V$E/ SAV integer           Point to start INDB
! MPJ        /S$A$V$E/ SAV integer           Point to start PJ
! MQC        /S$A$V$E/ SAV integer           Point to start QC
! MUT        /S$A$V$E/ SAV integer           Point to start UT
! MWA        /S$A$V$E/ SAV integer           Point to start WA
! MWJ        /S$A$V$E/ SAV integer           Point to start WJ
! MWLAST     /S$A$V$E/ SAV integer           Last value in WA
! MXB        /S$A$V$E/ SAV integer           Point to start XB
! MXP        /S$A$V$E/ SAV integer           Point to start XP
! MZP        /S$A$V$E/ SAV integer           Point to start ZP
! NALL       /S$A$V$E/ SAV integer           Sum of dimensions
! KP         /S$A$V$E/ SAV integer           Dummy option loop pointer
! LDFJAC     DUMMY-ARG     integer           Row dimension of FJAC
! LEVEL      /S$A$V$E/ SAV integer           Error processor status
! LP         /S$A$V$E/ SAV integer           Dummy option loop pointer
! LPDIFF     /S$A$V$E/ SAV integer           Dummy option loop diff
! MB         /S$A$V$E/ SAV integer           Point to start B
! MBB        /S$A$V$E/ SAV integer           Point to start BB
! MBLB       /S$A$V$E/ SAV integer           Point to start BLB
! MBUB       /S$A$V$E/ SAV integer           Point to start BUB
! MCON       DUMMY-ARG     integer           Number of constraints
! MDX        /S$A$V$E/ SAV integer           Point to start DX
! MDXL       /S$A$V$E/ SAV integer           Point to start DXL
! MEQUA      DUMMY-ARG     integer           Numer of least squares eqs
! MGR        /S$A$V$E/ SAV integer           Point to start GR
! NERR       /S$A$V$E/ SAV integer           Error processor number
! NMESS      /S$A$V$E/ SAV integer           Length of error message
! NOQUAD     /S$A$V$E/ SAV LOGICAL           Flag, suppress quad model
! NPMAX      /S$A$V$E/ SAV integer           Max number of quad terms
! NVARS      DUMMY-ARG     integer           Number of unknowns
! RDUM       /S$A$V$E/ SAV REAL              Dummy variable, error proc
! ROPT       DUMMY-ARG     REAL      ADJ-ARY Option data
! WA         DUMMY-ARG     REAL      ADJ-ARY Working array
! X          DUMMY-ARG     REAL      ADJ-ARY Values of the variables
! XMESS      /S$A$V$E/ SAV CHAR*128          Hold error message
!
  idum = 0
  rdum = 0.0D+00

  if ( iflag == 0 ) then

      noquad = .false.
      level = 1

      if ( mequa <= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'DQED - Fatal error!'
        write ( *, * ) '  MEQUA must be greater than 0.'
        write ( *, * ) '  Input value is ', mequa
        stop
      end if

      if ( nvars <= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'DQED - Fatal error!'
        write ( *, * ) '  NVARS must be greater than 0.'
        write ( *, * ) '  Input value is ', nvars
        stop
      end if

      if ( mcon < 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'DQED - Fatal error!'
        write ( *, * ) '  MCON must be greater than 0.'
        write ( *, * ) '  Input value is ', mcon
        stop
      end if

      do j = 1, nvars + mcon
         i = ind(j)
         if ( i < 1 .or. 4 < i ) then
           write ( *, * ) ' '
           write ( *, * ) 'DQED - Fatal error!'
           write ( *, * ) '  IND(J) must be greater than between 1 and 4.'
           write ( *, * ) '  Input value of IND(',j,') is ', ind(j)
           stop
         end if
      end do

      npmax = 5
!
!  LOOK THROUGH THE OPTION ARRAY FOR A CHANGE TO NPMAX,
!  THE AMOUNT OF ARRAY STORAGE USED FOR THE QUADRATIC PARAMETERS.
!
      lp = 1
      lpdiff = 0

   20 continue

      lp = lp + lpdiff
      lpdiff = 2
      kp = iopt(lp)
      jp = abs(kp)

      if ( jp == 99) then
        if ( kp > 0) go to 30
      end if
!
!  THIS IS THE ONLY OPTION WHERE THE PROCESSING POINTER
!  MUST BE CHANGED FROM THE VALUE 2.
!
      if ( jp == 13) then
        lpdiff = iopt(lp+1)
      end if
!
!  FOUND A CHANGE TO THE ARRAY SIZE FOR THE QUADRATIC MODEL.
!
      if ( jp == 15) then
          if ( kp > 0) npmax = iopt(lp+1)
      end if
!
!  SEE IF THE QUADRATIC MODEL IS SUPPRESSED.
!  THIS REQUIRES LESS STORAGE IN THE USER PROGRAM.
!
      if ( jp == 14) then
          if ( kp > 0) noquad = iopt(lp+1)  ==  1
      end if

      if ( jp < 1 .or. 17 < jp ) then
        write ( *, * ) ' '
        write ( *, * ) 'DQED - Fatal error!'
        write ( *, * ) '  Invalid option JP = ', jp
        stop
      end if

      go to 20

   30 continue

  end if

  if ( noquad) then
    nall = mcon + nvars + 2
  else
    nall = mcon + 2*nvars + npmax + 1
  end if
!
!  COMPUTE POINTERS INTO WORK SPACE FOR VARIOUS ARRAYS.
!
  mdx = 1
  mxb = mdx + nall + 2
  if ( mcon > 0)mxb=mxb+nall+2
  mb = mxb + nvars
  mbb = mb + nvars
  mblb = mbb + nvars
  if ( mcon > 0)mblb=mblb+nall
  mbub = mblb + nall
  if ( mcon > 0)mbub=mbub+nall
  mzp = mbub + nall
  mxp = mzp + mequa*npmax
  mqc = mxp + nvars*npmax
  mwj = mqc + max(mequa,nvars)*npmax
  mpj = mwj + nall* (nall+1)
  mgr = mpj + nvars + 1
  mdxl = mgr + nvars
  mwa = mdxl + nvars + nvars
  mwlast = mwa + 9* (mcon+1) + 13* (2*nvars+npmax+1)

  mindb = 3
  mind = mindb + nall + nvars
  milast = mind + 3* (mcon+1) + 4* (2*nvars+npmax+1)
!
!  CHECK LENGTHS OF ARRAYS ONCE PER PROBLEM.
!
  if ( iflag == 0) then

      iwaav = iwa(1)
      iiwaav = iwa(2)
!
!  RETURN THE ACTUAL AMOUNTS OF STORAGE REQD. FOR WA(*), IWA(*).
!
      iwa(1) = mwlast
      iwa(2) = milast

      if ( iwaav < mwlast) then
        write ( *, * ) ' '
        write ( *, * ) 'DQED - Fatal error!'
        write ( *, * ) '  Insufficient storage in WA.'
        write ( *, * ) '  Amount needed: ', mwlast
        write ( *, * ) '  Amount given:  ', iwaav
        stop
      end if

      if ( iiwaav < milast) then
        write ( *, * ) ' '
        write ( *, * ) 'DQED - Fatal error!'
        write ( *, * ) '  Insufficient storage in IWA.'
        write ( *, * ) '  Amount needed: ', milast
        write ( *, * ) '  Amount given:  ', iiwaav
        write ( *, * ) ' '
      end if

      iflag = 1

  end if

  call dqedmn ( dqedev, mequa, nvars, mcon, ind, bl, bu, x, fjac, &
    ldfjac, fnorm, igo, iopt, ropt, iwa(mind), wa(mwa), wa(mdx), wa(mxb), &
    wa(mb), wa(mbb), wa(mblb), wa(mbub), iwa(mindb), npmax, &
    wa(mzp), wa(mxp), wa(mqc), max(mequa,nvars), wa(mpj), &
    wa(mwj), nall, wa(mgr), wa(mdxl) )

  if ( 1 < igo ) then
    iflag = 0
  end if

  return
end
subroutine dqedev ( x, fj, ldfj, igo, iopt, ropt )

!***********************************************************************
!
!! DQEDEV evaluates functions being treated by DQED.
!
!  Discussion:
!
!    The user has NVARS variables X(I), and is trying to minimize
!    the square root of the sum of the squares of MEQUA functions
!    F(I)(X), subject to MCON constraints which have the form
!
!      BL(I) <= G(I)(X) <= BU(I)
!
!    where either the left or right bounding inequality may be dropped.
!
!  Parameters:
!
!    Input, REAL X(*), an array of length NVARS, containing
!    the values of the independent variables at which the
!    functions and partial derivatives should be evaluated.
!
!    Output, REAL FJ(LDFJ,NVARS+1).
!
!    If IGO is nonzero, then partial derivatives must
!    be placed in the first NVARS columns of FJ, as follows:
!
!      Rows I = 1 to MCON, and columns J = 1 to NVARS
!      should contain the partials of the I-th constraint
!      function G(I) with respect to the J-th variable.
!
!      Rows I=MCON+1 to MCON+MEQUA, and columns J = 1 to NVARS
!      should contain the partials of the (I-MCON)-th nonlinear
!      function F(I-MCON) with respect to the J-th variable.
!
!    Regardless of the value of IGO, column NVARS+1 must be
!    assigned the values of the constraints and functions, as follows:
!
!      Rows I = 1 to MCON, column J = NVARS+1, should contain
!      the value of G(I).
!
!      Rows I=MCON+1 to MCON+MEQUA, column J = NVARS+1, should
!      contain the value of F(I-MCON).
!
!    Input, integer LDFJ, the leading dimension of FJ, which
!    must be at least MCON+MEQUA.
!
!    Input/output, integer IGO.
!
!    On input, IGO tells the user whether the partial derivatives are needed.
!
!      0, the partial derivatives are not needed.
!      nonzero, the partial derivatives are needed.
!
!    On output, the user may reset the input value of IGO if one
!    of two situations is encountered:
!
!      99, the functions, constraints, or partial derivatives
!          could not be evaluated at the given input point X.  Request
!          that DQED reject that point, and try a different one.
!
!      Any other value, abort the run.
!
!    Input, integer IOPT(*), the integer option array.
!
!    Input, REAL ROPT(*), the double precision option array.
!
  implicit none

  integer ldfj

  real ( kind = 8 ) fj(ldfj,*)
  integer igo
  integer iopt(*)
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) x(*)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQEDEV - Fatal error!'
  write ( *, '(a)' ) '  DQEDEV must be written by the user.'

  stop
end
subroutine dqedgn ( mequa, nvars, mcon, ind, bl, bu, x, fjac, ldfjac, &
  fnorm, igo, iopt, ropt, iwa, wa )

!***********************************************************************
!
!! DQEDGN is a simplified version of the QED algorithm for the model problem.
!
!***BEGIN PROLOGUE  DQEDGN
!***REFER TO  DQED
!***ROUTINES CALLED  DQEDIP
!***END PROLOGUE  DQEDGN
!  DQEDGN:
! GLOSSARY OF VARIABLES. NOTATION:
! DUMMY-ARG A dummy argument, that is an argument to this prog. unit.
! /S$A$V$E/ SAV Denotes that this variable is local to the routine
!               and is saved between calls to it.
! integer, REAL, real, LOGICAL, CHARACTER
!               The types of the variables.
! ADJ-ARR An adjustable array, that is an argument to this prog. unit.
! Name      Memory Status  Type     Argument   Uses and comments.
!                                    Status
! ----      -------------  ----     --------   ------------------
! BL         DUMMY-ARG     REAL      ADJ-ARY Model lower bounds
! BU         DUMMY-ARG     REAL      ADJ-ARY Model upper bounds
! FJAC       DUMMY-ARG     REAL      ADJ-ARY Model Jacobian array
! FNORM      DUMMY-ARG     REAL              Model residual norm
! IGO        DUMMY-ARG     integer           direct model action
! IND        DUMMY-ARG     integer   ADJ-ARY Model bound indicators
! IOPT       DUMMY-ARG     integer   ADJ-ARY Option array
! IWA        DUMMY-ARG     integer   ADJ-ARY Working array
! LDFJAC     DUMMY-ARG     integer           Row dim of FJAC(*,*)
! MB         /S$A$V$E/ SAV integer           Pointer to B(*)
! MBB        /S$A$V$E/ SAV integer           Pointer to BB(*)
! MBLB       /S$A$V$E/ SAV integer           Pointer to BLB(*)
! MBUB       /S$A$V$E/ SAV integer           Pointer to BLB(*)
! MCON       DUMMY-ARG     integer           Number, model constraints
! MDX        /S$A$V$E/ SAV integer           Pointer to DX(*)
! MEQUA      DUMMY-ARG     integer           Number, model equations
! MINDB      /S$A$V$E/ SAV integer           Pointer to INDB(*)
! MIWA       /S$A$V$E/ SAV integer           Pointer to IWA(*)
! MWA        /S$A$V$E/ SAV integer           Pointer to WA(*)
! MXB        /S$A$V$E/ SAV integer           Pointer to XB(*)
! NALL       /S$A$V$E/ SAV integer           NVARS+MEQUA
! NVARS      DUMMY-ARG     integer           Number, user variables
! ROPT       DUMMY-ARG     REAL      ADJ-ARY Option array data
! WA         DUMMY-ARG     REAL      ADJ-ARY Working array
!
  implicit none

  integer ldfjac

  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) bu(*)
  real ( kind = 8 ) fjac(ldfjac,*)
  real ( kind = 8 ) fnorm
  integer igo
  integer ind(*)
  integer iopt(*)
  integer iwa(*)
  integer mb
  integer mbb
  integer mblb
  integer mbub
  integer mcon
  integer mdx
  integer mequa
  integer mindb
  integer miwa
  integer mwa
  integer mxb
  integer nall
  integer nvars
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) wa(*)
  real ( kind = 8 ) x(*)

  save
!
!  ALLOCATE BLOCKS OF WORKING STORAGE TO LOGICAL ARRAYS.
!
  nall = mcon + nvars
  mdx = 1
  mxb = mdx + 2*nall + 2
  mb = mxb + nvars
  mbb = mb + nvars
  mblb = mbb + nvars
  mbub = mblb + nall
  mwa = mbub + nall

  mindb = 1
  miwa = mindb + nall + nvars

  call dqedip(mequa,nvars,mcon,ind,bl,bu,x,fjac,ldfjac,fnorm,igo, &
    iopt,ropt,iwa(miwa),wa(mwa),wa(mdx),wa(mxb),wa(mb), &
    wa(mbb),wa(mblb),wa(mbub),iwa(mindb))
!
!  THESE DEFINE THE AMOUNT OF STORAGE FOR THE double precision AND
!  integer WORK ARRAYS, WA(*) AND IWA(*).
!
  mwa = mwa + 6*nvars + 5*mcon
  miwa = miwa + 2*nall
!
!  TOTAL WORKING STORAGE IN WA(*)=
!    9*NALL + 4*NVARS = 9*MCON + 13*NVARS.
!  TOTAL WORKING STORAGE IN IWA(*)=
!    3*NALL + NV      = 3*MCON +4*NVARS.
!
  return
end
subroutine dqedip ( mequa, nvars, mcon, ind, bl, bu, x, fjac, ldfjac, &
  fb, igo, iopt, ropt, iwa, wa, dx, xb, b, bb, blb, bub, indb )

!***********************************************************************
!
!! DQEDIP carries out the work of DQEDGN.
!
!
  implicit none

  integer ldfjac

  real ( kind = 8 ) alb
  real ( kind = 8 ) alfac
  real ( kind = 8 ) alpha
  real ( kind = 8 ) aub
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) bb(*)
  real ( kind = 8 ) bboost
  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) blb(*)
  real ( kind = 8 ) bold
  real ( kind = 8 ) bu(*)
  real ( kind = 8 ) bub(*)
  real ( kind = 8 ) c1516
  real ( kind = 8 ) chg
  real ( kind = 8 ) chgfac
  real ( kind = 8 ) colnrm
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dxnrm
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) fjac(ldfjac,*)
  real ( kind = 8 ) fl
  logical fulnwt
  real ( kind = 8 ) gval
  integer icase
  integer idamax
  integer :: iflag = 0
  integer igo
  integer igoiov
  integer igopoa
  integer igotfc
  integer igotnc
  integer ind(*)
  integer indb(*)
  integer iopt(*)
  integer ipls
  integer iprint
  integer iters
  integer itmax
  integer iwa(*)
  integer j
  integer jp
  integer k
  integer kl
  integer kp
  integer level
  integer lp
  integer lpdiff
  integer mcon
  integer mequa
  integer mode
  integer nall
  integer nerr
  logical newbst
  logical newopt
  integer nvars
  real ( kind = 8 ), parameter :: one = 1.0D+00
  logical passb
  real ( kind = 8 ) pb
  real ( kind = 8 ) pd
  real ( kind = 8 ) pv
  real ( kind = 8 ) rb
  real ( kind = 8 ) rdum
  logical retrea
  real ( kind = 8 ) rg
  real ( kind = 8 ) rnormc
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) semibg
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  logical term
  real ( kind = 8 ) told
  real ( kind = 8 ) tolf
  real ( kind = 8 ) tolp
  real ( kind = 8 ) tolsnr
  real ( kind = 8 ) tolx
  real ( kind = 8 ) wa(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xb(*)
  character ( len = 128 ) xmess

  save
!
!     OPTIONS:
!
!     1    SET THE PRINTED OUTPUT OFF/ON.  REQUIRES TWO ENTRIES
!          IN IOPT(*).  IF IOPT(*+1)=0, NO PRINTING; =1 PRINT.
!          (THIS PRINTOUT SHOWS VARIOUS QUANTITIES ASSOCIATED
!          WITH THE NONLINEAR PROBLEM BEING SOLVED. GOOD FOR
!          DEBUGGING A PROBLEM IN CASE OF TROUBLE.
!
!     2    SET THE MAXIMUM NUMBER OF INTERATIONS.  REQUIRES TWO ENTRIES
!          IN IOPT(*).  USE IOPT(*+1)= MAXIMUM NUMBER OF ITERATIONS.
!
!     3    PASS INITIAL BOUNDS FOR THE TRUST REGION TO THE NONLINEAR
!          SOLVER. REQUIRES TWO ENTRIES IN IOPT(*). USE IOPT(*+1) AS A
!          POINTER INTO ROPT(*) FOR START OF THE NVARS BOUNDS.
!
  rdum = 0.0D+00

  go to (50),iflag
!
!  PROCESS OPTION ARRAY
!
  assign 10 to igopoa
  go to 470
!
   10 continue
!
!  INITIALIZE OTHER VALUES
!
  assign 20 to igoiov
  go to 450

   20 continue
!
!  SET SO X(*)-DX(*) UPDATE WILL BE CORRECT FIRST TIME.
!
  dx(1:nvars) = 0.0D+00
  k = 0
!
!  FB = "INFINITY" ON THIS MACHINE.
!
  fb = huge ( fb )
  dxnrm = fb
  fl = 0.0D+00
!
!  LINEAR PROBLEM RESIDUAL.
!
  pv = 0.0D+00
  retrea = .false.
  fulnwt = .false.
  term = .false.

   30 continue

  if ( .not. retrea) iters = iters + 1
  if ( retrea) then
!
!  MUST RETREAT TO BEST X VALUES.
!
      k = 0
      kl = -1
      fl = fb
      call dcopy(nvars,xb,1,x,1)
  else
      kl = k
      x(1:nvars) = x(1:nvars) - dx(1:nvars)
      if ( term) then
          iflag = 0
          go to 390
      end if
  end if

  igo = 1
  iflag = 1
  go to 440

   50 continue
  fc = dnrm2(mequa,fjac(mcon+1,nvars+1),1)
!
!  TEST FOR CONVERGENCE
!
  assign 60 to igotfc
  go to 400

   60 continue

  if ( term) then
      iflag = 0
      go to 390
  end if

  newbst = fc  <  fb .or. (mcon > 0 .and.iters == 2)
  if ( newbst) k = 0
  if ( k == 0) then
      rg = 0.0D+00
      pb = 0.0D+00
      pd = huge ( pd )
!
!  WANT TO POSITION AT BEST X VALUES.
!
      fb = fc

      go to (70,90,110),2 - kl

      go to 120

   70     continue
!
!  IMMEDIATELY GOT A NEW BEST X.
!
      if ( t2 <= 0.25D+00 ) then
          bboost = 1.0D+00
          chg = max ( 4.0D+00 * t2, 0.1D+00 )
      end if

      do j = 1,nvars
         bb(j) = chg*bb(j)
      end do
!
!  THIS CODE FOR ALPHA HAS THE FOLLOWING EFFECT.
!  IF FC .EQ. PV, THEN ALPHA=ALFAC.
!  IF FC**2 .EQ. PV*FL THEN ALPHA=2.-1./ALFAC
!  IF FC**2 IS MUCH LARGER THAN PV*FL, THEN ALPHA=1.
!
      t = fc - pv
      if ( t == 0.0D+00 ) then
          alpha = alfac
      else
          alpha = (pv* (fl-pv))/ (fc+pv)/ (alfac-one)
          alpha = (abs(t)+alfac*alpha)/ (abs(t)+alpha)
      end if

      alfac = 1.5D+00 * alpha
      bboost = min(1.5D+00*alpha*bboost,semibg)
      go to 140

   90     continue
!
!  AT THE INITIAL X.
!
      alfac = 256.0D+00

      do j = 1,nvars
        if ( .not. passb) bb(j) = -x(j)
        if ( bb(j) == 0.0D+00 ) then
          colnrm = dnrm2(mequa,fjac(mcon+1,j),1)
          if ( colnrm/= 0.0D+00 ) bb(j) = -fc/colnrm
        end if

        if ( bb(j) == 0.0D+00 ) bb(j) = -1.0D+00
        xb(j) = x(j)
        b(j) = bb(j)
      end do

      alpha = 1.0D+00
      bboost = 0.5D+00
      go to 170

  110     continue
!
!  RETREAT TO BEST X.
!
      if ( alfac /= 256.0D+00 ) then
          alpha = min ( 1.0D+00 / alfac, 0.25D+00 )
          alfac = 1.25D+00
      else
          alpha = 0.25D+00*alpha
      end if

      bboost = 0.25D+00
      go to 140

  120     continue
!
!  NOT IMMEDIATELY A BEST X.
!
      rb = 0.0D+00
      do j = 1,nvars
         rb = max(rb,abs((xb(j)-x(j))/bb(j)))
      end do
      alpha = rb
      alfac = 2.0D+00
      bboost = ( 8.0D+00 / 7.0D+00 + rg )/ ( 2.0D+00 / 7.0D+00 + rg )

  140     continue

      do j = 1,nvars
         dx(j) = xb(j) - x(j)
         if ( dx(j) == 0.0D+00 ) then
             b(j) = alpha*bb(j)
         else
             xb(j) = x(j)
             b(j) = sign(alpha*bb(j),dx(j)) + bboost*dx(j)
         end if

        bb(j) = sign(min(sqrt ( huge ( bb(j) ) ),abs(b(j))),b(j))
     end do

  else
!
!  COMPUTE A GAUGE FOR RETREATING IF DID NOT GET A NEW BEST.
!
      if ( k == 1) then
          pb = pv
          pd = 1.5D+00 * (fb+pb* (pb/fb)) - 4.0D+00 * pb
      end if

      alpha = ( 0.5D+00 * fc+fl)/ (fc+fl)
      chg = min(alpha*chg,t2)
      chg = max ( chg, 0.1D+00 )

      do j = 1,nvars
         b(j) = chg*b(j)
         if ( k == 1) bb(j) = b(j)
      end do

  end if

  170 continue
!
!  TEST FOR CONVERGENCE
!
  assign 180 to igotfc
  go to 400

  180 continue

  if ( term) then
      iflag = 0
      go to 390
  end if

  k = k + 1
!
!  SOLVE LINEAR BOUNDED PROBLEM.
!
  do 240 j = 1,nvars

     if ( b(j) < 0.0D+00 ) then
         alb = b(j)
         if ( dx(j) == 0.0D+00 ) then
!
!  THIS CASE IS REQD. TO AVOID USING BUB(*) AT THE INITIAL PT.
!
             aub = -c1516*alb
         else
             aub = min(-c1516*alb,-dx(j)+bub(j))
         end if
     else
         aub = b(j)

         if ( dx(j) == 0.0D+00 ) then
             alb = -c1516*aub
         else
             alb = max(-c1516*aub,-dx(j)+blb(j))
         end if

     end if
!
!  RESTRICT THE STEP FURTHER IF USER GIVES BOUNDS.
!
     icase = ind(j)
     go to (190,200,210,220),icase

  190    continue
     aub = min(aub,x(j)-bl(j))
     go to 230
  200    continue
     alb = max(alb,x(j)-bu(j))
     go to 230
  210    continue
     aub = min(aub,x(j)-bl(j))
     alb = max(alb,x(j)-bu(j))
     go to 230
  220    continue
  230    continue
     blb(j) = alb
!
!  THIS NEXT LINE IS TO GUARANTEE THAT THE LOWER BOUND
!  IS .LE. THE UPPER BOUND.
!
     aub = max(aub,alb)
     bub(j) = aub
     indb(j) = 3

  240 continue
!
!  SEE IF USER HAS GIVEN GENERAL CONSTRAINTS.
!
  do j = nvars + 1,nall

     icase = ind(j)
     gval = fjac(j-nvars,nvars+1)

     go to (250,260,270,280),icase

  250    continue

     blb(j) = - (gval-bl(j))
     indb(j) = 1
     go to 290

  260    continue

     bub(j) = - (gval-bu(j))
     indb(j) = 2
     go to 290

  270    continue

     blb(j) = - (gval-bl(j))
     bub(j) = - (gval-bu(j))
     indb(j) = 3
     go to 290

  280    continue

     indb(j) = 4

  290    continue

  end do
!
!  SOLVE THE LEAST SQUARES PROBLEM WITH BOUNDS AND LINEAR
!  CONSTRAINTS.  THESE BOUNDS CAN COME FROM THE USER OR
!  THE ALGORITHM.
!
  call dbocls(fjac,ldfjac,mcon,mequa,nvars,blb,bub,indb,iopt(ipls), &
    dx,rnormc,pv,mode,wa,iwa)

  if ( iprint > 0) then
      write(*,9011) iters,fc,pv,k,kl,fb,alpha,bboost
      write(*,9001) '  +x=', (x(j),j=1,nvars)
      write(*,9001) ' +xb=', (xb(j),j=1,nvars)
      write(*,9001) ' +dx=', (dx(j),j=1,nall)
      write(*,9001) ' + b=', (b(j),j=1,nvars)
      write(*,9001) ' +lb=', (blb(j),j=1,nall)
      write(*,9001) ' +ub=', (bub(j),j=1,nall)
      write(*,'('' +///end of iteration.///'')')
  end if
!
!  TEST FOR NOISE IN LINEAR PROBLEM SOLN.
!
  term = ( mcon == 0 .and. (pv>=fc) )
  term=.false.
  if ( term) then
      if ( iprint > 0) then
          write(*,9021) pv,fc
      end if

      call dcopy(nvars,xb,1,x,1)
      igo = 4
      iflag = 0
      go to 390
  end if

  rg = max(rg, (pv-pb)/pd)
  if ( .not. retrea) then
      chg = one
      t2 = 0.0D+00

      do j = 1,nvars

         bold = b(j)
         t = dx(j)/bold
!
!  IF USER GIVES BOUNDS, AND THESE BOUNDS ARE HIT,
!  DO NOT DETRACT FROM DECLARING A FULL NEWTON STEP.
!
         icase = ind(j)
         go to (310,320,330,340),icase

  310        alb = (x(j)-bl(j))/bold
         aub = -semibg
         go to 350

  320        aub = (x(j)-bu(j))/bold
         alb = -semibg
         go to 350

  330        alb = (x(j)-bl(j))/bold
         aub = (x(j)-bu(j))/bold
         go to 350

  340        continue
         alb = -semibg
         aub = -semibg

  350        continue

         if ( t == 1.0D+00 ) then
             t2 = 1.0D+00
             b(j) = bold + bold
             chg = chg*chgfac
         else
             if ( abs(t) < 0.25D+00 .and. dx(j)/= 0.0D+00 ) then
                 b(j) = sign(0.25D+00 * bold,dx(j)) + 3.0D+00 * dx(j)
             else
                 b(j) = sign(bold,dx(j))
             end if
         end if
!
!  THIS TEST AVOIDS THE USER BOUNDS IN DECLARING A NEWTON STEP.
!
         if ( abs(alb-t)>=0.01D+00*abs(t) .and. &
              abs(aub-t) >= 0.01D+00*abs(t)) then
             if ( t > 0.0D+00 ) then
                 t2 = max(t2,t)
             else
                 t2 = max(t2,-t/c1516)
             end if
         end if

      end do

      fulnwt = t2  <  0.99D+00
      fl = fc
      dxnrm = abs(dx(idamax(nvars,dx,1)))
!
!  TEST FOR SMALL ABSOLUTE CHANGE IN X VALUES.
!
      term = dxnrm  <  told .and. fulnwt
      if ( term) then
          igo = 5

          go to 370

      end if

      term = dxnrm  <  dnrm2(nvars,x,1)*tolx .and. fulnwt
      term = term .and. (iters > 1)
      if ( term) then
          igo = 6

          go to 370

      end if

      go to 380

  370     continue
      go to 30

  end if

  380 continue
  go to 30

  390 continue
  go to 440
!
!  TEST FOR CONVERGENCE
!
  400 continue
!
!  TEST FOR SMALL FUNCTION NORM.
!
  term = fc <= tolf .or. term
!
!  IF HAVE CONSTRAINTS MUST ALLOW AT LEAST ONE MOVE.
!
  term = term .and. (mcon == 0 .or. iters > 1)
  if ( term) then
      igo = 2
      go to 420
  end if
!
!  TEST FOR NO CHANGE
!
  assign 410 to igotnc
  go to 430

  410 continue
  term = term .and. .not. retrea
  if ( term) then
      igo = 3
      go to 420
  end if

  term = iters >= itmax
  if ( term) then
      igo = 7
  end if

  420 continue
  go to igotfc
!
!  TEST FOR NO CHANGE
!
  430 continue
  t = sqrt(max( 0.0D+00, (fl-pv)* (fl+pv)))
  term = (abs(fb-fc)<=tolsnr*fb) .and. (t.le.pv*tolp)
  term = term .and. (abs(fc-fl)<=fb*tolsnr)
  term = term .and. fulnwt

  go to igotnc

  440 continue
  return
!
!  INITIALIZE OTHER VALUES
!
  450 continue

  iters = 0
  nall = mcon + nvars
  chgfac = 2.0D+00** (-one/ real ( nvars, kind = 8 ))
  c1516 = 15.0D+00 / 16.0D+00
  semibg = 1.0D+10
!
!  MAKE SURE THAT VARIABLES SATISFY THE BOUNDS AND CONSTRAINTS.
!
  do j = 1,nall
     blb(j) = bl(j)
     bub(j) = bu(j)
     indb(j) = ind(j)
  end do

  go to igoiov
!
!  PROCESS OPTION ARRAY
!
  470 continue
  iprint = 0
!
!  T = MACHINE REL. PREC.
!
  t = epsilon ( t )
  tolf = t
  tolx = tolf
  told = tolf
  tolsnr = 1.0D-03
  tolp = 1.0D-03
  itmax = 18
  passb = .false.
  level = 1
  ipls = 0
  lpdiff = 0
  lp = 1

  480 continue

  lp = lp + lpdiff
  lpdiff = 2
  kp = iopt(lp)
  newopt = kp  >  0
  jp = abs(kp)
!
!  SEE IF THIS IS THE LAST OPTION.
!
  if ( jp == 99) then
      if ( newopt) then
!
!  THE POINTER TO THE START OF OPTIONS FOR THE LINEAR
!  SOLVER MUST SATISFY THE REQUIREMENTS FOR THAT OPTION ARRAY.
!
          if ( ipls == 0) ipls = lp
          go to 490
      else
          lpdiff = 1
          go to 480
      end if
  end if
!
!  CHANGE PRINT OPTION.
!
  if ( jp == 1) then
      if ( newopt) iprint = iopt(lp+1)
      go to 480
  end if
!
!  SEE IF MAX. NUMBER OF ITERATIONS CHANGING.
!
  if ( jp == 2) then
      if ( newopt) itmax = iopt(lp+1)
      go to 480
  end if
!
!  SEE IF BOUNDS FOR THE TRUST REGION ARE BEING PASSED.
!
  if ( jp == 3) then
      if ( newopt) then
          call dcopy(nvars,ropt(iopt(lp+1)),1,bb,1)
          passb = .true.
      end if

      go to 480

  end if
!
!  CHANGE TOLERANCE ON THE LENGTH OF THE RESIDUALS.
!
  if ( jp == 4) then
      if ( newopt) tolf = ropt(iopt(lp+1))
      go to 480
  end if
!
!  CHANGE TOLERANCE ON THE NORM OF THE RELATIVE
!  CHANGE TO THE PARAMETERS.
!
  if ( jp == 5) then
      if ( newopt) tolx = ropt(iopt(lp+1))
      go to 480
  end if
!
!  CHANGE TOLERANCE ON ABSOLUTE CHANGE TO THE PARAMETERS.
!
  if ( jp == 6) then
      if ( newopt) told = ropt(iopt(lp+1))
      go to 480
  end if

  if ( jp == 7) then
!
!  CHANGE TOLERANCE FOR RELATIVE AGREEMENT BETWEEN
!  BEST FUNCTION NORM, LAST FUNCTION NORM AND THE
!  CURRENT FUNCTION NORM.
!
      if ( newopt) tolsnr = ropt(iopt(lp+1))
      go to 480
  end if

  if ( jp == 8) then
!
!  CHANGE TOLERANCE FOR AGREEMENT BETWEEN PREDICTED
!  VALUE OF RESIDUAL NORM AND THE PREVIOUS VALUE OF
!  THE RESIDUAL NORM.
!
      if ( newopt) tolp = ropt(iopt(lp+1))
      go to 480
  end if
!
!  CHANGE THE PRINT LEVEL IN THE ERROR PROCESSOR.
!
  if ( jp == 9) then
      if ( newopt) level = iopt(lp+1)
      go to 480

  end if
!
!  PASS AN OPTION ARRAY TO THE CONSTRAINED LINEAR SOLVER.
!  THIS OPTION IS A POINTER TO THE START OF THE OPTION
!  ARRAY FOR THE SUBPROGRAM.
!
  if ( jp == 10) then
      if ( newopt) ipls = iopt(lp+1)
      go to 480
  end if
!
!  MOVE THE PROCESSING POINTER BY THE VALUE IN THE
!  NEXT ENTRY OF THE OPTION ARRAY.  THIS DEVICE IS
!  INCLUDED SO THAT PASSING OPTIONS TO LOWER LEVEL
!  SUBROUTINES IS EASY TO DO.
!
  if ( jp == 11) then
      if ( newopt) lpdiff = iopt(lp+1)
      go to 480
  end if
!
!  SAW AN OPTION (OR GARBAGE) THAT IS NOT ON THE LIST.
!
  xmess = 'dqedip. invalid option processed. i1=iopt(*) entry. i2=iopt(i1).'
  nerr = 08
  igo = 16
  call xerrwv(xmess,nerr,level,2,lp,iopt(lp),0,rdum,rdum)
  iflag = 0

  go to 440

  490 continue
  go to igopoa

 9001 format (a4,1p10d12.4/ (4x,10d12.4))
 9011 format ('0+iter.=',i3,' fc=',1pd10.4,' pv=',1pd10.4,' k=',i4, &
   ' kl=',i4,' fb=',1pd10.4/12x,'al=',1pd14.4,' bb=',1pd14.4)
 9021 format (' linear residual>=current f. quitting.',1p2d12.5)
end
subroutine dqedmn ( dqedev, mequa, nvars, mcon, ind, bl, bu, x, fjac, &
  ldfjac, fb, igo, iopt, ropt, iwa, wa, dx, xb, b, bb, blb, bub, indb, &
  npmax, zp, xp, qc, mdqc, pj, wj, ldwj, gr, dxl )

!***********************************************************************
!
!! DQEDMN is the main solution routine called by DQED.
!
!***BEGIN PROLOGUE  DQEDMN
!***REFER TO  DQED
!***END PROLOGUE  DQEDMN
!
  implicit none

  integer ldfjac
  integer mdqc
  integer npmax
  integer nvars

  real ( kind = 8 ) ajn
  real ( kind = 8 ) alb
  real ( kind = 8 ) alfac
  real ( kind = 8 ) alpha
  real ( kind = 8 ) aub
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) bb(*)
  real ( kind = 8 ) bboost
  real ( kind = 8 ) bl(*)
  real ( kind = 8 ) blb(*)
  real ( kind = 8 ) bold
  real ( kind = 8 ) bu(*)
  real ( kind = 8 ) bub(*)
  real ( kind = 8 ) c1516
  real ( kind = 8 ) chg
  real ( kind = 8 ) chgfac
  real ( kind = 8 ) colnrm
  real ( kind = 8 ) cond
  real ( kind = 8 ) cosl
  real ( kind = 8 ) cosm
  real ( kind = 8 ) cosq
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dfn
  real ( kind = 8 ) dnrm2
  external dqedev
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dxl(*)
  real ( kind = 8 ) dxnrm
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) fjac(ldfjac,*)
  real ( kind = 8 ) fl
  logical fulnwt
  real ( kind = 8 ) gr(*)
  integer i
  integer icase
  integer idamax
  integer :: iflag = 0
  integer igo
  integer igoelm
  integer igoeqm
  integer igotfc
  integer igow
  integer ind(*)
  integer indb(*)
  integer iopt(*)
  integer ipls
  integer iprint
  integer iters
  integer itmax
  integer iwa(*)
  integer j
  logical jactri
  integer jk
  integer jp
  integer k
  integer kl
  integer kp
  integer l
  integer ldwj
  integer level
  logical linmod
  integer lk
  integer lp
  integer lpdiff
  integer mcon
  integer mconst
  integer me
  integer mequa
  integer mk
  logical mustcn
  integer nall
  integer nerr
  logical newbst
  logical newopt
  integer nit
  logical noquad
  integer np
  integer nt
  integer ntterm
  integer nv
  real ( kind = 8 ), parameter :: one = 1.0D+00
  logical passb
  real ( kind = 8 ) pb
  real ( kind = 8 ) pd
  real ( kind = 8 ) pj(*)
  real ( kind = 8 ) pv
  real ( kind = 8 ) pvl
  real ( kind = 8 ) qc(mdqc,npmax)
  real ( kind = 8 ) rb
  real ( kind = 8 ) rc
  real ( kind = 8 ) rcond
  real ( kind = 8 ) rdum
  logical retrea
  logical revers
  real ( kind = 8 ) rg
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) sc
  real ( kind = 8 ) semibg
  real ( kind = 8 ) ss
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  logical term
  real ( kind = 8 ) told
  real ( kind = 8 ) tolf
  real ( kind = 8 ) tolp
  real ( kind = 8 ) tolsnr
  real ( kind = 8 ) toluse
  real ( kind = 8 ) tolx
  real ( kind = 8 ) tt
  logical useq
  logical useql
  real ( kind = 8 ) wa(*)
  real ( kind = 8 ) wj(ldwj,*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xb(*)
  character xmess*128
  real ( kind = 8 ) xp(nvars,npmax)
  real ( kind = 8 ) zn
  real ( kind = 8 ) zp(mequa,npmax)

  save

  rdum = 0.0D+00

  if ( iflag /= 0 ) then
    go to 50
  end if

  lk = min ( mequa, nvars+1 )
  nt = min ( nvars+1, mequa-1 )
!
!  PROCESS OPTION ARRAY
!
  go to 1100
!
   10 continue
!
!  INITIALIZE OTHER VALUES
!
  go to 1030

20 continue
!
!  SET SO X(*)-DX(*) UPDATE WILL BE CORRECT FIRST TIME.
!
  dx(1:nvars) = 0.0D+00
  k = 0
!
!  Set "INFINITY" ON THIS MACHINE.
!
  fb = huge ( fb )
  dxnrm = fb
  fl = 0.0D+00
!
!  MODEL PROBLEM RESIDUAL.
!
  pv = 0.0D+00
  pvl = 0.0D+00
  retrea = .false.
  fulnwt = .false.
  term = .false.

30 continue

  iters = iters + 1

  if ( retrea ) then
    x(1:nvars) = xb(1:nvars)
    k = 0
    kl = -1
    fl = fb
  else
    kl = k
    x(1:nvars) = x(1:nvars) - dx(1:nvars)
  end if

  if ( term) then
    iflag = 0
    go to 840
  end if

  iflag = 1
  igo = 1
  if ( np == npmax-1 .and. np < nvars ) then
    igo = -1
  end if
!
!  THERE ARE TWO POSSIBLE WAYS TO GET FUNCTION AND DERIVATIVE
!  VALUES FROM THE USER.  THE OPTIONAL WAY IS REVERSE COMMUNICATION.
!  THE NOMINAL WAY IS THROUGH FORWARD COMMUNICATION.
!
  if ( revers) then
    go to 1020
  else
    call dqedev ( x, fjac, ldfjac, igo, iopt, ropt )
  end if

   50 continue
!
!  IF IGO HAS BEEN CHANGED BY THE USER TO A VALUE .GT. 1, THEN
!  THIS IS AN ABORT SIGNAL.  STOP UNLESS IT = 99.
!
  if ( igo == 99) then
!
!  IF IGO = 99 THE EVALUATION CAN'T BE PERFORMED.
!  WE FORCE A RETREAT AND RESTART IN THIS CASE.
!
      do i = mcon + 1, mcon + mequa
         fjac(i,nvars+1) = fc
         fjac(i,1:nvars) = 0.0D+00
      end do
!
!  A RETREAT IS FORCED TO OCCUR WITH THIS ASSIGNMENT.
!
      retrea = .true.

  end if

  fc = dnrm2(mequa,fjac(mcon+1,nvars+1),1)

  if ( igo > 1 .and. igo /= 99 ) then
    iflag = 0
    x(1:nvars) = xb(1:nvars)
    go to 1020
  end if
!
!  SAVE PAST FUNCTION AND VARIABLE VALUES.
!  DO NOT UPDATE THE PAST POINTS UNLESS THERE IS A
!  SIGNIFICANT CHANGE IN THE X(*) VALUES.
!
  if ( np >= 0 ) then
      if ( dxnrm > toluse*dnrm2(nvars,x,1)) then
          lp = nvars
          if ( .not. noquad) np = min(np,npmax-1,lp) + 1
          do j = np - 1,1,-1
             call dcopy(nvars,xp(1,j),1,xp(1,j+1),1)
             call dcopy(mequa,zp(1,j),1,zp(1,j+1),1)
          end do
      end if
  end if
!
!  PUT IN THE PRESENT VALUES OF THE VARIABLES and functions.
!
  xp(1:nvars,1) = x(1:nvars)

  call dcopy(mequa,fjac(mcon+1,nvars+1),1,zp(1,1),1)
!
!  THIS STATEMENT HAS THE EFFECT OF A FIRST TIME FLAG.
!
  np = max ( np, 0 )
!
!  COMPUTE THE COSINES OF THE PAST MOVES WITH THE MOST CURRENT MOVE.
!
  do l = 2,np
     qc(1:nvars,l) = xp(1:nvars,l) - xp(1:nvars,1)
  end do

  l = 3

  do while ( l <= np )
!
!  CALCULATE THE DIRECTION COSINES OF THE PAST MOVES.
!
    t = dot_product ( qc(1:nvars,2), qc(1:nvars,l) )

    tt = dnrm2 ( nvars, qc(1,2), 1 ) * dnrm2 ( nvars, qc(1,l), 1 )

    if ( tt > 0.0D+00 ) then
      t = t / tt
    else
      t = 1.0D+00
    end if

    if ( iprint > 0 ) then
      write (*, &
        '('' past move number, cosine of move'',i3,2x,f6.2)') l - 2, t
    end if

    if ( abs ( t ) > 0.98D+00 ) then
!
!  DISCARD PAST INFORMATION ASSOCIATED WITH THIS MOVE IF CLOSE TO
!  A PAST MOVE.
!
      do j = l,np - 1
        call dcopy(mequa,zp(1,j+1),1,zp(1,j),1)
        call dcopy(nvars,xp(1,j+1),1,xp(1,j),1)
        call dcopy(nvars,qc(1,j+1),1,qc(1,j),1)
      end do

      np = np - 1
      cycle
    end if

    l = l + 1

  end do
!
!  COMPUTE FUNCTION DIFFERENCES IN QC.
!
  do j = 1,np - 1
     do i = 1,mequa
        qc(i,j+1) = zp(i,j+1) - zp(i,1)
     end do
  end do
!
!  NOW HAVE F(PAST)-F(CURRENT) IN QC( , ), COLS. 2,...,NP USED.
!  COMPUTE NORM OF DIFFERENCE OF FUNCTION VALUES.
!
  if ( np > 1) then
    dfn = dnrm2(mequa,qc(1,2),1)
  else
    dfn = 0.0D+00
  end if
!
!  NEXT ADD PRODUCT OF JACOBIAN AND PAST VARIABLE DIFFERENCES.
!
  do i = 1,np - 1
     do j = 1,nvars
        call daxpy(mequa,- (xp(j,i+1)-xp(j,1)),fjac(mcon+1,j),1,qc(1,i+1),1)
     end do
  end do

  250 continue
!
!  COMPUTE THE SYMMETRIC MATRIX WHOSE ENTRIES ARE THE
!  SQUARES OF THE DOT PRODUCTS OF THE PAST DIRECTIONS.
!  THIS MATRIX IS REQUIRED TO OBTAIN THE QUADRATIC TERMS
!  ASSOCIATED WITH INTERPOLATING TO PAST FUNCTION VALUES.
!
  do l = 2,np
    do j = l,np
      t = 0.0D+00
      do i = 1,nvars
        t = t + (xp(i,j)-xp(i,1)) * (xp(i,l)-xp(i,1))
      end do
      wj(l-1,j-1) = t
      wj(j-1,l-1) = t
    end do
  end do
!
!  COMPUTE NORM OF REMAINDER INCLUDING LINEAR TERMS,
!  USING THE LAST MOVE.
!
  useq = np  >  1 .and. .not. retrea

  zn = 1.0D+00

  if ( np > 1) then
      zn = dnrm2(mequa,qc(1,2),1)
!
!  COMPUTE RATIO OF Z TERMS TO CURRENT F VALUE..
!
      if ( useq) then
        useq = (zn  > 1.0D-04*dfn .and. zn  <  dfn*0.75D+00 ) .or. useql
      end if

      if ( dfn > 0.0D+00 ) then
        zn = zn / dfn
      else
        zn = 1.0D+00
      end if

      if ( iprint > 0) then
        call dvout(1,zn,'('' ratio of z term to past df norm'')',4)
      end if
!
!  SCALE THE MATRIX (MATRIX := D*MATRIX*D, WHERE D**2 = RECIPROCAL
!  OF THE DIAGONAL TERMS OF THE MATRIX.
!
      do i = 1,np - 1
         dxl(i) = wj(i,i)
         if ( dxl(i) == 0.0D+00 ) then
             np = i
             go to 250
         else
             dxl(i) = 1.0D+00 / dxl(i)
         end if
      end do

      do i = 1,np - 1
         do j = 1,np - 1
            wj(i,j) = (wj(i,j)*dxl(i))* (wj(i,j)*dxl(j))
         end do
      end do
!
!  USE THE LINPACK ROUTINES DGECO(), DGESL() TO OBTAIN
!  THE COEFFICIENTS OF THE QUADRATIC TERMS, ONE ROW AT A TIME.
!
      call dgeco(wj,ldwj,np-1,iwa,rcond,wa)

      if ( iprint > 0) then
        write(*,'('' rcond from dgeco() = '',2x,       1pd15.4)') rcond
      end if

      if ( cond*rcond < 1.0D+00 ) then
        np = np - 1
        go to 250
      end if
!
!  COPY A ROW OF THE INTERPOLATED DATA TO A WORKING ARRAY.
!  USE THIS ARRAY TO OBTAIN A ROW OF THE QUADRATIC TERMS.
!
!  SCALE THE RIGHT HAND SIDE DATA.
!
!  RESCALE THE SOLUTION DATA.
!
!  THE SIGN CHANGE COMES FROM A CHANGE
!  OF SIGN IN THE INNER LOOP MODEL PROBLEM.
!
      do i = 1,mequa

         call dcopy(np-1,qc(i,2),mdqc,wa,1)

         do j = 1,np - 1
           wa(j) = wa(j)*dxl(j)
         end do

         call dgesl(wj,ldwj,np-1,iwa,wa,0)

         do j = 1,np - 1
           wa(j) = -2.0D+00 * wa(j)*dxl(j)
         end do

         call dcopy(np-1,wa,1,qc(i,2),mdqc)

      end do

  end if

  350 continue
!
!  NOW HAVE THE QUADRATIC TERMS COMPUTED.
!  NEXT WILL TRIANGULARIZE THE JACOBIAN TO SAVE SPACE
!  WHEN USING THE QUADRATIC MODEL.
!
!  CONSTRUCT AND THEN APPLY PLANE ROTATIONS
!  TO ACHIEVE UPPER TRIANGULAR FORM.  THIS LOOP
!  AFFECTS THE JACOBIAN AND RIGHT HAND SIDE.
!
!  APPLY THE TRANSFORMATION TO THE QUADRATIC TERMS.
!
  if ( jactri) then

    do j = 1,nt
      do i = j + 1,mequa
        call drotg(fjac(mcon+j,j),fjac(mcon+i,j),sc,ss)
        call drot(nvars-j+1,fjac(mcon+j,j+1),ldfjac,fjac(mcon+i,j+1), &
          ldfjac,sc,ss)
        call drot(np-1,qc(j,2),mdqc,qc(i,2),mdqc,sc,ss)
        fjac(mcon+i,j) = 0.0D+00
      end do
    end do
!
!  NOW WE FINISH TRIANGULARIZING THE QUADRATIC TERMS.
!  NOTE THAT THIS DOES NOT AFFECT THE RIGHT HAND SIDE.
!
    do l = 1,np - 1
      do i=nvars+l+2,mequa
         call drotg(qc(nvars+l+1,l+1),qc(i,l+1),sc,ss)
         call drot(np-l-1,qc(nvars+l+1,min(l+2,npmax)),mdqc, &
                          qc(i,min(l+2,npmax)),mdqc,sc,ss)
         qc(i,l+1) = 0.0D+00
      end do
    end do

  end if
!
!  COMPUTE CURRENT NORM OF J**T*F(X).
!
  do j = 1,nvars

     if ( jactri) then
       jk = j
     else
       jk = mequa
     end if

     pj(j) = ddot ( jk, fjac(mcon+1,j), 1, fjac(mcon+1,nvars+1), 1 )

  end do

  ajn = dnrm2(nvars,pj,1)
!
!  SAVE J**T*F FOR DIRECTION TESTING WITH LINEAR AND QUADRATIC MOVES.
!
  if ( ajn > 0.0D+00 ) then
    call dscal(nvars,one/ajn,pj,1)
  end if

  call dcopy(nvars,pj,1,gr,1)
  newbst = fc  <  fb
  if ( newbst) k = 0
!
!  WANT TO POSITION AT BEST X VALUES.
!
  if ( k == 0) then

      pb = 0.0D+00
      pd = huge ( pd )

      if ( .not. retrea) then
          fb = fc
      end if

      go to (410,430,450),2 - kl

      go to 470

  410     continue
!
!  IMMEDIATELY GOT A NEW BEST X.
!
      rg = 0.0D+00
      if ( t2 <= 0.25D+00 ) then
          bboost = one
          chg = max ( 4.0D+00 * t2, 0.1D+00 )
      end if

      do j = 1,nvars
         bb(j) = chg*bb(j)
      end do
      t = 0.25D+00 / (alfac-1.0D+00 )
      alpha = (zn+alfac*t)/ (zn+t)
      alfac = 1.5D+00 * alpha
      bboost = min(1.5D+00 *alpha*bboost,semibg)
      go to 490

  430     continue
      useql = .false.
      rg = 0.0D+00
!
!  AT THE INITIAL X.
!
      alfac = 256.0D+00

      do j = 1,nvars

         if ( .not. passb) then
             bb(j) = -x(j)
         end if

         if ( bb(j) == 0.0D+00 ) then

             if ( jactri) then
                 jk = j
             else
                 jk = mequa
             end if

             colnrm = dnrm2(jk,fjac(mcon+1,j),1)

             if ( colnrm/= 0.0D+00 ) then
                 bb(j) = ddot(jk,fjac(mcon+1,j),1,fjac(mcon+1,nvars+1),1)
                 bb(j) = -max(abs(bb(j))/colnrm/colnrm,fc/colnrm)
             else
                 bb(j) = -one
             end if

         end if

         xb(j) = x(j)
         b(j) = bb(j)

      end do

      alpha = one
      bboost = 0.5D+00
      go to 520

  450     continue
!
!  RETREAT TO BEST X.
!
      if ( alfac /= 256.0D+00 ) then
        alpha = min ( 4.0D+00 / alfac, 0.25D+00 )
        alfac = 1.25D+00
      else
        alpha = 0.25D+00 * alpha
      end if

      bboost = 0.25D+00
      useql = .false.

      do j = 1,nvars
!
!  THE NEXT LINES HAVE THE EFFECT OF REDUCING THE BOUNDS
!  AT THE CURRENT BEST TO ABOUT 1/16 THEIR CURRENT VALUES
!  PROVIDED THE CURRENT BOUNDS ARE RELATIVELY SMALL.  IF
!  THE CURRENT BOUNDS ARE RELATIVELY LARGE, THE BOUNDS AT
!  THE BEST ARE LEFT ABOUT THE SAME.
!
         t = abs(b(j))
         tt = abs(bb(j))
         t = ( t + 0.25D+00 * tt ) / ( t + 4.0D+00 * tt )
         b(j) = t*bb(j)
         bb(j) = b(j)
         dx(j) = 0.0D+00

      end do

      go to 520

  470     continue
!
!  NOT IMMEDIATELY A BEST X.
!
      rb = 0.0D+00
      do j = 1,nvars
         rb = 0.125D+00 * max(rb,abs((xb(j)-x(j))/bb(j)))
      end do

      alpha = rb
      alfac = 2.0D+00
      bboost = ( 1.0D+00 + rg )/ ( 0.25D+00 + 2.0D+00 * rg )

  490     continue

      do j = 1,nvars
         dx(j) = xb(j) - x(j)
         if ( dx(j) == 0.0D+00 ) then
             b(j) = alpha*bb(j)
         else
             xb(j) = x(j)
             b(j) = sign(alpha*bb(j),dx(j)) + bboost*dx(j)
         end if
         bb(j) = sign(min(sqrt( huge ( bb(j) ) ),abs(b(j))),b(j))
      end do

  else
!
!  MUST MAKE SURE THAT PD GETS SET TO A REASONABLE VALUE.
!  COMPUTE A GAUGE FOR RETREATING IF DID NOT GET A NEW BEST.
!
      alpha = ( 0.5D+00 * zn + 1.0D+00 )/ ( zn + 1.0D+00 )
      chg = alpha*chg

      if ( k == 1) then
          chg = min(chg,t2)
          chg = max(chg, 0.1D+00 )
          pb = pv
          pd = 1.5D+00 * (fb+pb* (pb/fb)) - 4.0D+00 *pb
      end if

      do j = 1,nvars
         b(j) = chg*b(j)
         if ( k == 1) bb(j) = b(j)
      end do

  end if

  520 continue
!
!  TEST FOR CONVERGENCE
!
  assign 530 to igotfc
  go to 980

  530 continue

  if ( term) then
      iflag = 0
      go to 840
  end if

  k = k + 1
!
!  SOLVE MODEL BOUNDED PROBLEM.
!
  do j = 1,nvars

     if ( b(j) < 0.0D+00 ) then
         alb = b(j)
         if ( dx(j) == 0.0D+00 ) then
!
!  THIS CASE IS REQD. TO AVOID USING BUB(*) AT THE INITIAL PT.
!
             aub = -c1516*alb
         else
             aub = min(-c1516*alb,-dx(j)+bub(j))
         end if
     else
         aub = b(j)

         if ( dx(j) == 0.0D+00 ) then
             alb = -c1516*aub
         else
             alb = max(-c1516*aub,-dx(j)+blb(j))
         end if

     end if
!
!  THIS NEXT CODE, ENDING WITH ***, POINTS THE BOX TOWARDS THE BEST
!  VALUE OF X WHEN NOT AT A NEW BEST.
!
     if ( k>=2) then
         if ( xb(j) > x(j)) then
             bub(j) = bub(j) * 0.25D+00

             if ( x(j)-blb(j) > xb(j)) then
                 blb(j) = blb(j) * 0.75D+00
             else
                 blb(j) = min(blb(j),0.75D+00 * (x(j)-xb(j)))
             end if

         else

             blb(j) = blb(j) * 0.25D+00

             if ( x(j)-bub(j) < xb(j)) then
                 bub(j) = bub(j) * 0.75D+00 
             else
                 bub(j) = max(bub(j),0.75D+00 * (x(j)-xb(j)))
             end if

         end if

     end if
!
!  RESTRICT THE STEP FURTHER IF USER GIVES BOUNDS.
!
     icase = ind(j)
     go to (540,550,560,570),icase
  540    aub = min(aub,x(j)-bl(j))
     go to 580
  550    alb = max(alb,x(j)-bu(j))
     go to 580
  560    aub = min(aub,x(j)-bl(j))
     alb = max(alb,x(j)-bu(j))
     go to 580
  570    continue
  580    blb(j) = alb
!
!  THIS NEXT LINE IS TO GUARANTEE THAT THE LOWER BOUND
!  IS .LE. THE UPPER BOUND.
!
     aub = max(aub,alb)
     bub(j) = aub
     indb(j) = 3

  end do
!
!  COMPUTE JACOBIAN*DX AND COMPARE NORM WITH CURRENT FUNCTION.
!
  if ( np > 1) then

      pj(1:lk) = 0.0D+00

      do j = 1,nvars

         if ( jactri) then
             jk = j
         else
             jk = mequa
         end if

         call daxpy(jk,dx(j),fjac(mcon+1,j),1,pj,1)

      end do

      t = dnrm2(lk,pj,1)
!
!  THIS TEST SAYS TO USE THE QUADRATIC MODEL IF
!  THE LAST STEP IS APPROXIMATELY IN THE NULL SPACE OF THE JACOBIAN.
!
      useq = useq .or. t  <  dfn * 0.75D+00

      if ( dfn > 0.0D+00 ) then
          dfn = t/dfn
      else
          dfn = 0.0D+00
      end if

  end if

  if ( iprint > 0) then
      call dvout(1,dfn,'('' ratio of j*dx norm to past df norm'')',4)
  end if
!
!  CHECK IF QUAD. MODEL IS BEING SUPPRESSED.
!
  useq = useq .and. .not. noquad
!
!  START THE PROCESS USING THE LINEAR MODEL.
!
  linmod = .true.
  mconst = mcon

  610 continue
!
!  COMPUTE THE REQUIRED DIMENSIONS OF THE MODEL PROBLEM.
!
  if ( linmod) then
      mk = 0
      me = min(mequa,nvars+1)
!
!  SET THE INITIAL VALUES FOR THE LINEAR MODEL PROBLEM.
!
      dx(1:nvars) = 0.0D+00
  else if ( useq) then
      mk = min(mequa,nvars+np+1)
      me = nvars + mk
      call dcopy(mk,fjac(mcon+1,nvars+1),1,dx(nvars+1),1)
  else
      go to 730
  end if

  nv = nvars + mk
!
!  NOTE THAT THE RESIDUALS ARE FREE VARIABLES.
!
  do i = nvars + 1,nv
     indb(i) = 4
  end do

  nit = 0
!
!  THE JACOBIAN, RIGHT SIDE, QUAD. TERMS ARE AN
!  UPPER TRAPeZOIDAL DATA ARRAY.  THIS WILL MAKE SOLVING
!  THE MODEL PROBLEM MORE EFFICIENT.
!
  630 continue
!
!  CALL A SIMPLIFIED VERSION OF THE ALGORITHM TO SOLVE
!  THE MODEL PROBLEM.  THE FUNCTION AND JACOBIAN
!  ARE COMPUTED FOR THIS SUBPROBLEM DURING THE REVERSE
!  COMMUNICATION REQUESTS.
!
  call dqedgn(me,nv,mconst,indb,blb,bub,dx,wj,ldwj,pv,igow, &
    iopt(ipls),ropt,iwa,wa)
!
!  CHECK FOR AN ERROR THAT WAS SEEN IN THE LOW-LEVEL NONLINEAR SOLVER.
!
  if ( igow > 7) then
    igo = igow
    iflag = 0
    go to 1020
  end if
!
!  CLEAR OUT THE WJ(*,*) ARRAY THAT HOLDS
!  THE JACOBIAN FOR THE INNER LOOP PROBLEM.
!
  wj(1,1) = 0.0D+00
  call dcopy(ldwj* (nv+1),wj,0,wj,1)
  if ( useq .and. .not. linmod ) wj(mconst+1,nvars+1) = one
!
!  PUT IN A UNIT MATRIX FOR THE PARTIALS
!  WITH RESPECT TO THE RESIDUALS.
!
  call dcopy(mk,wj(mconst+1,nvars+1),0,wj(mconst+1,nvars+1),ldwj+1)
!
!  THE FORM OF THE UPDATE BEING COMPUTED IS X(*)-DX(*).
!  THE VALUE OF DX IS ITSELF COMPUTED AS THE SOLUTION
!  TO A NONLINEAR PROBLEM WITH UPDATES OF THE FORM
!  DX(*)-D(DX)(*).
!
  do i = 1,mcon
     call dcopy(nvars,fjac(i,1),ldfjac,wj(i,1),ldwj)
     wj(i,nv+1) = fjac(i,nvars+1) - ddot(nvars,dx,1,wj(i,1),ldwj)
  end do
!
!  SEE IF USER HAS GIVEN GENERAL CONSTRAINTS.
!  USE THESE CONSTRAINTS TO PLACE EQUIVALENT CONSTRAINTS
!  ON THE CHANGES BEING COMPUTED.
!
  do j = 1,mcon
     blb(nv+j) = bl(j+nvars)
     bub(nv+j) = bu(j+nvars)
     indb(nv+j) = ind(j+nvars)
  end do
!
!  EVALUATE LINEAR MODEL
!

  assign 660 to igoelm
  go to 940

  660 continue
!
!  EVALUATE QUADRATIC MODEL
!
  assign 670 to igoeqm
  go to 850

  670 continue


  if ( igow > 1) then
      pv = dnrm2(me,wj(mconst+1,nv+1),1)
      if ( linmod) then
          pvl = pv
          dxl(1:nvars) = dx(1:nvars)
          linmod = .false.
          go to 610
!
!  IF THE PREDICTED NORM IS GREATER THAN THE CURRENT
!  RESIDUAL NORM, DROP THE QUADRATIC MODEL AND USE THE
!  LINEAR MODEL.
!
      else if ( (pv>=fc) .and. useq) then
          if ( iprint > 0) then
            write(*,*) 'Abandon quadratic model.'
          end if
          useq = .false.
      end if
      go to 730
  end if
!
!  FOR EITHER CASE TRANSFER THE JACOBIAN FOR THE MODEL
!  PROBLEM.  THE TRANSPOSE OF THIS MATRIX IS THE PARTIALS
!  WITH RESPECT TO THE RESIDUALS.
!
  do j = 1,nvars
     call dcopy(mk,wj(mconst+1,j),1,wj(mconst+mk+j,nvars+1),ldwj)
  end do
!
!  NOW UPDATE THE RESIDUALS FOR BOTH SETS OF MODEL EQUATIONS.
!  IN ROWS 1,...,MK THIS INVOLVES ADDING DX(NVARS+I) TO ROW I.
!  FOR ROWS MK+1,...,ME THIS REQUIRES ADDING MULTIPLES OF THE
!  COLS. OF THE TRANSPOSED JACOBIAN.
!
  do i = 1,mk
     t = dx(nvars+i)
     wj(mconst+i,nv+1) = wj(mconst+i,nv+1) + t
  end do
!
!  SYMMETRIZE THE SECOND DERIVATIVE MATRIX.  THIS
!  IS NOT REQUIRED WHEN THE MODEL IS LINEAR, BUT IT
!  DOES NOT HURT THEN.
!
  if ( useq .and. .not. linmod) then
      do j = 1,nvars
         do i = j,nvars
            wj(mconst+mk+i,j) = wj(mconst+mk+j,i)
         end do
      end do
  end if
!
!  COMPUTE RESIDUALS ON THE EQUATIONS K*R = 0.
!
  do j = nvars + 1,nv
     call daxpy(nvars,dx(j),wj(mconst+mk+1,j),1,wj(mconst+mk+1,nv+1),1)
  end do

  nit = nit + 1
  go to 630

  730 continue

!
!  COMPUTE THE ANGLES BETWEEN THE LINEAR AND QUADRATIC STEP.
!  TAKE THE ONE, IF THERE IS A CHOICE, CLOSEST TO THE GRADIENT.
!  IF THE QUADRATIC MOVE IS QUITE CLOSE TO THE GRADIENT, TAKE
!  THAT MOVE IN PREFERENCE TO THE LINEAR MOVE.
!
  cosl = dot_product ( gr(1:nvars), dxl(1:nvars) )

  t = dnrm2(nvars,dxl,1)
  if ( t > 0.0D+00 ) cosl = cosl/t
  cosq = -one
  cosm = -one
  tt = 0.0D+00


  if ( useq) then
      cosq = dot_product ( gr(1:nvars), dx(1:nvars) )
      tt = dnrm2(nvars,dx,1)
      if ( tt > 0.0D+00 ) cosq = cosq/tt
!
!  COMPUTE THE COSINE OF THE ANGLE BETWEEN THE QUAD. AND
!  LINEAR MOVES.
!
      if ( t > 0.0D+00 .and. tt > 0.0D+00 ) then
        cosm = dot_product ( dx(1:nvars), dxl(1:nvars) ) / t / tt
      end if

  end if

  if ( iprint > 0) then
    write(*,*)'cos of quad. move and grad., cos of lin. move ' &
      //'and grad., cosof each move'
    write(*,*)'flag for trying quad. move.'
    write(*,'(1p3d12.4,l6)') cosq,cosl,cosm,useq
  end if

  if ( iprint > 0)then
    write(*,*)'length of quad., then linear moves'
    write(*,'(1p2d12.4)') tt,t
  end if
!
!  CHOOSE MOVE PARTIALLY BASED ON ANGLE MOVES MAKE WITH EACH OTHER.
!
  useq = useq .and. (cosm > 0.0D+00 .or. cosl < cosm) .and. cosq > 0.0D+00
  useql = useq

  if ( .not. useq) then
      pv = pvl
      call dcopy(nvars,dxl,1,dx,1)
      ntterm = 0
  else
      ntterm = np - 1
  end if
!
!  TEST FOR NOISE IN MODEL PROBLEM SOLN.
!
  term = (pv>=fc) .and. .not. retrea .and. .not. useq
  term = term .and. mcon  ==  0
  term = .false.
  if ( term) then
      if ( iprint > 0) then
          write(*,9021) pv,fc
      end if
!
!  VALUE MEANS MODEL RES. .GE. NONLINEAR FUNCTION VALUE.
!
      igo = 5
      iflag = 0
      go to 840
  end if

  if ( pv > pb .and. pv + pd /= 0.0D+00 ) then
    rc = 4.0D+00 * (pv-pb)/ (pv+pd)
  else
    rc = 0.0D+00
  end if
!
!  IF USING A QUADRATIC MODEL AND RETREATING SEEMS TO BE
!  NECESSARY, SEE IF RETREATING WOULD BE NEEDED WITH A
!  LINEAR MODEL.  ONLY THEN RETREAT.
!
  if ( rc<=one .or. .not. useq) go to 750
!
!  EVALUATE LINEAR MODEL
!
  nv = nvars
  assign 740 to igoelm
  go to 940

  740 continue

  pvl = dnrm2(min(mequa,nvars+1),wj(mconst+1,nv+1),1)

  if ( pvl > pb ) rc = 4.0D+00 * (pvl-pb)/ (pvl+pd)

 750 continue

  rg = max(rg,rc)

  if ( iprint > 0) then
      write(*,9011) iters,fc,pv,rc,ajn,k,kl,fb,alpha,bboost,nit,useq,ntterm
      write(*,9001) '  x=', (x(j),j=1,nvars)
      write(*,9001) ' dx=', (dx(j),j=1,nvars)
      write(*,9001) '  b=', (b(j),j=1,nvars)
      write(*,9001) ' lb=', (blb(j),j=1,nall)
      write(*,9001) ' ub=', (bub(j),j=1,nall)
      write(*,'('' ///end of iteration.///'')')
  end if

  retrea = rc  >  1

  if ( .not. retrea) then
      chg = one
      t2 = 0.0D+00
      do 810 j = 1,nvars
         bold = b(j)
         t = dx(j)/bold
         alb = 2.0D+00
         aub = 2.0D+00
!
!  IF USER GIVES BOUNDS, AND THESE BOUNDS ARE HIT,
!  DO NOT DETRACT FROM DECLARING A FULL NEWTON STEP.
!
         icase = ind(j)

         if ( ind(j) == 1 ) then
           alb = (x(j)-bl(j))/bold
           aub = -semibg
         else if ( ind(j) == 2 ) then
           aub = (x(j)-bu(j))/bold
           alb = -semibg
         else if ( ind(j) == 3 ) then
           alb = (x(j)-bl(j))/bold
           aub = (x(j)-bu(j))/bold
         else if ( ind(j) == 4 ) then
           alb = -semibg
           aub = -semibg
         end if

         if ( t == one) then
             t2 = one
             b(j) = bold + bold
             chg = chg*chgfac
         else
             if ( abs(t) < 0.25D+00 .and. dx(j)/= 0.0D+00 ) then
                 b(j) = sign ( 0.25D+00 * bold,dx(j)) + 3.0D+00 * dx(j)
             else
                 b(j) = sign(bold,dx(j))
             end if
         end if
!
!  THIS TEST AVOIDS THE USER BOUNDS IN DECLARING A NEWTON STEP.
!
         if ( abs(alb-t)>=0.01D+00*abs(t) .and. &
              abs(aub-t) >= 0.01D+00*abs(t)) then
             if ( t > 0.0D+00 ) then
                 t2 = max(t2,t)
             else
                 t2 = max(t2,-t/c1516)
             end if
         end if

  810     continue

      fulnwt = t2  <  0.99D+00
      dxnrm = abs(dx(idamax(nvars,dx,1)))
!
!  TEST FOR SMALL ABSOLUTE CHANGE IN X VALUES.
!
      term = dxnrm  <  told .and. fulnwt
      if ( term) then
          igo = 6
!
!  VALUE MEANS CHANGE (IN PARAMETERS) WAS SMALL AND A
!  FULL STEP (NOT HITTING TRUST CONSTRAINTS) TAKEN.
!
          go to 820
      end if

      term = dxnrm  <  dnrm2(nvars,x,1)*tolx .and. fulnwt
      term = term .and. (iters > 1)
      if ( term) then
          igo = 7
!
!  VALUE MEANS RELATIVE CHANGE IN PARAMETERS WAS SMALL AND A
!  FULL STEP (NOT HITTING CONSTRAINTS WITH AT LEAST 2 ITERATIONS)
!  WAS TAKEN.
!
          go to 820
      end if

      go to 830

  820     continue

      go to 30

  end if

  830 continue
  fl = fc

  go to 30

  840 continue
  go to 1020
!
!  EVALUATE QUADRATIC MODEL
!
  850 continue
!
!  IF THE MODEL IS GENUINELY QUADRATIC, ADD IN THE EXTRA
!  TERMS AND COMPUTE THE SECOND DERIVATIVE INFORMATION.
!
  if ( useq .and. .not. linmod ) then

!
!  COMPUTE THE DOT PRODUCT OF CURRENT PROPOSED STEP AND
!  PAST DIRECTIONS REPRESENTED IN THE MODEL.
!
      do l = 1,np - 1
         t = 0.0D+00
         do j = 1,nvars
            t = t + dx(j)* (xp(j,l+1)-xp(j,1))
         end do
         pj(l) = t
      end do
!
!  STORAGE LAYOUT, WITH K = J**T, OF WJ(*,*).
!    [J    :  I    : F+R ]
!    [H    :  K    : K*R ]
!  ADD IN THE QUADRATIC TERMS FOR THE FUNCTION.
!
      do l = 1,np - 1
         jk = min(nvars+l+1,mequa)
         call daxpy(jk,0.5D+00*pj(l)**2,qc(1,l+1),1,wj(mconst+1,nv+1),1)
      end do
!
!  ADD THE LINEAR TERMS TO THE INNER LOOP JACOBIAN.
!
      do l = 1,np - 1
         jk = min(nvars+l+1,mequa)
         do j = 1,nvars
            call daxpy(jk,pj(l)* (xp(j,l+1)-xp(j,1)),qc(1,l+1),1, &
                       wj(mconst+1,j),1)
         end do
      end do
!
!  COMPUTE THE UPPER TRIANGULAR PART OF THE SECOND DERIVATIVE TERMS.
!
      do i = 1,nvars
         do j = i,nvars
            do l = 1,np - 1
               jk = min(nvars+l+1,mequa)
               wj(mconst+mk+i,j) = wj(mconst+mk+i,j) + &
                 (xp(j,l+1)-xp(j,1))*(xp(i,l+1)-xp(i,1))* &
                 ddot (jk,dx(nvars+1),1,qc(1,l+1),1)
            end do
          end do
        end do

  end if

  go to igoeqm
!
!  EVALUATE LINEAR MODEL
!
  940 continue
!
!  TRANSFER THE JACOBIAN THAT WOULD RESULT FROM
!  USING JUST A LINEAR MODEL.
!
  do j = 1,nvars
    if ( jactri) then
      jk = j
    else
      jk = mequa
    end if

    call dcopy(jk,fjac(mcon+1,j),1,wj(mconst+1,j),1)
  end do
!
!  TRANSFER THE PRESENT VALUES OF THE FUNCTION.
!
  call dcopy(min(mequa,nvars+1),fjac(mcon+1,nvars+1),1,wj(mconst+1,nv+1),1)
!
!  CHANGE SIGN FOR THE MODEL PROBLEM.
!
  do i = 1,min(mequa,nvars+1)
    wj(mconst+i,nv+1) = -wj(mconst+i,nv+1)
  end do
!
!  COMPUTE THE LINEAR TERM OF THE MODEL.
!
  do j = 1,nvars
     if ( jactri) then
         jk = j
     else
         jk = mequa
     end if

     call daxpy(jk,dx(j),wj(mconst+1,j),1,wj(mconst+1,nv+1),1)
  end do

  go to igoelm
!
!  TEST FOR CONVERGENCE
!
  980 continue

  term = iters >= itmax
  if ( term) then
      igo = 8
!
!  VALUE MEANS THAT MAX. NUMBER OF ALLOWED ITERATIONS TAKEN.
!
      go to 1000
  end if
!
!  TEST FOR SMALL FUNCTION NORM.
!
  term = fc <= tolf
!
!  IF HAVE CONSTRAINTS MUST ALLOW AT LEAST ONE MOVE.
!
  term = term .and. (mcon == 0 .or. iters > 1)
  if ( term) then
      igo = 2
!
!  VALUE MEANS FUNCTION NORM WAS SMALL.
!
      go to 1000
  end if
!
!  TEST FOR NO CHANGE
!
  go to 1010

  990 continue

  term = term .and. .not. retrea
  if ( term) then
      igo = 3
!
!  VALUE MEANS THE FUNCTION IS PROBABLY REACHING A LOCAL MINIMUM
!  BUT MOVES ARE STILL HITTING TRUST REGION CONSTRAINTS.
!
      if ( fulnwt) igo = 4
!
!  VALUE MEANS THAT FUNCTION IS REACHING A LOCAL MINIMUM
!  AND MOVES ARE NOT HITTING THE TRUST REGION CONSTRAINTS.
!
      if ( igo == 3) term = term .and. .not. mustcn
      go to 1000
  end if

 1000 continue
  go to igotfc
!
!  TEST FOR NO CHANGE
!
 1010 continue

  term = (abs(fb-pv)<=tolsnr*fb) .and. (abs(fc-pv).le.fb*tolp)
  term = term .and. (abs(fc-fl)<=fb*tolsnr)
  term = term .and. (abs(pvl-pv)<=fb*tolsnr)

  go to 990

 1020 continue

  return
!
!  INITIALIZE OTHER VALUES
!
 1030 continue
!
!  THE NUMBER OF PAST DIFFERENCES USED IN THE QUADRATIC MODEL.
!
  np = 0
!
!  IF NO MORE EQUATIONS THAN VARIABLES, NO NEED TO
!  PRETRIANGULARIZE THE JACOBIAN MATRIX.
!
  jactri = ( nvars < mequa )
!
!  MAKE SURE THAT VARIABLES SATISFY CONSTRAINTS.
!  GENERALLY THIS MAY TAKE A CALL TO DBOCLS().
!  AS LONG AS THE FUNCTIONS ARE DEFINED AT POINTS
!  THAT DO NOT SATISFY THE CONSTRAINTS, THE FIRST
!  ALGORITHM STEP WILL BRING IT ONTO THE CONSTRAINTS.
!
  do j = 1, nvars

    if ( ind(j) == 1 ) then
      x(j) = max ( x(j), bl(j) )
    else if ( ind(j) == 2 ) then
      x(j) = min ( x(j), bu(j) )
    else if ( ind(j) == 3 ) then
      x(j) = max(x(j),bl(j))
      x(j) = min(x(j),bu(j))
    else

    end if

  end do

  iters = 0
  nall = mcon + nvars
  chgfac = 2.0D+00** ( -one / real ( nvars, kind = 8 ) )
  c1516 = 15.0D+00 / 16.0D+00
  semibg = 1.0D+10

  go to 20
!
!  PROCESS OPTION ARRAY
!
 1100 continue
  iprint = 0
!
!  D1MACH(4)=RELPR=MACHINE REL. PREC.
!
  t = epsilon ( t )
  tolf = sqrt(t)
  toluse = tolf
  tolx = tolf
  told = tolf
  tolsnr = 1.0D-05
  tolp = 1.0D-05
  cond = 30.0D+00
  itmax = 75
  level = 1
  ipls = 0
  passb = .false.
  noquad = .false.
  revers = .false.
  mustcn = .false.
  lp = 1
  lpdiff = 0

 1110 continue

  lp = lp + lpdiff
  lpdiff = 2
  kp = iopt(lp)
  newopt = kp  >  0
  jp = abs(kp)
!
!  SEE IF THIS IS THE LAST OPTION..
!
!  THE POINTER TO THE START OF OPTIONS FOR THE INNER LOOP
!  SOLVER MUST SATISFY THE REQUIREMENTS FOR THAT OPTION ARRAY.
!
  if ( jp == 99) then
      if ( newopt) then
          if ( ipls == 0) ipls = lp
          go to 1120
      else
          lpdiff = 1
          go to 1110
      end if
  end if
!
!  CHANGE PRINT OPTION.
!
  if ( jp == 1) then
      if ( newopt) iprint = iopt(lp+1)
      go to 1110
  end if
!
!  SEE IF MAX. NUMBER OF ITERATIONS CHANGING.
!
  if ( jp == 2) then
      if ( newopt) itmax = iopt(lp+1)
      go to 1110
  end if
!
!  SEE IF BOUNDS FOR THE TRUST REGION ARE BEING PASSED.
!
  if ( jp == 3) then
      if ( newopt) then
          call dcopy(nvars,ropt(iopt(lp+1)),1,bb,1)
          passb = .true.
      end if
      go to 1110
  end if
!
!  CHANGE TOLERANCE ON THE LENGTH OF THE RESIDUALS.
!
  if ( jp == 4) then
      if ( newopt) tolf = ropt(iopt(lp+1))
      go to 1110
  end if
!
!  CHANGE TOLERANCE ON THE NORM OF THE RELATIVE
!  CHANGE TO THE PARAMETERS.
!
  if ( jp == 5) then
      if ( newopt) tolx = ropt(iopt(lp+1))
      go to 1110
  end if
!
!  CHANGE TOLERANCE ON ABSOLUTE CHANGE TO THE PARAMETERS.
!
  if ( jp == 6) then
      if ( newopt) told = ropt(iopt(lp+1))
      go to 1110
  end if
!
!  CHANGE TOLERANCE FOR RELATIVE AGREEMENT BETWEEN
!  BEST FUNCTION NORM, LAST FUNCTION NORM AND THE
!  CURRENT FUNCTION NORM.
!
  if ( jp == 7) then
      if ( newopt) tolsnr = ropt(iopt(lp+1))
      go to 1110
  end if
!
!  CHANGE TOLERANCE FOR AGREEMENT BETWEEN PREDICTED
!  VALUE OF RESIDUAL NORM AND THE PREVIOUS VALUE OF
!  THE RESIDUAL NORM.
!
  if ( jp == 8) then
      if ( newopt) tolp = ropt(iopt(lp+1))
      go to 1110
  end if
!
!  CHANGE TOLERANCE SUCH THAT RELATIVE CHANGES IN THE
!  VALUES OF THE PARAMETERS IMPLY THAT THE PREVIOUS
!  VALUE OF THE FUNCTION WILL NOT BE USED IN THE
!  QUADRATIC MODEL.
!
  if ( jp == 9) then
    if ( newopt) toluse = ropt(iopt(lp+1))
    go to 1110
  end if
!
!  CHANGE THE LARGEST CONDITION NUMBER TO ALLOW WHEN
!  SOLVING FOR THE QUADRATIC COEFFICIENTS OF THE MODEL.
!
  if ( jp == 10) then
      if ( newopt) cond = ropt(iopt(lp+1))
      go to 1110
  end if
!
!  CHANGE THE PRINT LEVEL IN THE ERROR PROCESSOR.
!
  if ( jp == 11) then
      if ( newopt) level = iopt(lp+1)
      go to 1110
  end if
!
!  PASS AN OPTION ARRAY TO THE CONSTRAINED LINEAR SOLVER.
!  THIS OPTION IS A POINTER TO THE START OF THE OPTION
!  ARRAY FOR THE SUBPROGRAM.
!
  if ( jp == 12) then
      if ( newopt) ipls = iopt(lp+1)
      go to 1110
  end if
!
!  MOVE THE PROCESSING POINTER BY THE VALUE IN THE
!  NEXT ENTRY OF THE OPTION ARRAY.  THIS DEVICE IS
!  INCLUDED SO THAT PASSING OPTIONS TO LOWER LEVEL
!  SUBROUTINES IS EASY TO DO.
!
  if ( jp == 13) then
      if ( newopt) lpdiff = iopt(lp+1)
      go to 1110
  end if
!
!  OPTION TO SUPPRESS USING THE QUADRATIC MODEL, EVER.
!
  if ( jp == 14) then
      if ( newopt) noquad = iopt(lp+1)  ==  1
      go to 1110
  end if
!
!  MORE STORAGE WAS GIVEN FOR THE QUADRATIC MODEL ARRAYS.
!  THIS OPTION WAS PROCESSED BY THE INTERFACE UNIT.
!
  if ( jp == 15) go to 1110
!
!  USE FORWARD COMMUNICATION TO GET THE DERIVATIVES
!  AND FUNCTION VALUES.
!
  if ( jp == 16) then
      if ( newopt) revers = iopt(lp+1)  ==  1
      go to 1110
  end if
!
!  FORCE A FULL NEWTON STEP WHEN NEAR THE MINIMUM.
!  DO NOT ALLOW CONVERGENCE CLAIMS WHEN HITTING BOUNDS.
!
  if ( jp == 17) then
      if ( newopt) mustcn = iopt(lp+1)  ==  1
      go to 1110
  end if
!
!  SAW AN OPTION (OR GARBAGE) THAT IS NOT ON THE LIST.
!
  xmess ='dqedmn. invalid option processed. i1=iopt(*) entry. i2=iopt(i1).'
  nerr = 07
  igo = 15
  call xerrwv(xmess,nerr,level,2,lp,iopt(lp),0,rdum,rdum)
  iflag = 0
  go to 1020
 1120 continue
  go to 10

 9001 format (a4,1p10d12.4/ (4x,10d12.4))
 9011 format ('0iter.=',i3,' fc=',1pd10.4,' pv=',1pd10.4,03x,' rc=', &
    1pd10.4,' j**t*f=',1pd10.4,/,' k=',i4,' kl=',i4,/10x,' fb=', &
    1pd10.4,' al=',1pd10.4,' bb=',1pd12.4/' inner iterations =',i5, &
    ' use quad. model?=',l5,' num. of terms =',i5)
 9021 format (' model residual>=current f. quitting.',1p2d12.5)
end
subroutine drot ( n, x, incx, y, incy, c, s )

!*******************************************************************************
!
!! DROT applies a plane rotation.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to be rotated.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to be rotated.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
!    Input, real ( kind = 8 ) C, S, parameters (presumably the cosine and sine 
!    of some angle) that define a plane rotation.
!
  implicit none

  real ( kind = 8 ) c
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real ( kind = 8 ) s
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    do i = 1, n
      stemp = c * x(i) + s * y(i)
      y(i) = c * y(i) - s * x(i)
      x(i) = stemp
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = c * x(ix) + s * y(iy)
      y(iy) = c * y(iy) - s * x(ix)
      x(ix) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine drotg ( sa, sb, c, s )

!*******************************************************************************
!
!! DROTG constructs a Givens plane rotation.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) SA, SB, ...
!
!    Output, real ( kind = 8 ) C, S, ...
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) r
  real ( kind = 8 ) roe
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) scale
  real ( kind = 8 ) z

  if ( abs ( sa ) > abs ( sb ) ) then
    roe = sa
  else
    roe = sb
  end if

  scale = abs ( sa ) + abs ( sb )

  if ( scale == 0.0D+00 ) then
    c = 1.0D+00
    s = 0.0D+00
    r = 0.0D+00
  else
    r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
    r = sign ( 1.0D+00, roe ) * r
    c = sa / r
    s = sb / r
  end if

  if ( abs ( c ) > 0.0D+00 .and. abs ( c ) <= s ) then
    z = 1.0D+00 / c
  else
    z = s
  end if

  sa = r
  sb = z

  return
end
subroutine dscal ( n, sa, x, incx )

!*******************************************************************************
!
!! DSCAL scales a vector by a constant.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries of X.
!
  implicit none

  integer i
  integer incx
  integer ix
  integer m
  integer n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine dswap ( n, x, incx, y, incy )

!*******************************************************************************
!
!! DSWAP interchanges two vectors.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to swap.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to swap.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
  implicit none

  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      stemp = x(i)
      x(i) = y(i)
      y(i) = stemp
    end do

    do i = m+1, n, 3

      stemp = x(i)
      x(i) = y(i)
      y(i) = stemp

      stemp = x(i + 1)
      x(i + 1) = y(i + 1)
      y(i + 1) = stemp

      stemp = x(i + 2)
      x(i + 2) = y(i + 2)
      y(i + 2) = stemp

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = x(ix)
      x(ix) = y(iy)
      y(iy) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine dvout ( n, dx, ifmt, idigit )

!***********************************************************************
!
!! DVOUT prints double precision vectors.
!
!  Example:
!
!    PRINT AN ARRAY CALLED (COSTS OF PURCHASES) OF LENGTH 100 SHOWING
!    6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
!    SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
!
!      double precision COSTS(100)
!      N = 100
!      IDIGIT = -6
!      CALL DVOUT(N,COSTS,'(''1COSTS OF PURCHASES'')',IDIGIT)
!
!  Author:
!
!    JOHN A. WISNIEWSKI and RICHARD J. HANSON
!    SANDIA LABS ALBUQUERQUE.
!
!  INPUT..
!
!  N,DX(*) PRINT THE double precision ARRAY DX(I),I=1,...,N, ON
!          OUTPUT UNIT *. THE HEADING IN THE FORTRAN FORMAT
!          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST
!          STEP. THE COMPONENTS DX(I) ARE INDEXED, ON OUTPUT,
!          IN A PLEASANT FORMAT.
!  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT
!          UNIT * WITH THE VARIABLE FORMAT FORTRAN STATEMENT
!                write(*,IFMT)
!  IDIGIT  PRINT AT LEAST IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
!          THE SUBPROGRAM WILL CHOOSE THAT integer 6,14,20 OR 28
!          WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF
!          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
!          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY DX(*). (THIS
!          CAN BE USED ON MOST TIME-SHARING TERMINALS). IF
!          IDIGIT.GE.0, 133 PRINTING COLUMNS ARE UTILIZED. (THIS CAN
!          BE USED ON MOST LINE PRINTERS).
!
  implicit none

  real ( kind = 8 ) dx(*)
  integer i
  integer idigit
  character ( len = * ) ifmt
  integer k1
  integer k2
  integer n
  integer ndigit

  write ( *, ifmt )

  if ( n <= 0 ) then
    return
  end if

  ndigit = idigit

  if ( idigit == 0) ndigit = 6

  if ( idigit < 0 ) then

    ndigit = -idigit

    if ( ndigit <= 6 ) then

      do k1=1,n,4
        k2 = min(n,k1+3)
        write(*,1000) k1,k2,(dx(i),i=k1,k2)
      end do

    else if ( ndigit <= 14 ) then

      do k1=1,n,2
        k2 = min0(n,k1+1)
        write(*,1001) k1,k2,(dx(i),i=k1,k2)
      end do

    else if ( ndigit <= 20 ) then

      do k1=1,n,2
        k2=min0(n,k1+1)
        write(*,1002) k1,k2,(dx(i),i=k1,k2)
      end do

    else

      do k1=1,n
        k2 = k1
        write(*,1003) k1,k2,(dx(i),i=k1,k2)
      end do

    end if

  else

    if ( ndigit <= 6 ) then

      do k1=1,n,8
        k2 = min0(n,k1+7)
        write(*,1000) k1,k2,(dx(i),i=k1,k2)
      end do

    else if ( ndigit <= 14 ) then

      do k1=1,n,5
        k2 = min(n,k1+4)
        write(*,1001) k1,k2,(dx(i),i=k1,k2)
      end do

    else if ( ndigit <= 20 ) then

      do k1=1,n,4
        k2 = min(n,k1+3)
        write(*,1002) k1,k2,(dx(i),i=k1,k2)
      end do

    else

      do k1=1,n,3
        k2 = min(n,k1+2)
        write(*,1003) k1,k2,(dx(i),i=k1,k2)
      end do

    end if

  end if

  return
 1000 format(1x,i4,' - ',i4,1x,1p8e14.5)
 1001 format(1x,i4,' - ',i4,1x,1p5e22.13)
 1002 format(1x,i4,' - ',i4,1x,1p4e28.19)
 1003 format(1x,i4,' - ',i4,1x,1p3e36.27)
end
subroutine timestamp ( )

!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine xerrwv ( xmess, nerr, level, ni, i1, i2, nr, r1, r2 )

!***********************************************************************
!
!! XERRWV is an error output message routine.
!
  implicit none

  integer i1
  integer i2
  integer level
  integer nerr
  integer ni
  integer nr
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  character ( len = * ) xmess
!
!     if ( level < 1) return
!
  write(*,'(1x,a)') trim ( xmess )
  write(*,'('' error number = '',i5,'', message level = '',i5)')nerr,level

  if ( ni == 1 .or. ni == 2 )then
    write(*,'('' i1 = '',i8)') i1
  end if

  if ( ni == 2)then
    write(*,'('' i2 = '',i8)') i2
  end if

  if ( nr == 1 .or. nr == 2)then
    write(*,'('' r1 = '',1pe15.7)') r1
  end if

  if ( nr == 2)then
    write(*,'('' r2 = '',1pe15.7)') r2
  end if

  return
end
