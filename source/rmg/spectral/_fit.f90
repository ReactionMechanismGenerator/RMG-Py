!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	RMG - Reaction Mechanism Generator
!
!	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
!	RMG Team (rmg_dev@mit.edu)
!
! 	Permission is hereby granted, free of charge, to any person obtaining a
! 	copy of this software and associated documentation files (the 'Software'),
! 	to deal in the Software without restriction, including without limitation
! 	the rights to use, copy, modify, merge, publish, distribute, sublicense,
! 	and/or sell copies of the Software, and to permit persons to whom the
! 	Software is furnished to do so, subject to the following conditions:
!
! 	The above copyright notice and this permission notice shall be included in
! 	all copies or substantial portions of the Software.
!
! 	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! 	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! 	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! 	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! 	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! 	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! 	DEALINGS IN THE SOFTWARE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Items placed in this module will be exposed to Python.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setParams(p_nT, p_Tlist, p_Cvlist, p_mcon, p_mequa, p_nvars, &
	p_nvib, p_nrot)
	! Set the module parameters based on the values in the parameter list.

	integer, intent(in) :: p_nT								! The number of temperatures at which Cv data is given
	real(8), dimension(1:p_nT), intent(in) :: p_Tlist		! The list of temperatures in K at which Cv data is given
	real(8), dimension(1:p_nT), intent(in) :: p_Cvlist		! The list of heat capacity data (dimensionless) at each temperature
	integer, intent(in) :: p_mcon
	integer, intent(in) :: p_mequa
	integer, intent(in) :: p_nvars							! The number of variables used in the fitting
	integer, intent(in) :: p_nvib							! The number of harmonic oscillators to fit
	integer, intent(in) :: p_nrot							! The number of hindered rotors to fit
	
	call set_params(p_nT, p_Tlist, p_Cvlist, p_mcon, p_mequa, p_nvars, &
		p_nvib, p_nrot)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cleanup()
	! Clean up arrays allocated when setting up the module parameters.

	call clean_up()
	
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fitModes(x0, nx, bl, bu, ind, maxiter, xout, igo)
	! Execute a constrained nonlinear optimization to fit molecular degrees of
	! freedom to heat capacity data.

	implicit none

	integer, intent(in) :: nx
	real(8), dimension(1:nx), intent(in) :: x0			! The initial guess to the solver for each variable
	real(8), dimension(1:nx), intent(in) :: bl			! The lower bounds for each variable
	real(8), dimension(1:nx), intent(in) :: bu			! The lower bounds for each variable
	integer, dimension(1:nx), intent(in) :: ind			! The type of bounds for each variable
	integer, intent(in) :: maxiter						! The maximum number of iterations to try
	real(8), dimension(1:nx), intent(out) :: xout		! The returned solution vector
	integer, intent(out) :: igo							! The returned solution status

	call fit_modes(x0, nx, bl, bu, ind, maxiter, xout, igo)

end subroutine

