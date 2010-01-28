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

subroutine caseDirect(x, fj, ldfj, igo, iopt, ropt)
	! Objective function and its Jacobian for the case of directly-fitted
	! oscillators and rotors.

	use params
	implicit none

	! Input/output parameters
	integer, intent(in) :: ldfj
	real(8), dimension(1:nvars), intent(out) :: x
	real(8), dimension(1:ldfj,1:nvars+1), intent(out) :: fj
	integer, intent(inout) :: igo
	integer, dimension(1:24), intent(in) :: iopt
	real(8), dimension(1:1), intent(in) :: ropt

	real(8), dimension(1:nT) :: diff
	real(8), dimension(1:Nvib+Nrot, 1:nT) :: Cv
	real(8), dimension(1:Nvib+2*Nrot, 1:nT) :: dCv

	integer i, j

	! Calculate the heat capacities and derivatives at each temperature
	do i = 1, Nvib
		call harmonicOscillator_heatCapacity(Tlist, nT, x(i), Cv(i,:))
		call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(i), dCv(i,:))
	end do
	do i = 1, Nrot
		call hinderedRotor_heatCapacity(Tlist, nT, x(Nvib+2*i-1), x(Nvib+2*i), Cv(Nvib+i,:))
		call hinderedRotor_d_heatCapacity_d_freq(Tlist, nT, x(Nvib+2*i-1), x(Nvib+2*i), dCv(Nvib+2*i-1,:))
		call hinderedRotor_d_heatCapacity_d_barr(Tlist, nT, x(Nvib+2*i-1), x(Nvib+2*i), dCv(Nvib+2*i,:))
	end do

	! There are no constraints other than the bounds, so we don't need to
	! calculate the values of the constraints or their Jacobian

	! There is a minimization objective for each temperature
	! First calculate the least-squares objective at each temperature
	! This goes in the last column of the fj matrix
	diff = -Cvlist
	do i = 1, Nvib+Nrot
		diff = diff + Cv(i,:)
	end do
	fj(mcon+1:mcon+mequa,nvars+1) = diff * diff

	! We also need the Jacobian of the objective functions if igo is nonzero
	if (igo /= 0) then
		do j = 1, mequa
			do i = 1, Nvib+2*Nrot
				fj(mcon+i,j) = 2.0 * diff(j) * dCv(i,j)
			end do
		end do
	end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine caseNvibNrot(x, fj, ldfj, igo, iopt, ropt)
	! Objective function and its Jacobian for the case of many oscillators and
	! many rotors.
	
	use params
	implicit none

	! Input/output parameters
	integer, intent(in) :: ldfj
	real(8), dimension(1:nvars), intent(out) :: x
	real(8), dimension(1:ldfj,1:nvars+1), intent(out) :: fj
	integer, intent(inout) :: igo
	integer, dimension(1:24), intent(in) :: iopt
	real(8), dimension(1:1), intent(in) :: ropt

	real(8), dimension(1:nT) :: diff
	real(8), dimension(1:nT) :: Cv1, Cv2, Cv3, Cv4
	real(8), dimension(1:nT) :: dCv1, dCv2, dCv3, dCv4

	integer i

	! Calculate the heat capacities and derivatives at each temperature
	call harmonicOscillator_heatCapacity(Tlist, nT, x(2), Cv1)
	call harmonicOscillator_heatCapacity(Tlist, nT, x(3), Cv2)
	call hinderedRotor_heatCapacity(Tlist, nT, nu_low, x(5), Cv3)
	call hinderedRotor_heatCapacity(Tlist, nT, nu_high, x(6), Cv4)
	call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(2), dCv1)
	call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(3), dCv2)
	call hinderedRotor_d_heatCapacity_d_barr(Tlist, nT, nu_low, x(5), dCv3)
	call hinderedRotor_d_heatCapacity_d_barr(Tlist, nT, nu_high, x(6), dCv4)

	! There are no constraints other than the bounds, so we don't need to
	! calculate the values of the constraints or their Jacobian

	! There is a minimization objective for each temperature
	! First calculate the least-squares objective at each temperature
	! This goes in the last column of the fj matrix
	diff = x(1) * Cv1 + (nvib - x(1)) * Cv2 + x(4) * Cv3 + (nrot - x(4)) * Cv4 - Cvlist
	fj(mcon+1:mcon+mequa,nvars+1) = diff * diff
	
	! We also need the Jacobian of the objective functions if igo is nonzero
	if (igo /= 0) then
		do i = 1, mequa
			! The first harmonic oscillator degeneracy
			fj(mcon+i,1) = 2.0 * diff(i) * (Cv1(i) - Cv2(i))
			! The first harmonic oscillator pseudo-frequency
			fj(mcon+i,2) = 2.0 * diff(i) * (x(1) * dCv1(i))
			! The second harmonic oscillator pseudo-frequency
			fj(mcon+i,3) = 2.0 * diff(i) * ((real(nvib) - x(1)) * dCv2(i))
			! The first hindered rotor degeneracy
			fj(mcon+i,4) = 2.0 * diff(i) * (Cv3(i) - Cv4(i))
			! The first hindered rotor pseudo-barrier
			fj(mcon+i,5) = 2.0 * diff(i) * (x(4) * dCv3(i))
			! The second hindered rotor pseudo-barrier
			fj(mcon+i,5) = 2.0 * diff(i) * ((real(nrot) - x(4)) * dCv4(i))
		end do
	end if

end subroutine
