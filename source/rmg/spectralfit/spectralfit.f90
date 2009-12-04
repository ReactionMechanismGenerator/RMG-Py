!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   RMG - Reaction Mechanism Generator
!
!   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
!   RMG Team (rmg_dev@mit.edu)
! 
!   Permission is hereby granted, free of charge, to any person obtaining a
!   copy of this software and associated documentation files (the 'Software'),
!   to deal in the Software without restriction, including without limitation
!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!   and/or sell copies of the Software, and to permit persons to whom the
!   Software is furnished to do so, subject to the following conditions:
! 
!   The above copyright notice and this permission notice shall be included in
!   all copies or substantial portions of the Software.
! 
!   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!   DEALINGS IN THE SOFTWARE.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This module should be compiled with f2py to provide an interface to the
! frankie Fortran code from Python. It depends on calc_freq_code.f90 and
! dqed.f90.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fitSpectralData(Cv, Tlist, Ntemp, Nvib, Nhind, vib, hind)
    ! Estimate the spectroscopic degrees of freedom by fitting parameters to
    ! heat capacity data.
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `Cv`       in     A list of heat capacities Cv/R at various temperatures
    !                   for the unknown degrees of freedom (i.e. the known 
    !                   degrees of freedom should have already been removed)
    ! `Tlist`    in     The temperatures corresponding to the heat capacities
    ! `Ntemp`    in     The number of temperatures and heat capacities provided
    ! `Nvib`     in     The number of 1D quantum harmonic oscillators to fit
    ! `Nhind`    in     The number of 1D Pitzer hindered rotors to fit
    ! `vib`      out    A vector of fitted 1D quantum harmonic oscillator 
    !                   frequencies in cm^-1
    ! `hind`     out    A matrix of fitted 1D Pitzer hindered rotor frequency-
    !                   barrier pairs, both in cm^-1
    ! ========== ====== ========================================================

    integer, intent(in) :: Ntemp
    integer, intent(in) :: Nvib
    integer, intent(in) :: Nhind
    real(8), dimension(1:Ntemp), intent(in) :: Cv
    real(8), dimension(1:Ntemp), intent(in) :: Tlist
    real(8), dimension(1:Nvib), intent(out) :: vib
    real(8), dimension(1:Nhind,1:2), intent(out) :: hind
    
    if (Ntemp /= 7) then
        write (*,'(a)'), 'The number of temperatures must be exactly 7.'
        stop
    end if

    ! Fit the harmonic oscillator and hindered rotor modes to the heat capacity
    ! The function fitSpectralDataToHeatCapacity is in calc_freq_code.f90
    call fitSpectralDataToHeatCapacity(Cv, Tlist, Ntemp, Nvib, Nhind, vib, hind)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fitSpectralDataNoRotors(Cv, Tlist, Ntemp, Nvib, vib)
    ! Estimate the spectroscopic degrees of freedom by fitting parameters to
    ! heat capacity data.
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `Cv`       in     A list of heat capacities Cv/R at various temperatures
    !                   for the unknown degrees of freedom (i.e. the known
    !                   degrees of freedom should have already been removed)
    ! `Tlist`    in     The temperatures corresponding to the heat capacities
    ! `Ntemp`    in     The number of temperatures and heat capacities provided
    ! `Nvib`     in     The number of 1D quantum harmonic oscillators to fit
    ! `vib`      out    A vector of fitted 1D quantum harmonic oscillator
    !                   frequencies in cm^-1
    ! ========== ====== ========================================================

    integer, intent(in) :: Ntemp
    integer, intent(in) :: Nvib
    real(8), dimension(1:Ntemp), intent(in) :: Cv
    real(8), dimension(1:Ntemp), intent(in) :: Tlist
    real(8), dimension(1:Nvib), intent(out) :: vib

    real(8), dimension(1:0,1:2) :: hind

    if (Ntemp /= 7) then
        write (*,'(a)'), 'The number of temperatures must be exactly 7.'
        stop
    end if

    ! Fit the harmonic oscillator and hindered rotor modes to the heat capacity
    ! The function fitSpectralDataToHeatCapacity is in calc_freq_code.f90
    call fitSpectralDataToHeatCapacity(Cv, Tlist, Ntemp, Nvib, 0, vib, hind)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fitSpectralDataNoOscillators(Cv, Tlist, Ntemp, Nhind, hind)
	! Estimate the spectroscopic degrees of freedom by fitting parameters to
	! heat capacity data.
	!
	! ========== ====== ========================================================
	! Parameter  Intent Description
	! ========== ====== ========================================================
	! `Cv`       in     A list of heat capacities Cv/R at various temperatures
	!                   for the unknown degrees of freedom (i.e. the known
	!                   degrees of freedom should have already been removed)
	! `Tlist`    in     The temperatures corresponding to the heat capacities
	! `Ntemp`    in     The number of temperatures and heat capacities provided
	! `Nhind`    in     The number of 1D Pitzer hindered rotors to fit
	! `hind`     out    A matrix of fitted 1D Pitzer hindered rotor frequency-
	!                   barrier pairs, both in cm^-1
	! ========== ====== ========================================================

	integer, intent(in) :: Ntemp
	integer, intent(in) :: Nhind
	real(8), dimension(1:Ntemp), intent(in) :: Cv
	real(8), dimension(1:Ntemp), intent(in) :: Tlist
	real(8), dimension(1:Nhind,1:2), intent(out) :: hind

	real(8), dimension(1:0,1:2) :: vib

	if (Ntemp /= 7) then
		write (*,'(a)'), 'The number of temperatures must be exactly 7.'
		stop
	end if

	! Fit the harmonic oscillator and hindered rotor modes to the heat capacity
	! The function fitSpectralDataToHeatCapacity is in calc_freq_code.f90
	call fitSpectralDataToHeatCapacity(Cv, Tlist, Ntemp, 0, Nhind, vib, hind)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

