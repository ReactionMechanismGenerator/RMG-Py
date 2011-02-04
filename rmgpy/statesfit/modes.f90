!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   RMG - Reaction Mechanism Generator
!
!   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
!   RMG Team (rmg_dev@mit.edu)
!
!   Copyright (c) 2009 by Josh Allen (chemejosh@gmail.com)
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

! Contains functions for evaluating the density of states, partition function,
! heat capacity, and heat capacity derivatives for various molecular degrees of
! freedom. This code is written in Fortran for speed and for accessibility to
! Fortran routines, such as those used to fit degrees of freedom to heat
! capacity data.
!
! Additional details on the functions themselves, including the formulas
! implemented, can be found in the Python module most directly associated with
! this file, modes.py.
!
! The functions in this file are grouped by mode: translation, free rotor,
! harmonic oscillator, hindered rotor.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine translation_densityOfStates(Elist, nE, mass, dim, V, rho)

    implicit none
    integer, intent(in) :: nE                           ! The number of energies
    real(8), dimension(1:nE), intent(in) :: Elist       ! The energies at which to evaluate in J/mol
    real(8), intent(in) :: mass                         ! The species mass in kg/mol
    integer, intent(in) :: dim                          ! The dimensionality of the translation (2 or 3)
    real(8), intent(in) :: V                            ! The system volume in m^3
    real(8), dimension(1:nE), intent(out) :: rho        ! The evaluated density of states in mol/J

    integer i
    real(8) qt

    qt = ((2 * 3.141592654 * mass) / (6.626e-34 * 6.626e-34))**(dim/2.0) * V
    if (dim == 2) then
        ! Dimensionality is 2
        do i = 1, nE
            rho(i) = qt
        end do
    elseif (dim == 3) then
        ! Dimensionality is 3
        do i = 1, nE
            rho(i) = qt * sqrt(Elist(i)) / sqrt(3.141592694) / 2    ! qt * E^0.5 / 0.5!
        end do
    else
        ! Unexpected dimensionality, so return zero
        write (*,*) 'Unexpected dimensionality encountered:', dim
        rho = 0 * rho
    end if

end subroutine

subroutine translation_partitionFunction(Tlist, nT, mass, dim, V, Q)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: mass                         ! The species mass in kg/mol
    integer, intent(in) :: dim                          ! The dimensionality of the translation (2 or 3)
    real(8), intent(in) :: V                            ! The system volume in m^3
    real(8), dimension(1:nT), intent(out) :: Q          ! The evaluated partition function in m^3

    integer i
    real(8) qt

    qt = ((2 * 3.141592654 * mass) / (6.626e-34 * 6.626e-34))**(dim/2.0) * V

    do i = 1, nT
        Q(i) = qt * (8.314472 * Tlist(i))**(dim/2.0)
    end do

end subroutine

subroutine translation_heatCapacity(Tlist, nT, mass, dim, V, Cvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: mass                         ! The species mass in kg/mol
    integer, intent(in) :: dim                          ! The dimensionality of the translation (2 or 3)
    real(8), intent(in) :: V                            ! The system volume in m^3
    real(8), dimension(1:nT), intent(out) :: Cvlist     ! The evaluated heat capacity in J/mol*K

    integer i

    do i = 1, nT
        Cvlist(i) = 0.5 * dim
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine freeRotor_densityOfStates(Elist, nE, freq, nfreq, linear, rho)

    implicit none
    integer, intent(in) :: nE                           ! The number of energies
    integer, intent(in) :: nfreq                        ! The number of rotor frequencies
    real(8), dimension(1:nE), intent(in) :: Elist       ! The list of energies to evaluate in J/mol
    real(8), dimension(1:nfreq), intent(in) :: freq     ! The list of frequencies in cm^-1
    integer, intent(in) :: linear                       ! 0 = nonlinear, 1 = linear
    real(8), dimension(1:nE), intent(out) :: rho        ! The evaluated density of states in mol/J

    integer i
    real(8) theta

    if (linear /= 0) then
        theta = 6.626e-34 * 2.9979e10 * freq(1) * 6.022e23
        do i = 1, nE
            rho(i) = 1.0 / theta
        end do
    else
        theta = 1.0
        do i = 1, nfreq
            theta = theta * 6.626e-34 * 2.9979e10 * freq(i) * 6.022e23
        end do
        do i = 1, nE
            rho(i) = 2.0 * sqrt(Elist(i) / theta)
        end do
    end if

end subroutine

subroutine freeRotor_partitionFunction(Tlist, nT, freq, nfreq, linear, Q)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    integer, intent(in) :: nfreq                        ! The number of rotor frequencies
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), dimension(1:nfreq), intent(in) :: freq     ! The list of frequencies in cm^-1
    integer, intent(in) :: linear                       ! 0 = nonlinear, 1 = linear
    real(8), dimension(1:nT), intent(out) :: Q          ! The evaluated partition function in m^3

    integer i
    real(8) theta

    if (linear /= 0) then
        theta = 6.626e-34 * 2.9979e10 * freq(1) / 1.381e-23
        do i = 1, nT
            Q(i) = Tlist(i) / theta
        end do
    else
        theta = 1.0
        do i = 1, nfreq
            theta = theta * 6.626e-34 * 2.9979e10 * freq(i) / 1.381e-23
        end do
        do i = 1, nT
            Q(i) = sqrt(3.141592654 * Tlist(i)**nfreq / theta)
        end do
    end if

end subroutine

subroutine freeRotor_heatCapacity(Tlist, nT, freq, nfreq, linear, Cvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    integer, intent(in) :: nfreq                        ! The number of rotor frequencies
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), dimension(1:nfreq), intent(in) :: freq     ! The list of frequencies in cm^-1
    integer, intent(in) :: linear                       ! 0 = nonlinear, 1 = linear
    real(8), dimension(1:nT), intent(out) :: Cvlist     ! The evaluated heat capacity in J/mol*K

    integer i

    if (linear /= 0) then
        do i = 1, nT
            Cvlist(i) = 1.0
        end do
    else
        do i = 1, nT
            Cvlist(i) = 1.5
        end do
    end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine harmonicOscillator_partitionFunction(Tlist, nT, freq, Q)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The oscillator frequency in cm^-1
    real(8), dimension(1:nT), intent(out) :: Q          ! The evaluated partition function in m^3

    integer i
    real(8) xi

    do i = 1, nT
        xi = 6.626e-34 * 2.9979e10 * freq / (1.381e-23 * Tlist(i))
        Q(i) = 1.0 / (1 - exp(-xi))
    end do

end subroutine

subroutine harmonicOscillator_heatCapacity(Tlist, nT, freq, Cvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The oscillator frequency in cm^-1
    real(8), dimension(1:nT), intent(out) :: Cvlist     ! The evaluated heat capacity in J/mol*K

    integer i
    real(8) x, exp_x, one_minus_exp_x

    do i = 1, nT
        x = freq / (0.695039 * Tlist(i))        ! kB = 0.695039 cm^-1/K
        
        exp_x = exp(x)
        one_minus_exp_x = 1.0 - exp_x

        Cvlist(i) = x * x * exp_x / one_minus_exp_x / one_minus_exp_x
    end do

end subroutine

subroutine harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, freq, dCvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The oscillator frequency in cm^-1
    real(8), dimension(1:nT), intent(out) :: dCvlist    ! The evaluated heat capacity partial derivative in J/mol*K/cm^-1

    integer i
    real(8) x, exp_x, one_minus_exp_x

    do i = 1, nT
        x = freq / (0.695039 * Tlist(i))        ! kB = 0.695039 cm^-1/K

        exp_x = exp(x)
        one_minus_exp_x = 1.0 - exp_x

        dCvlist(i) = x * exp_x / one_minus_exp_x / one_minus_exp_x * &
            (2.0 + x + 2.0 * x * exp_x / one_minus_exp_x) * x / freq
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hinderedRotor_densityOfStates(Elist, nE, freq, barr, rho)

    integer, intent(in) :: nE                           ! The number of energies
    real(8), dimension(1:nE), intent(in) :: Elist       ! The energies at which to evaluate in J/mol
    real(8), intent(in) :: freq                         ! The rotor frequency in cm^-1
    real(8), intent(in) :: barr                         ! The rotor barrier height in cm^-1
    real(8), dimension(1:nE), intent(out) :: rho        ! The evaluated density of states in mol/J

    integer i
    real(8) pre, V0
    real(8) cellipk

    ! The following is only valid in the classical limit

    pre = 2.0 / 3.141592654 / (6.626e-34 * 2.9979e10 * freq * 6.022e23)
    V0 = 6.626e-34 * 2.9979e10 * barr * 6.022e23

    do i = 1, nE
        if (Elist(i) / V0 < 1) then
            rho(i) = pre * cellipk(Elist(i) / V0)
        else
            rho(i) = pre * sqrt(V0 / Elist(i)) * cellipk(V0 / Elist(i))
        end if
    end do

end subroutine

subroutine hinderedRotor_partitionFunction(Tlist, nT, freq, barr, Q)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The rotor frequency in cm^-1
    real(8), intent(in) :: barr                         ! The rotor barrier height in cm^-1
    real(8), dimension(1:nT), intent(out) :: Q          ! The evaluated partition function in m^3

    integer i
    real(8) x, z
    real(8) besseli0

    do i = 1, nT
        x = 6.626e-34 * 2.9979e10 * freq / 1.381e-23 / Tlist(i)
        z = 0.5 * 6.626e-34 * 2.9979e10 * barr / 1.381e-23 / Tlist(i)
        ! The following is only valid in the classical limit
        Q(i) = sqrt(2.0 * 3.141592654 * z) / x * exp(-z) * besseli0(z)
    end do

end subroutine

subroutine hinderedRotor_heatCapacity(Tlist, nT, freq, barr, Cvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The rotor frequency in cm^-1
    real(8), intent(in) :: barr                         ! The rotor barrier height in cm^-1
    real(8), dimension(1:nT), intent(out) :: Cvlist     ! The evaluated heat capacity in J/mol*K

    integer i
    real(8) x, z, exp_x, one_minus_exp_x, BB
    real(8) besseli1, besseli0

    do i = 1, nT
        x = 6.626e-34 * 2.9979e10 * freq / 1.381e-23 / Tlist(i)
        z = 0.5 * 6.626e-34 * 2.9979e10 * barr / 1.381e-23 / Tlist(i)
        
        exp_x = exp(x)
        one_minus_exp_x = 1.0 - exp_x
        BB = besseli1(z) / besseli0(z)
        
        Cvlist(i) = x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + &
            z * (z - BB - z * BB * BB)
            
    end do

end subroutine

subroutine hinderedRotor_d_heatCapacity_d_freq(Tlist, nT, freq, barr, dCvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The rotor frequency in cm^-1
    real(8), intent(in) :: barr                         ! The rotor barrier height in cm^-1
    real(8), dimension(1:nT), intent(out) :: dCvlist    ! The evaluated heat capacity partial derivative in J/mol*K/cm^-1

    integer i
    real(8) x, exp_x, one_minus_exp_x
    
    do i = 1, nT
        x = 6.626e-34 * 2.9979e10 * freq / 1.381e-23 / Tlist(i)
        
        exp_x = exp(x)
        one_minus_exp_x = 1.0 - exp_x
        
        dCvlist(i) = x * exp_x / one_minus_exp_x / one_minus_exp_x * &
            (2 + x + 2 * x * exp_x / one_minus_exp_x) * x / freq
    end do

end subroutine

subroutine hinderedRotor_d_heatCapacity_d_barr(Tlist, nT, freq, barr, dCvlist)

    implicit none
    integer, intent(in) :: nT                           ! The number of temperatures
    real(8), dimension(1:nT), intent(in) :: Tlist       ! The list of temperatures to evaluate in K
    real(8), intent(in) :: freq                         ! The rotor frequency in cm^-1
    real(8), intent(in) :: barr                         ! The rotor barrier height in cm^-1
    real(8), dimension(1:nT), intent(out) :: dCvlist    ! The evaluated heat capacity partial derivative in J/mol*K/cm^-1

    integer i
    real(8) z, BB
    real(8) besselratio
    
    do i = 1, nT
        z = 0.5 * 6.626e-34 * 2.9979e10 * barr / 1.381e-23 / Tlist(i)

        BB = besselratio(z)      !besseli(1, z) / besseli(0, z)

        dCvlist(i) = z * (1 - 2 * z * BB + BB * BB + 2 * z * BB * BB * BB) * z / barr
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kRotor_densityOfStates(E, Ngrains, rho)
    ! Determine the density of states for an arbitrary active K-rotor at the
    ! specified energies. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     The energies to determine the density of states at in
    !                   cm^-1
    ! `Ngrains`  in     The number of energy grains
    ! `rho`      out    The density of states at the specified energies in
    !                   (cm^-1)^-1
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer :: Ngrains
    real(8), dimension(1:Ngrains), intent(in) :: E
    real(8), dimension(1:Ngrains), intent(out) :: rho
    
    integer r

    do r = 1, Ngrains
        if (E(r) == 0) then
            rho(r) = 0.0
        else
            rho(r) = 1.0 / sqrt(1.0 * E(r))
        end if
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
