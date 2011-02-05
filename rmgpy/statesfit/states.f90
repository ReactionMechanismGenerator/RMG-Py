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

subroutine densityOfStates(E, Ngrains, vib, Nvib, rot, Nrot, hind, Nhind, &
    symm, linear, rho, msg)
    ! Determine the density of states at the specified energies. The parameters
    ! are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     The energies to determine the density of states at in
    !                   J/mol
    ! `Ngrains`  in     The number of energy grains
    ! `vib`      in     An array of harmonic oscillator frequencies in cm^-1
    ! `Nvib`     in     The number of harmonic oscillator modes
    ! `rot`      in     An array of rigid rotor frequencies in cm^-1
    ! `Nrot`     in     The number of rigid rotor modes
    ! `hind`     in     An array of hindered rotor frequencies and barrier
    !                   heights in cm^-1
    ! `Nhind`    in     The number of hindered rotor modes
    ! `symm`     in     The combined external + internal symmetry number
    ! `linear`   in     1 if the molecule is linear, 0 if nonlinear
    ! `rho`      out    The density of states at the specified energies in
    !                   mol/J
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: Ngrains
    integer, intent(in) :: Nvib
    integer, intent(in) :: Nrot
    integer, intent(in) :: Nhind
    real(8), dimension(1:Ngrains), intent(in) :: E
    real(8), dimension(1:Nvib), intent(in) :: vib
    real(8), dimension(1:Nrot), intent(in) :: rot
    real(8), dimension(1:Nhind,1:2), intent(in) :: hind
    integer, intent(in) :: symm
    integer, intent(in) :: linear
    real(8), dimension(1:Ngrains), intent(out) :: rho
    character(len=128), intent(out) :: msg

    ! The energy grain size
    real(8) :: dE
    ! A temporary density of states array for single modes to be convolved
    real(8), dimension(1:Ngrains) :: rho0
    ! Some integer indices
    integer r, i
    
    real(8) :: dE_vib
    integer :: mult, Ngrains_vib
    real(8), dimension(:), allocatable :: E_vib, rho_vib

    ! Set msg to successful; will be changed later if an error occurs
    msg = ""

    ! Zero the density of states output
    do r = 1, Ngrains
        rho(r) = 0.0
    end do

    ! Set the energy grain size
    dE = E(2) - E(1)

    ! Rigid rotor modes
    if (Nrot > 0) then
        call freeRotor_densityOfStates(E, Ngrains, rot, Nrot, linear, rho0)
        if (msg(1:1) /= ' ') return
        call convolve(rho, rho0, E, Ngrains)
    else
        ! If no free rotors than add an active K-rotor
        ! The rotational constant cancels out in subsequent calculations, so we
        ! arbitrarily choose 1.0
        ! This also helps to smooth the density of states
        call kRotor_densityOfStates(E, Ngrains, rho0)
        if (msg(1:1) /= ' ') return
        call convolve(rho, rho0, E, Ngrains)
    end if

    ! Hindered rotor modes
    do i = 1, Nhind
        call hinderedRotor_densityOfStates(E, Ngrains, hind(i,1), hind(i,2), rho0)
        if (msg(1:1) /= ' ') return
        call convolve(rho, rho0, E, Ngrains)
    end do
    
    ! Symmetry number
    do r = 1, Ngrains
        rho(r) = rho(r) / symm
    end do

    ! Harmonic oscillator modes
    if (Nvib > 0) then

        ! Must use grain size < 10 cm^-1 for convolution of vibrational modes
        ! into the density of states to ensure reasonable accuracy of the
        ! Beyer-Swinehart procedure
        dE_vib = dE
        mult = 1
        do while (dE_vib > 100.0)
            dE_vib = dE_vib / 2.0
            mult = mult * 2
        end do
        Ngrains_vib = (Ngrains-1)*mult + 1

        allocate( E_vib(1:Ngrains_vib), rho_vib(1:Ngrains_vib) )
        do r = 1, Ngrains_vib
            E_vib(r) = (r - 1) * dE_vib + minval(E)
            i = (r - 1) / mult + 1
            if (i == 1) then
                rho_vib(r) = 0.0
            elseif (r == mult + 1) then
                rho_vib(r) = rho(i)
            elseif (i >= Ngrains) then
                rho_vib(r) = rho(Ngrains)
            else
                rho_vib(r) = rho(i) * ( rho(i+1) / rho(i)  ) ** ((E_vib(r) - E(i)) / dE)
            end if
        end do

        ! The Beyer-Swinehart algorithm is an efficient way to convolve in
        ! the vibrational modes
        call beyerSwinehart(E_vib, Ngrains_vib, vib, Nvib, rho_vib, msg)
        if (msg(1:1) /= ' ') return

        rho(1:Ngrains) = rho_vib(1:Ngrains_vib:mult)

        deallocate( E_vib, rho_vib )

    end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine beyerSwinehart(E, Ngrains, vib, Nvib, rho, msg)
    ! Convolve vibrational modes into a density of states vector using the
    ! Beyer-Swinehart algorithm. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     The energies to determine the density of states at in
    !                   cm^-1
    ! `Ngrains`  in     The number of energy grains
    ! `vib`      in     An array of harmonic oscillator frequencies in cm^-1
    ! `Nvib`     in     The number of harmonic oscillator modes
    ! `rho`      in/out The density of states at the specified energies in
    !                   (cm^-1)^-1
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: Ngrains
    integer, intent(in) :: Nvib
    real(8), dimension(1:Ngrains), intent(in) :: E
    real(8), dimension(1:Nvib), intent(in) :: vib
    real(8), dimension(1:Ngrains), intent(inout) :: rho
    character(len=128), intent(out) :: msg

    real(8) dE
    integer i, n, dn

    msg = ''

    ! Set grain size
    dE = E(2) - E(1)

    ! The Beyer-Swinehart algorithm
    do i = 1, Nvib
        dn = nint(vib(i) * 6.626e-34 * 2.9979e10 * 6.022e23 / dE)
        do n = dn+1, Ngrains
            rho(n) = rho(n) + rho(n-dn)
        end do
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine beyerSwinehartSR(E, Ngrains, vib, Nvib, rho, msg)
    ! Convolve vibrational modes into a density of states vector using the
    ! Beyer-Swinehart algorithm with the Stein-Rabinovitch modification. The
    ! Stein-Rabinovitch modification is more accurate, but more expensive. The
    ! parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     The energies to determine the density of states at in
    !                   cm^-1
    ! `Ngrains`  in     The number of energy grains
    ! `vib`      in     An array of harmonic oscillator frequencies in cm^-1
    ! `Nvib`     in     The number of harmonic oscillator modes
    ! `rho`      in/out The density of states at the specified energies in
    !                   (cm^-1)^-1
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: Ngrains
    integer, intent(in) :: Nvib
    real(8), dimension(1:Ngrains), intent(in) :: E
    real(8), dimension(1:Nvib), intent(in) :: vib
    real(8), dimension(1:Ngrains), intent(inout) :: rho
    character(len=128), intent(out) :: msg

    real(8) dE
    integer i, n, dn, dn0
    real(8), dimension(1:Ngrains) :: rho0

    msg = ''
    
    ! Set grain size
    dE = E(2) - E(1)
    
    rho0 = rho

    ! The Beyer-Swinehart algorithm with the Stein-Rabinovitch modification
    do i = 1, Nvib
        dn0 = nint(vib(i) * 6.626e-34 * 2.9979e10 * 6.022e23 / dE)
        dn = dn0
        do while (dn < Ngrains)
            do n = 1, Ngrains-dn
                rho0(n+dn) = rho0(n+dn) + rho(n)
            end do
            dn = dn + dn0
        end do
        rho = rho0
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine convolve(rho1, rho2, E, Ngrains)
    ! Convolve two density of states vectors. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `rho1`     in     The first density of states vector in (cm^-1)^-1
    ! `rho2`     in     The second density of states vector in (cm^-1)^-1
    ! `E`        in     The energies in cm^-1
    ! `Ngrains`  in     The number of energy grains
    ! `convolve` out    The convolved density of states in (cm^-1)^-1
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: Ngrains
    real(8), dimension(1:Ngrains), intent(inout) :: rho1
    real(8), dimension(1:Ngrains), intent(in) :: rho2
    real(8), dimension(1:Ngrains), intent(in) :: E
    
    integer i, j
    integer found1, found2
    real(8) dE

    real(8), dimension(1:Ngrains) :: rho

    ! Set grain size
    dE = E(2) - E(1)

    ! Zero output
    do i = 1, Ngrains
        rho(i) = 0.0
    end do

    ! Check that both vectors have nonzero components
    found1 = 0
    found2 = 0
    do i = 1, Ngrains
        if (rho1(i) > 0.0) found1 = 1
        if (rho2(i) > 0.0) found2 = 1
    end do

    if (found1 == 0 .and. found2 == 0) then
        return
    elseif (found1 /= 0 .and. found2 == 0) then
        rho = rho1
    elseif (found1 == 0 .and. found2 /= 0) then
        rho = rho2
    else

        ! Do convolution
        do i = 1, Ngrains
            do j = 1, i
                rho(i) = rho(i) + rho2(i-j+1) * rho1(j) * dE
            end do
        end do

    end if

    rho1 = rho

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
