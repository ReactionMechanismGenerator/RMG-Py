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

module params

    implicit none

    integer :: nT
    real(8), dimension(:), pointer :: Tlist
    real(8), dimension(:), pointer :: Cvlist

    integer :: mcon
    integer :: mequa
    integer :: nvars

    integer :: nvib
    integer :: nrot
    
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_params(p_nT, p_Tlist, p_Cvlist, p_mcon, p_mequa, p_nvars, &
    p_nvib, p_nrot)
    ! Set the module parameters based on the values in the parameter list.

    use params
    implicit none
    
    integer, intent(in) :: p_nT                             ! The number of temperatures at which Cv data is given
    real(8), dimension(1:p_nT), intent(in) :: p_Tlist       ! The list of temperatures in K at which Cv data is given
    real(8), dimension(1:p_nT), intent(in) :: p_Cvlist      ! The list of heat capacity data (dimensionless) at each temperature
    integer, intent(in) :: p_mcon
    integer, intent(in) :: p_mequa
    integer, intent(in) :: p_nvars                          ! The number of variables used in the fitting
    integer, intent(in) :: p_nvib                           ! The number of harmonic oscillators to fit
    integer, intent(in) :: p_nrot                           ! The number of hindered rotors to fit
    
    nT = p_nT
    allocate( Tlist(1:nT), Cvlist(1:nT) )
    Tlist = p_Tlist
    Cvlist = p_Cvlist

    mcon = p_mcon
    mequa = p_mequa
    nvars = p_nvars
    
    nvib = p_nvib
    nrot = p_nrot

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine clean_up()
    ! Clean up arrays allocated when setting up the module parameters.

    use params
    implicit none
    
    deallocate( Tlist, Cvlist )
    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fit_modes(x0, nx, bl, bu, ind, maxiter, xout, igo)
    ! Execute a constrained nonlinear optimization to fit molecular degrees of
    ! freedom to heat capacity data.

    use params
    implicit none

    integer, intent(in) :: nx
    real(8), dimension(1:nx), intent(in) :: x0          ! The initial guess to the solver for each variable
    real(8), dimension(1:nx), intent(in) :: bl          ! The lower bounds for each variable
    real(8), dimension(1:nx), intent(in) :: bu          ! The lower bounds for each variable
    integer, dimension(1:nx), intent(in) :: ind         ! The type of bounds for each variable
    integer, intent(in) :: maxiter                      ! The maximum number of iterations to try
    real(8), dimension(1:nx), intent(out) :: xout       ! The returned solution vector
    integer, intent(out) :: igo                         ! The returned solution status

    real(8), dimension(1:nx) :: x

    external caseDirect
    external casePseudo
    external casePseudoRot

    ! These variables are required by the nonlinear solver
    integer, parameter :: lwork = 3000
    integer, parameter :: liwork = 1500
    real(8) :: fj(nT,nvars+1)
    real(8) :: fnorm
    integer :: iopt(24)
    integer :: iwork(liwork)
    real(8) :: ropt(1)
    real(8) :: work(lwork)
    integer :: ldfj

    ldfj = mcon + mequa

    ! Tell how much storage we gave the solver
    iwork(1) = lwork
    iwork(2) = liwork

    ! Additional solver options
    iopt(1)=4         ! Set the option to change the value of TOLF
    iopt(2)=1         ! Where in ROPT to look for the new value
    ropt(1)=1.E-5     ! New value for TOLF
    iopt(3)=2         ! Change the number of interations
    iopt(4)=maxiter   ! Maximum number of iterations
    iopt(5)=17        ! Do not allow the flag IGO to return the value IGO=3
    iopt(6)=1         ! Forces a full model step
    iopt(7)=99        ! No further options are changed

    ! x is used to store the input to dqed and is overwritten with the output
    ! from dqed, so it needs to be a temporary variable
    x = x0

    ! Execute the optimization
    ! We use a different helper function based on the number of oscillators
    ! and rotors we are trying to fit
    if (nvib + 2 * nrot < nT) then
        call dqed(caseDirect, mequa, nvars, mcon, ind, bl, bu, &
            x, fj, ldfj, fnorm, igo, iopt, ropt, iwork, work)
    elseif (nvib + 2 < nT) then
        call dqed(casePseudoRot, mequa, nvars, mcon, ind, bl, bu, &
            x, fj, ldfj, fnorm, igo, iopt, ropt, iwork, work)
    else
        call dqed(casePseudo, mequa, nvars, mcon, ind, bl, bu, &
            x, fj, ldfj, fnorm, igo, iopt, ropt, iwork, work)
    end if

    ! Store the results in the output vector
    ! We defer interpreting them until after returning to Python
    xout = x

end subroutine

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
                fj(mcon+j,i) = 2.0 * diff(j) * dCv(i,j)
            end do
        end do
    end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine casePseudo(x, fj, ldfj, igo, iopt, ropt)
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
    real(8), dimension(1:nT) :: dCv1, dCv2, dCv3, dCv4, dCv5

    integer i

    ! Calculate the heat capacities and derivatives at each temperature
    call harmonicOscillator_heatCapacity(Tlist, nT, x(1), Cv1)
    call harmonicOscillator_heatCapacity(Tlist, nT, x(3), Cv2)
    call harmonicOscillator_heatCapacity(Tlist, nT, x(4), Cv3)
    call hinderedRotor_heatCapacity(Tlist, nT, x(5), x(6), Cv4)
    call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(1), dCv1)
    call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(3), dCv2)
    call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(4), dCv3)
    call hinderedRotor_d_heatCapacity_d_freq(Tlist, nT, x(5), x(6), dCv4)
    call hinderedRotor_d_heatCapacity_d_barr(Tlist, nT, x(5), x(6), dCv5)

    ! There are no constraints other than the bounds, so we don't need to
    ! calculate the values of the constraints or their Jacobian

    ! There is a minimization objective for each temperature
    ! First calculate the least-squares objective at each temperature
    ! This goes in the last column of the fj matrix
    diff = Cv1 + x(2) * Cv2 + (nvib - x(2) - 1) * Cv3 + Nrot * Cv4 - Cvlist
    fj(mcon+1:mcon+mequa,nvars+1) = diff * diff
    
    ! We also need the Jacobian of the objective functions if igo is nonzero
    if (igo /= 0) then
        do i = 1, mequa
            ! The first harmonic oscillator frequency
            fj(mcon+i,1) = 2.0 * diff(i) * dCv1(i)
            ! The second harmonic oscillator degeneracy
            fj(mcon+i,2) = 2.0 * diff(i) * (Cv2(i) - Cv3(i))
            ! The second harmonic oscillator pseudo-frequency
            fj(mcon+i,3) = 2.0 * diff(i) * (x(2) * dCv2(i))
            ! The third harmonic oscillator pseudo-frequency
            fj(mcon+i,4) = 2.0 * diff(i) * ((real(nvib) - x(2) - 1) * dCv3(i))
            ! The hindered rotor pseudo-frequency
            fj(mcon+i,5) = 2.0 * diff(i) * Nrot * dCv4(i)
            ! The hindered rotor pseudo-barrier
            fj(mcon+i,6) = 2.0 * diff(i) * Nrot * dCv5(i)
        end do
    end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine casePseudoRot(x, fj, ldfj, igo, iopt, ropt)
    ! Objective function and its Jacobian for the case of many real oscillators
    ! and many pseudo-rotors.

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
    real(8), dimension(1:Nvib+1,1:nT) :: Cv
    real(8), dimension(1:Nvib+2,1:nT) :: dCv

    integer i, j

    ! Calculate the heat capacities and derivatives at each temperature
    do i = 1, Nvib
        call harmonicOscillator_heatCapacity(Tlist, nT, x(i), Cv(i,:))
        call harmonicOscillator_d_heatCapacity_d_freq(Tlist, nT, x(i), dCv(i,:))
    end do
    call hinderedRotor_heatCapacity(Tlist, nT, x(Nvib+1), x(Nvib+2), Cv(Nvib+1,:))
    call hinderedRotor_d_heatCapacity_d_freq(Tlist, nT, x(Nvib+1), x(Nvib+2), dCv(Nvib+1,:))
    call hinderedRotor_d_heatCapacity_d_barr(Tlist, nT, x(Nvib+1), x(Nvib+2), dCv(Nvib+2,:))
    Cv(Nvib+1,:) = Nrot * Cv(Nvib+1,:)
    dCv(Nvib+1,:) = Nrot * dCv(Nvib+1,:)
    dCv(Nvib+2,:) = Nrot * dCv(Nvib+2,:)

    ! There are no constraints other than the bounds, so we don't need to
    ! calculate the values of the constraints or their Jacobian

    ! There is a minimization objective for each temperature
    ! First calculate the least-squares objective at each temperature
    ! This goes in the last column of the fj matrix
    diff = -Cvlist
    do i = 1, Nvib+1
        diff = diff + Cv(i,:)
    end do
    fj(mcon+1:mcon+mequa,nvars+1) = diff * diff

    ! We also need the Jacobian of the objective functions if igo is nonzero
    if (igo /= 0) then
        do j = 1, mequa
            do i = 1, Nvib+2
                fj(mcon+j,i) = 2.0 * diff(j) * dCv(i,j)
            end do
        end do
    end if

end subroutine