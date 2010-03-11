!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   MEMURN - Master Equation Model for Unimolecular Reaction Networks
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

subroutine collisionMatrix(T, P, E, collFreq, densStates, E0, dEdown, Ngrains, &
    Mcoll, msg)
    ! Calculates the collisional transfer rate matrix for a given isomer.
    !
    ! ============ ====== ======================================================
    ! Parameter    Intent Description
    ! ============ ====== ======================================================
    ! `T`          in     The temperature to evaluate k(T,P) at in K
    ! `P`          in     The pressure to evaluate k(T,P) at in Pa
    ! `E`          in     A 1D array of energies in J/mol
    ! `collFreq`   in     The (modified) collision frequency for the isomer in
    !                     Hz
    ! `densStates` in     The density of states for the isomer
    ! `E0`         in     The ground-state energy for the isomer in J/mol
    ! `dEdown`     in     The average energy transferred in a deactivating
    !                     collision in J/mol
    ! `Ngrains`    in     The number of energy grains being used
    ! `Mcoll`      out    The collision matrix for the isomer
    ! `msg`        out    If the subroutine was unsuccessful, this string will
    !                     contain a brief message describing the error; the
    !                     string will be empty if the subroutine was successful
    ! ============ ====== ======================================================

    ! Type definitions of parameters
    integer, intent(in) :: Ngrains
    real(8), intent(in) :: T
    real(8), intent(in) :: P
    real(8), dimension(1:Ngrains), intent(in) :: E
    real(8), intent(in) :: collFreq
    real(8), dimension(1:Ngrains), intent(in) :: densStates
    real(8), intent(in) :: E0
    real(8), intent(in) :: dEdown
    real(8), dimension(1:Ngrains,1:Ngrains), intent(out) :: Mcoll
    character(len=128), intent(out) :: msg
    
    real(8), dimension(:), allocatable      ::  C       ! Vector of normalization coefficients

    real(8) dE

    integer         :: halfbandwidth
    integer         :: lb, ub, start
    integer         :: r, s

    ! Initialize msg to empty string
    msg = ''
    
    ! Allocate/zero vector of normalization coefficients
    allocate( C(1:Ngrains) )
    do i = 1, Ngrains
        C(i) = 0.0
    end do
    
    ! Get minimum and maximum energies and grain size
    dE = E(2) - E(1)

    ! Determine bandwidth (at which transfer probabilities are so low that they can be truncated
    ! with negligible error)
    halfbandwidth = ceiling(16 * dEdown / dE)
    
    ! Determine start grain (corresponding to isomer ground-state energy)
    start = 0
    do r = 1, Ngrains
        if (start == 0 .and. densStates(r) > 0) then
            start = r
        end if
    end do
    if (start == 0) then
        msg = 'Unable to determine starting energy grain.'
        return
    end if

    ! Determine unnormalized entries in collisional tranfer probability matrix for the current isomer
    do r = 1, Ngrains
        do s = 1, Ngrains
            lb = max(r - halfbandwidth, start)
            ub = min(r + halfbandwidth, Ngrains)
            if (r >= start .and. s >= lb .and. s <= ub) then
                call transferRate(s, r, E(s), E(r), dEdown, E0, densStates, Ngrains, T, Mcoll(s,r))
            else
                Mcoll(s,r) = 0.0
            end if
        end do
    end do

    ! Normalize using detailed balance
    do r = start, Ngrains
        C(r) = (1 - sum(Mcoll(start:r-1,r))) / sum(Mcoll(r:Ngrains,r))
        ! Check for normalization consistency (i.e. all numbers are positive)
        if (C(r) <= 0) then
            msg = 'Error normalizing collisional transfer probabilities matrix.'
            return
        end if
        Mcoll(r,r+1:Ngrains) = Mcoll(r,r+1:Ngrains) * C(r)
        Mcoll(r:Ngrains,r) = Mcoll(r:Ngrains,r) * C(r)
        Mcoll(r,r) = Mcoll(r,r) - 1
    end do

    ! Normalize using detailed balance
    ! This is the way described by Pilling and Holbrook
!    do r = Ngrains, start, -1
!        C(r) = (1 - sum(Mcoll(r+1:Ngrains,r))) / sum(Mcoll(1:r,r))
!        ! Check for normalization consistency (i.e. all numbers are positive)
!        if (C(r) <= 0) then
!            msg = 'Error normalizing collisional transfer probabilities matrix.'
!            return
!        end if
!        ! Apply normalization condition
!        Mcoll(r,1:r-1) = Mcoll(r,1:r-1) * C(r)
!        Mcoll(1:r,r) = Mcoll(1:r,r) * C(r)
!        Mcoll(r,r) = Mcoll(r,r) - 1
!    end do

    ! DEBUG: Check that both our constraints are satisfied
!    do r = 1, Ngrains
!        if (sum(Mcoll(:,r)) > 0.0001) then
!            msg = 'Error: Column not normalized properly!'
!            return
!        end if
!    end do
!    do r = start, Ngrains
!        do s = start, r - 1
!            if (abs(Mcoll(r,s) * densStates(s) * exp(-E(s) / 8.314472 / T) &
!                - Mcoll(s,r) * densStates(r) * exp(-E(r) / 8.314472 / T)) > &
!                0.0001 * Mcoll(r,s) * densStates(s) * exp(-E(s) / 8.314472 / T)) then
!                msg = 'Error: Detailed balance not satisfied!'
!                return
!            end if
!        end do
!    end do

    ! Multiply by collision frequency to determine collision rates
    Mcoll = collFreq * Mcoll

    ! Clean up
    deallocate( C )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transferRate(i, j, Ei, Ej, alpha, E0, rho, Ngrains, T, rate)
    ! Computes the unnormalized probability of a collisional energy transfer
    ! from state j (with energy Ej) to state i (with energy Ei) using a
    ! single-exponental-down model.
    !
    ! ============ ====== ======================================================
    ! Parameter    Intent Description
    ! ============ ====== ======================================================
    ! `i`          in     index of destination state
    ! `j`          in     index of source state
    ! `Ei`         in     energy of destination state i in cm^-1
    ! `Ej`         in     energy of source state j in cm^-1
    ! `alpha`      in     parameter in exponential-down model in (cm^-1)^-1
    ! `E0`         in     electronic + zero-point energy of ground state in
    !                     cm^-1
    ! `rho`        in     density of states
    ! `Ngrains`    in     number of energy grains in density of states
    ! `T`          in     temperature of interest in K
    ! `rate`       out    collisional transfer probability in s^-1
    ! ============ ====== ======================================================

    ! Type definitions of parameters
    integer, intent(in)                 ::  Ngrains
    integer, intent(in)                 ::  i
    integer, intent(in)                 ::  j
    real(8), intent(in)                 ::  Ei
    real(8), intent(in)                 ::  Ej
    real(8), intent(in)                 ::  alpha
    real(8), intent(in)                 ::  E0
    real(8), dimension(1:Ngrains), intent(in)   ::  rho
    real(8), intent(in)                 ::  T
    real(8), intent(out)                ::  rate

    real(8) R

    R = 8.314472

    ! Evaluate collisional transfer probability - theoretically correct way
    if (Ej < E0 .or. Ei < E0 .or. rho(j) .eq. 0) then
       rate = 0.
    elseif (Ej >= Ei) then
       rate = exp(-(Ej - Ei) / alpha)
    else
       rate = exp(-(Ei - Ej) / alpha) * rho(i) / rho(j) * exp( -(Ei - Ej) / (R * T))
    end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fullMEMatrix(E, E0, Mcoll0, Kij, Gnj, Fim, indices, &
    nRows, nGrains, nIsom, nReac, nProd, Mcoll, Mrxn, msg)
    ! Construct the full master equation matrix. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     A 1D array of energies in J/mol
    ! `E0`       in     The active-state energy cutoffs for each isomer in J/mol
    ! `Mcoll0`   in     The collision matrices for each isomer
    ! `Kij`      in     The microcanonical isomerization rate coefficients in
    !                   s^-1
    ! `Fim`      in     The microcanonical association rates (rate coefficients
    !                   times bimolecular equilibrium distributions) in s^-1
    ! `Gnj`      in     The microcanonical dissociation rate coefficients in
    !                   s^-1
    ! `indices`  in     The indexing scheme to use
    ! `nRows`    in     The number of rows and columns in the full matrix
    ! `nGrains`  in     The number of energy grains being used
    ! `nIsom`    in     The number of isomers in the network
    ! `nReac`    in     The number of reactant channels in the network (both A + B <=> C)
    ! `nProd`    in     The number of product channels in the network (A -> B + C only)
    ! `Mcoll`    out    The full master equation matrix - collision terms
    ! `Mrxn`     out    The full master equation matrix - reaction terms
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: nRows
    integer, intent(in) :: nIsom
    integer, intent(in) :: nReac
    integer, intent(in) :: nProd
    integer, intent(in) :: nGrains
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), dimension(1:nIsom+nReac), intent(in) :: E0
    real(8), dimension(1:nIsom,1:nGrains,1:nGrains), intent(in) :: Mcoll0
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nReac,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    integer, dimension(1:nGrains,1:nIsom), intent(in) :: indices
    real(8), dimension(1:nRows,1:nRows), intent(out) :: Mcoll
    real(8), dimension(1:nRows,1:nRows), intent(out) :: Mrxn
    character(len=128), intent(out) :: msg

    integer i, j, n, r, s, u, v

    ! Initialize msg to empty string
    msg = ''

    ! Allocate and zero collisional and reactive matrices
    do i = 1, nRows
        do j = 1, nRows
            Mcoll(i,j) = 0.0
            Mrxn(i,j) = 0.0
        end do
    end do

    ! Construct collisional matrix
    do i = 1, nIsom
        do r = 1, nGrains
            if (indices(r,i) > 0) then
                do s = 1, nGrains
                    if (indices(s,i) > 0) then
                        Mcoll(indices(r,i), indices(s,i)) = Mcoll0(i,r,s)
                    end if
                end do
            end if
        end do
    end do
    
    ! Construct reactive matrix - isomerization reactions
    do i = 1, nIsom
        do j = 1, i-1
            if (Kij(i,j,nGrains) > 0 .or. Kij(j,i,nGrains) > 0) then
                do r = 1, nGrains
                    u = indices(r,i)
                    v = indices(r,j)
                    if (u > 0 .and. v > 0) then
                        Mrxn(v, u) = Kij(j,i,r)
                        Mrxn(u, u) = Mrxn(u, u) - Kij(j,i,r)
                        Mrxn(u, v) = Kij(i,j,r)
                        Mrxn(v, v) = Mrxn(v, v) - Kij(i,j,r)
                    end if
                end do
            end if
        end do
    end do
    ! Construct reactive matrix - dissocation/association reactions
    do i = 1, nIsom
        do n = 1, nReac+nProd
            if (Gnj(n,i,nGrains) > 0 .or. Fim(i,n,nGrains) > 0) then
                do r = 1, nGrains
                    u = indices(r,i)
                    v = nRows - nReac + n
                    if (u > 0) then
                        Mrxn(u, u) = Mrxn(u, u) - Gnj(n,i,r)
                        if (n <= nReac) then
                            Mrxn(v, u) = Gnj(n,i,r)
                            Mrxn(u, v) = Fim(i,n,r)
                            Mrxn(v, v) = Mrxn(v, v) - Fim(i,n,r)
                        end if
                    end if
                end do
            end if
        end do
    end do
    
    ! Clean up

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
