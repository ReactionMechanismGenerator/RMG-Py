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

subroutine estimateRateCoefficients_RS(T, P, E, Mcoll, densStates, E0, Eres, &
Kij, Fim, Gnj, dEdown, nIsom, nReac, nProd, nGrains, K, msg)
    ! Estimate the phenomenological rate coefficients using the (modified) strong
    ! collision method. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `T`        in     The temperature to evaluate k(T,P) at in K
    ! `P`        in     The pressure to evaluate k(T,P) at in Pa
    ! `E`        in     A 1D array of energies in J/mol
    ! `Mcoll`    in     The collision matrix for each isomer
    ! `densStates` in     The density of states for each isomer
    ! `E0`       in     The ground-state energy for each isomer in J/mol
    ! `Eres`     in     The active-state energy cutoffs for each isomer in J/mol
    ! `Kij`      in     The microcanonical isomerization rate coefficients in
    !                   s^-1
    ! `Fim`      in     The microcanonical association rates (rate coefficients
    !                   times bimolecular equilibrium distributions) in s^-1
    ! `Gnj`      in     The microcanonical dissociation rate coefficients in
    !                   s^-1
    ! `dEdown`   in     The average energy transferred in a deactivating
    !                   collision in J/mol
    ! `nIsom`    in     The number of isomers in the network
    ! `nReac`    in     The number of reactant channels in the network (both A + B <=> C)
    ! `nProd`    in     The number of product channels in the network (A -> B + C only)
    ! `nGrains`  in     The number of energy grains being used
    ! `K`        out    The matrix of phenomenological rate coefficients k(T,P)
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    real(8), intent(in) :: T
    real(8), intent(in) :: P
    integer, intent(in) :: nIsom
    integer, intent(in) :: nReac
    integer, intent(in) :: nProd
    integer, intent(in) :: nGrains
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), dimension(1:nIsom,1:nGrains,1:nGrains), intent(in) :: Mcoll
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    real(8), dimension(1:nIsom), intent(in) :: E0
    real(8), dimension(1:nIsom+nReac+nProd), intent(in) :: Eres
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nReac,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    real(8), intent(in) :: dEdown
    real(8), dimension(1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K
    character(len=128), intent(out) :: msg

    ! Number of reservoir and active-state energy grains for each isomer
    integer, dimension(:), allocatable          ::  nRes, nAct

    ! Pseudo-steady state grain populations
    real(8), dimension(:,:,:), allocatable  ::  pa
    
    ! Indices
    integer i, n, r

    ! Determine reservoir cutoff grains for each unimolecular isomer
    allocate( nRes(1:nIsom), nAct(1:nIsom) )
    call reservoirCutoffs(E0(1:nIsom), Eres, nIsom, E, nGrains, dEdown, densStates, nRes)
    do i = 1, nIsom
        nAct(i) = nGrains - nRes(i)
    end do

    ! Determine pseudo-steady state populations of active state
    allocate( pa(1:nGrains, 1:nIsom+nReac, 1:nIsom) )
    pa = 0 * pa
    if (nIsom == 1) then
        ! Not worth it to call banded solver when only one well
        call activeStateFull(T, P, E, Mcoll, densStates, Kij, Fim, Gnj, &
            nIsom, nReac, nProd, nGrains, nRes, nAct, pa, msg)
        if (msg(1:1) /= ' ') return
    else
        ! Very worth it to call banded solver when more than one well
        call activeStateBanded(T, P, E, Mcoll, densStates, Kij, Fim, Gnj, dEdown, &
            nIsom, nReac, nProd, nGrains, nRes, nAct, pa, msg)
        if (msg(1:1) /= ' ') return
    end if

    ! Check that PSSA populations are all nonnegative; fail if not
    do i = 1, nIsom
        do n = 1, nIsom+nReac
            do r = nRes(i), nGrains
                if (pa(r,n,i) < 0.0) then
                    msg = 'One or more negative steady-state populations encountered.'
                    return
                end if
            end do
        end do
    end do

    ! Initialize phenomenological rate coefficient matrix
    do i = 1, nIsom+nReac+nProd
        do n = 1, nIsom+nReac+nProd
            K(i,n) = 0
        end do
    end do

    ! Determine phenomenological rate coefficients
    do i = 1, nIsom
        do n = 1, nIsom+nReac
            do r = 1, nRes(i)
                K(i,n) = K(i,n) + sum(Mcoll(i,r,nRes(i)+1:nGrains) * pa(nRes(i)+1:nGrains,n,i))
            end do
        end do
    end do
    
    do j = 1, nIsom+nReac
        do n = 1, nReac+nProd
            do i = 1, nIsom
                K(n+nIsom,j) = K(n+nIsom,j) + sum(Gnj(n,i,nRes(i)+1:nGrains) * pa(nRes(i)+1:nGrains,j,i))
            end do
        end do
    end do

    do n = 1, nIsom+nReac+nProd
        K(n,n) = -sum(K(1:n-1,n)) - sum(K(n+1:nIsom+nReac+nProd,n))
    end do

    ! Clean up
    deallocate( nRes, nAct, pa )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reservoirCutoffs(E0, Eres, nIsom, E, nGrains, dEdown, densStates, nRes)
    ! Determines the grain below which the reservoir approximation will be used
    ! and above which the pseudo-steady state approximation will be used by
    ! examining the energies of the transition states connected to each
    ! unimolecular isomer. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E0`       in     The ground-state energy for each isomer in J/mol
    ! `Eres`     in     The active-state energy cutoffs for each isomer in J/mol
    ! `nIsom`    in     The number of isomers in the network
    ! `E`        in     A 1D array of energies in J/mol
    ! `nGrains`  in     The number of energy grains being used
    ! `dEdown`   in     The average energy transferred in a deactivating
    !                   collision in J/mol
    ! `nRes`     out    The number of reservoir grains for each isomer
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: nIsom
    integer, intent(in) :: nGrains
    real(8), dimension(1:nIsom), intent(in) :: E0
    real(8), dimension(1:nIsom), intent(in) :: Eres
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), intent(in) :: dEdown
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    integer, dimension(1:nIsom), intent(out) :: nRes

    real(8) Emin, dE
    integer, dimension(1:nIsom) :: start
    integer i, r

    Emin = minval(E)
    dE = E(2) - E(1)

    ! Determine reservoir cutoffs by looking at transition state energies
    do i = 1, nIsom
        ! Find the ground-state energy grain for this isomer
        start(i) = 0
        do r = 1, nGrains
            if (densStates(i,r) > 0 .and. start(i) == 0) then
                start(i) = r
            end if
        end do
        ! Now find the reservoir cutoff grain for this isomer
        nRes(i) = ceiling((Eres(i) - 10 * dEdown - Emin) / dE)
        ! Sometimes the above will result in nRes(i) < start (i.e. a cutoff
        ! below the ground state); we need to handle this if it happends
        if (nRes(i) < start(i)) then
            ! First try to place the cutoff a few grains below the lowest
            ! transition state energy
            nRes(i) = ceiling((Eres(i) - Emin) / dE) - 4
            if (nRes(i) < start(i)) then
                ! If this is still too low, then just put the cutoff a few
                ! grains above the ground state
                nRes(i) = start(i) + 4
            end if
        end if
    end do

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine activeStateFull(T, P, E, Mcoll, densStates, Kij, Fim, Gnj, &
    nIsom, nReac, nProd, nGrains, nRes, nAct, pa, msg)
    ! Determine the pseudo-steady state populations for the active state
    ! grains using a full matrix linear solve.
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     A 1D array of energies in J/mol
    ! `Mcoll`    in     The collision matrix for each isomer
    ! `densStates`   in     The density of states for each isomer
    ! `Kij`      in     The microcanonical isomerization rate coefficients in
    !                   s^-1
    ! `Fim`      in     The microcanonical association rates (rate coefficients
    !                   times bimolecular equilibrium distributions) in s^-1
    ! `Gnj`      in     The microcanonical dissociation rate coefficients in
    !                   s^-1
    ! `nIsom`    in     The number of isomers in the network
    ! `nReac`    in     The number of reactant channels in the network
    ! `nGrains`  in     The number of energy grains being used
    ! `nRes`     out    The number of reservoir grains for each isomer
    ! `nAct`     out    The number of active-state grains for each isomer
    ! `pa`       out    The steady state populations
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: nIsom
    integer, intent(in) :: nReac
    integer, intent(in) :: nProd
    integer, intent(in) :: nGrains
    real(8), intent(in) :: T
    real(8), intent(in) :: P
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), dimension(1:nIsom,1:nGrains,1:nGrains), intent(in) :: Mcoll
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nReac,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    integer, dimension(1:nIsom), intent(in) :: nRes
    integer, dimension(1:nIsom), intent(in) :: nAct
    real(8), dimension(1:nGrains,1:nIsom+nReac,1:nIsom), intent(out) :: pa
    character(len=128), intent(out) :: msg

    !! Accounting matrix
    integer, dimension(:,:), allocatable        ::  indices

    ! Active state grain matrix and RHS vectors
    real(8), dimension(:,:), allocatable    ::  L
    real(8), dimension(:,:), allocatable    ::  Z

    ! Variables for LAPACK
    integer, dimension(:), allocatable              ::  iPiv
    integer                                         ::  info

    ! Indices
    integer i, n, r, s

    ! Construct accounting matrix
    ! Row is grain number, column is well number, value is index into active-state matrix
    allocate( indices(1:nGrains, 1:nIsom ) )
    call accountingMatrix(nGrains, nIsom, nRes, indices)

    ! Create and zero active-state matrix and RHS vectors
    allocate( L(1:sum(nAct), 1:sum(nAct)), Z(1:sum(nAct), nIsom+nReac) )
    do i = 1, sum(nAct)
        do j = 1, sum(nAct)
            L(i,j) = 0.0
        end do
    end do
    do i = 1, sum(nAct)
        do j = 1, nIsom+nReac
            Z(i,j) = 0.0
        end do
    end do

    ! Collisional terms in active-state matrix and RHS vectors
    do i = 1, nIsom
        do r = nRes(i)+1, nGrains
            do s = nRes(i)+1, nGrains
                L(indices(r,i), indices(s,i)) = Mcoll(i,r,s)
            end do
            Z(indices(r,i), i) = sum(Mcoll(i,r,1:nRes(i)) * densStates(i,1:nRes(i)) * exp(-E(1:nRes(i)) / 8.314472 / T))
        end do
    end do

    ! Isomerization terms in active-state matrix and RHS vectors
    do i = 1, nIsom
        do j = 1, i-1
            do r = max(nRes(i), nRes(j))+1, nGrains
                L(indices(r,j), indices(r,i)) = Kij(j,i,r)
                L(indices(r,i), indices(r,i)) = L(indices(r,i), indices(r,i)) - Kij(j,i,r)
                L(indices(r,i), indices(r,j)) = Kij(i,j,r)
                L(indices(r,j), indices(r,j)) = L(indices(r,j), indices(r,j)) - Kij(i,j,r)
            end do
        end do
    end do

    ! Dissociation/association terms in active-state matrix and RHS vectors
    do i = 1, nIsom
        do n = 1, nReac+nProd
            do r = nRes(i)+1, nGrains
                L(indices(r,i), indices(r,i)) = L(indices(r,i), indices(r,i)) - Gnj(n,i,r)
            end do
        end do
        do n = 1, nReac
            do r = nRes(i)+1, nGrains
                Z(indices(r,i), n+nIsom) = Fim(i,n,r)
            end do
        end do
    end do

    Z = -Z

    ! Solve for pseudo-steady state populations of active state
    allocate( iPiv(1:sum(nAct)) )
    call DGESV(sum(nAct), nIsom+nReac, L, sum(nAct), iPiv, Z, sum(nAct), info)
    if (info /= 0) then
        msg = 'Active-state matrix is singular.'
        return
    end if
    deallocate( iPiv )

    ! Convert solution to pseudo-steady state populations
    do r = minval(nRes)+1, nGrains
        do n = 1, nIsom+nReac
            do i = 1, nIsom
                if (indices(r,i) > 0) pa(r,n,i) = Z(indices(r,i), n)
            end do
        end do
    end do

    ! Clean up
    deallocate( indices, L, Z )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine activeStateBanded(T, P, E, Mcoll, densStates, Kij, Fim, Gnj, dEdown, &
    nIsom, nReac, nProd, nGrains, nRes, nAct, pa, msg)
    ! Determine the pseudo-steady state populations for the active state
    ! grains using a banded matrix linear solve.
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `E`        in     A 1D array of energies in J/mol
    ! `Mcoll`    in     The collision matrix for each isomer
    ! `eqDist`   in     The normalized equilibrium distributions for each isomer
    ! `Kij`      in     The microcanonical isomerization rate coefficients in
    !                   s^-1
    ! `Fim`      in     The microcanonical association rates (rate coefficients
    !                   times bimolecular equilibrium distributions) in s^-1
    ! `Gnj`      in     The microcanonical dissociation rate coefficients in
    !                   s^-1
    ! `dEdown`   in     The average energy transferred in a deactivating
    !                   collision in J/mol
    ! `nIsom`    in     The number of isomers in the network
    ! `nProd`    in     The number of reactant/product channels in the network
    ! `nGrains`  in     The number of energy grains being used
    ! `nRes`     out    The number of reservoir grains for each isomer
    ! `nAct`     out    The number of active-state grains for each isomer
    ! `pa`       out    The steady state populations
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    ! Type definitions of parameters
    integer, intent(in) :: nIsom
    integer, intent(in) :: nReac
    integer, intent(in) :: nProd
    integer, intent(in) :: nGrains
    real(8), intent(in) :: T
    real(8), intent(in) :: P
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), dimension(1:nIsom,1:nGrains,1:nGrains), intent(in) :: Mcoll
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nReac,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    real(8), intent(in) :: dEdown
    integer, dimension(1:nIsom), intent(in) :: nRes
    integer, dimension(1:nIsom), intent(in) :: nAct
    real(8), dimension(1:nGrains,1:nIsom+nReac,1:nIsom), intent(out) :: pa
    character(len=128), intent(out) :: msg

    ! Accounting matrix
    integer, dimension(:,:), allocatable        ::  indices

    ! Active state grain matrix and RHS vectors
    real(8), dimension(:,:), allocatable        ::  L
    real(8), dimension(:,:), allocatable        ::  Z

    integer                                     ::  bandwidth, halfbandwidth
    real(8)                                     ::  dE

    ! Variables for LAPACK
    integer, dimension(:), allocatable          ::  iPiv
    integer                                     ::  info

    ! Indices
    integer i, j, n, r, s

    ! Construct accounting matrix
    ! Row is grain number, column is well number, value is index into active-state matrix
    allocate( indices(1:nGrains, 1:nIsom ) )
    indices = 0 * indices
    call accountingMatrix(nGrains, nIsom, nRes, indices)

    ! Determine bandwidth (at which transfer probabilities are so low that they can be truncated
    ! with negligible error)
    dE = E(2) - E(1)
    halfbandwidth = ceiling(16 * dEdown / dE) * nIsom
    bandwidth = 2 * halfbandwidth + 1

    ! Create and zero active-state matrix and RHS vectors
    allocate( L(3 * halfbandwidth + 1, 1:sum(nAct)), Z(1:sum(nAct), 1:nIsom+nReac) )
    do i = 1, 3 * halfbandwidth + 1
        do j = 1, sum(nAct)
            L(i,j) = 0.0
        end do
    end do
    do i = 1, sum(nAct)
        do j = 1, nIsom+nReac
            Z(i,j) = 0.0
        end do
    end do

    ! Collisional terms in active-state matrix and RHS vectors
    halfbandwidth = halfbandwidth / nIsom
    do i = 1, nIsom
        do s = nRes(i)+1, nGrains
            do r = max(nRes(i)+1, s - halfbandwidth), min(nGrains, s + halfbandwidth)
                L(bandwidth + indices(r,i) - indices(s,i), indices(s,i)) = Mcoll(i,r,s)
            end do
            Z(indices(s,i), i) = sum(Mcoll(i,s,1:nRes(i)) * densStates(i,1:nRes(i)) * exp(-E(1:nRes(i)) / 8.314472 / T))
        end do
    end do
    halfbandwidth = halfbandwidth * nIsom

    ! Isomerization terms in active-state matrix and RHS vectors
    do i = 1, nIsom
        do j = 1, i-1
            do r = max(nRes(i), nRes(j))+1, nGrains
                L(bandwidth + indices(r,i) - indices(r,j), indices(r,j)) = Kij(i,j,r)
                L(bandwidth, indices(r,j)) = L(bandwidth, indices(r,j)) - Kij(i,j,r)
                L(bandwidth + indices(r,j) - indices(r,i), indices(r,i)) = Kij(j,i,r)
                L(bandwidth, indices(r,i)) = L(bandwidth, indices(r,i)) - Kij(j,i,r)
            end do
        end do
    end do
    
    ! Dissociation/association terms in active-state matrix and RHS vectors
    do i = 1, nIsom
        do n = 1, nReac+nProd
            do r = nRes(i)+1, nGrains
                L(bandwidth, indices(r,i)) = L(bandwidth, indices(r,i)) - Gnj(n,i,r)
            end do
        end do
        do n = 1, nReac
            do r = nRes(i)+1, nGrains
                Z(indices(r,i), n+nIsom) = Fim(i,n,r)
            end do
        end do
    end do
    
    Z = -Z

    ! Solve for pseudo-steady state populations of active state
    allocate( iPiv(1:sum(nAct)) )
    call DGBSV(sum(nAct), halfbandwidth, halfbandwidth, nIsom+nReac, &
        L, 3 * halfbandwidth + 1, iPiv, Z, sum(nAct), info)
    if (info /= 0) then
        msg = 'Active-state matrix is singular.'
        return
    end if
    deallocate( iPiv )

    ! Convert solution to pseudo-steady state populations
    do r = minval(nRes)+1, nGrains
        do n = 1, nIsom+nReac
            do i = 1, nIsom
                if (indices(r,i) > 0) pa(r,n,i) = Z(indices(r,i), n)
            end do
        end do
    end do

    ! Clean up
    deallocate( indices, L, Z )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine accountingMatrix(nGrains, nIsom, nRes, indices)
    ! Determine a mapping of isomer and grain indices to matrix element indices.
    ! This inner grouping is isomer/product index, while the outer grouping is
    ! energy grain index; this results in a much more tightly banded matrix than
    ! the reverse. The output is a matrix where row is grain number, column is
    ! isomer/product number, and the value is the index into the matrix.

    integer, intent(in)                     ::  nGrains
    integer, intent(in)                     ::  nIsom
    integer, dimension(1:nIsom), intent(in) ::  nRes
    integer, dimension(1:nGrains,1:nIsom), intent(out)  :: indices

    integer i, r, row

    ! Construct accounting matrix
    row = 1
    do r = 1, nGrains
        do i = 1, nIsom
            if (r > nRes(i)) then
                indices(r,i) = row
                row = row + 1
            else
                indices(r,i) = 0
            end if
        end do
    end do

end subroutine

