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

subroutine estimateRateCoefficients_MSC(T, P, E, collFreq, densStates, Eres, &
    Kij, Fim, Gnj, nIsom, nReac, nProd, nGrains, K, msg)
    ! Estimate the phenomenological rate coefficients using the (modified) strong
    ! collision method. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `T`        in     The temperature to evaluate k(T,P) at in K
    ! `P`        in     The pressure to evaluate k(T,P) at in Pa
    ! `E`        in     A 1D array of energies in J/mol
    ! `collFreq` in     The (modified) collision frequencies for each isomer in
    !                   Hz
    ! `eqDist`   in     The normalized equilibrium distributions for each isomer
    ! `Eres`     in     The active-state energy cutoffs for each isomer in J/mol
    ! `Kij`      in     The microcanonical isomerization rate coefficients in
    !                   s^-1
    ! `Fim`      in     The microcanonical association rates (rate coefficients
    !                   times bimolecular equilibrium distributions) in s^-1
    ! `Gnj`      in     The microcanonical dissociation rate coefficients in
    !                   s^-1
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
    real(8), dimension(1:nIsom), intent(in) :: collFreq
    real(8), dimension(1:nIsom+nReac+nProd), intent(in) :: Eres
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nReac,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    real(8), dimension(1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K
    character(len=128), intent(out) :: msg

    ! Steady-state populations
    real(8), dimension(:,:,:), allocatable      ::  pa
    
    ! Steady-state matrix and vector
    real(8), dimension(:,:), allocatable        ::  A
    real(8), dimension(:,:), allocatable        ::  b
    ! Indices i and j represent sums over unimolecular wells
    integer                                     ::  i, j
    ! Indices m and n represent sums over bimolecular sources/sinks
    integer                                     ::  n
    ! Indices r and s represent sums over energy grains
    integer                                     ::  r
    ! Variables for BLAS and LAPACK
    integer, dimension(:), allocatable          ::  iPiv
    integer                                     ::  info
    ! The isomer or channel of current interest
    integer                                     ::  src
    ! The active-state energy grain of the current isomer
    integer                                     ::  start
    ! A temporary used to store double-precision values
    real(8)                                     ::  val

    ! Allocate matrices and vectors used in this subroutine
    allocate( A(1:nIsom, 1:nIsom) )
    allocate( b(1:nIsom, 1:nIsom+nReac) )
    allocate( iPiv(1:nIsom) )
    allocate( pa(1:nGrains, 1:nIsom, 1:nIsom+nReac) )

    ! Zero the phenomenological rate coefficient matrix
    do i = 1, nIsom + nReac + nProd
        do j = 1, nIsom + nReac + nProd
            K(i,j) = 0.0
        end do
    end do

    ! Set msg to successful; will be changed later if an error occurs
    msg = ""

    ! Determine the starting grain for the calculation based on the
    ! active-state cutoff energy
    start = ceiling((minval(Eres) - minval(E)) / (E(2) - E(1))) + 1
    if (start < 1 .or. start > nGrains) then
        msg = 'Unable to determine starting grain; check active-state energies.'
        return
    end if
        
    ! Zero the steady-state population vector for grains below the cutoff
    do r = 1, start-1
        do i = 1, nIsom
            do n = 1, nIsom+nReac
                pa(r,i,n) = 0.0
            end do
        end do
    end do

    ! Iterate over the grains, calculating the PSSA concentrations
    do r = start, nGrains

        ! Zero A matrix and b vector
        do i = 1, nIsom
            do j = 1, nIsom
                A(i,j) = 0.0
            end do
            do n = 1, nIsom+nReac
                b(i,n) = 0.0
            end do
        end do

        ! Collisional deactivation
        do i = 1, nIsom
            A(i,i) = A(i,i) - collFreq(i)
        end do

        ! Isomerization reactions
        do i = 1, nIsom
            do j = 1, i - 1
                A(i,j) = Kij(i,j,r)
                A(j,j) = A(j,j) - Kij(i,j,r)
                A(j,i) = Kij(j,i,r)
                A(i,i) = A(i,i) - Kij(j,i,r)
            end do
        end do

        ! Dissociation reactions
        do n = 1, nReac+nProd
            do j = 1, nIsom
                A(j,j) = A(j,j) - Gnj(n,j,r)
            end do
        end do

        ! Populate RHS vectors, one per isomer and reactant
        do n = 1, nIsom + nReac
            if (n <= nIsom) then
                ! Thermal activation via collisions
                b(n,n) = collFreq(n) * densStates(n, r) * exp(-E(r) / 8.314472 / T)
            else
                ! Chemical activation via association reaction
                do j = 1, nIsom
                    b(j,n) = Fim(j,n-nIsom,r)
                end do
            end if
        end do

        ! Solve for steady-state population
        call DGESV( nIsom, nIsom+nReac, A, nIsom, iPiv, b, nIsom, info )
        if (info > 0) then
            msg = "A singular matrix was encountered."
            return
        end if
        pa(r,:,:) = -b

        ! Check that our populations are all positive
        do i = 1, nIsom
            do n = 1, nIsom+nReac
                if (pa(r,i,n) < 0.0) then
                    msg = "A negative steady-state concentration was encountered."
                    return
                end if
            end do
        end do

    end do


    do src = 1, nIsom + nReac

        ! Calculate stabilization rates (i.e.) R + R' --> Ai or M --> Ai
        do i = 1, nIsom
            if (i /= src) then
                val = collFreq(i) * sum(pa(:,i,src))
                K(i,src) = K(i,src) + val
                K(src,src) = K(src,src) - val
            end if
        end do

        ! Calculate dissociation rates (i.e.) R + R' --> Bn + Cn or M --> Bn + Cn
        do n = 1, nReac+nProd
            do j = 1, nIsom
                if (n+nIsom /= src) then
                    val = sum(Gnj(n,j,:) * pa(:,j,src))
                    K(n+nIsom,src) = K(n+nIsom,src) + val
                    K(src,src) = K(src,src) - val
                end if
            end do
        end do

    end do

    ! Clean up
    deallocate( A, b, iPiv, pa )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
