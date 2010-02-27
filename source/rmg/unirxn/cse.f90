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

subroutine estimateRateCoefficients_CSE(T, P, E, Mcoll0, E0, densStates, &
    eqRatios, Kij, Fim, Gnj, nIsom, nProd, nGrains, K, msg)
    ! Estimate the phenomenological rate coefficients using the chemically-
    ! significant eigenvalues method. The parameters are:
    !
    ! ========== ====== ========================================================
    ! Parameter  Intent Description
    ! ========== ====== ========================================================
    ! `T`        in     The temperature to evaluate k(T,P) at in K
    ! `P`        in     The pressure to evaluate k(T,P) at in Pa
    ! `E`        in     A 1D array of energies in J/mol
    ! `Mcoll0`   in     The collision matrix for each isomer
    ! `densStates` in     The density of states for each isomer
    ! `eqRatios` in     The ratio of each isomer/product channel at equilibrium
    ! `Kij`      in     The microcanonical isomerization rate coefficients in
    !                   s^-1
    ! `Fim`      in     The microcanonical association rates (rate coefficients
    !                   times bimolecular equilibrium distributions) in s^-1
    ! `Gnj`      in     The microcanonical dissociation rate coefficients in
    !                   s^-1
    ! `nIsom`    in     The number of isomers in the network
    ! `nProd`    in     The number of reactant/product channels in the network
    ! `nGrains`  in     The number of energy grains being used
    ! `K`        out    The matrix of phenomenological rate coefficients k(T,P)
    ! `msg`      out    If the subroutine was unsuccessful, this string will
    !                   contain a brief message describing the error; the
    !                   string will be empty if the subroutine was successful
    ! ========== ====== ========================================================

    implicit none

    ! Type definitions of parameters
    real(8), intent(in) :: T
    real(8), intent(in) :: P
    integer, intent(in) :: nIsom
    integer, intent(in) :: nProd
    integer, intent(in) :: nGrains
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), dimension(1:nIsom,1:nGrains,1:nGrains), intent(in) :: Mcoll0
    real(8), dimension(1:nIsom+nProd), intent(in) :: E0
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    real(8), dimension(1:nIsom+nProd), intent(in) :: eqRatios
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nProd,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    real(8), dimension(1:nIsom+nProd,1:nIsom+nProd), intent(out) :: K
    character(len=128), intent(out) :: msg

    ! The size of the master equation matrix
    integer nRows
    
    ! The indices of the master equation matrix
    integer, dimension(1:nGrains,1:nIsom) :: indices
    
    ! Intermediate matrices
    real(8), dimension(:,:), allocatable :: Mcoll, Mrxn, M, p0, dXij
    real(8), dimension(:), allocatable :: S, Sinv, V, Cj
    
    ! BLAS and LAPACK temporary variables
    real(8), dimension(:), allocatable :: work
    integer info

    ! Indices
    integer i, j, r, index, start

    ! Construct accounting matrix
	! Row is grain number, column is well number, value is index into matrix
    nRows = 0
    do r = 1, nGrains
        do i = 1, nIsom
            if (densStates(i,r) > 0) then
                nRows = nRows + 1
                indices(r,i) = nRows
            else
                indices(r,i) = 0
            end if
        end do
    end do
    ! Product wells are appended to matrix, one row each
    nRows = nRows + nProd

    ! Allocate intermediate matrices
    allocate( Mcoll(1:nRows, 1:nRows), Mrxn(1:nRows, 1:nRows), M(1:nRows, 1:nRows) )
    allocate( S(1:nRows), Sinv(1:nRows) )

	! Generate full unsymmetrized master equation matrix (dimension nRows x nRows)
    call fullMEMatrix(E, E0, Mcoll0, Kij, Gnj, Fim, indices, &
        nRows, nGrains, nIsom, nProd, Mcoll, Mrxn, msg)
    
	! Generate symmetrization matrix and its inverse
    do r = 1, nGrains
        do i = 1, nIsom
            index = indices(r,i)
            if (index > 0) then
                Sinv(index) = sqrt(densStates(i,r) * exp(-E(r) / 8.314472 / T) * eqRatios(i))
                S(index) = 1.0 / Sinv(index)
            end if
        end do
    end do
    do i = 1, nProd
        index = nRows - nProd + i
        Sinv(index) = sqrt(1.0 * eqRatios(nIsom+i))
        S(index) = 1.0 / Sinv(index)
    end do

    ! Symmetrize master equation matrix: Msymm = S * M * Sinv
    ! Since S and Sinv are diagonal we can do this very efficiently
    do i = 1, nRows
        do j = 1, nRows
            M(i,j) = S(i) * (Mcoll(i,j) + Mrxn(i,j)) * Sinv(j)
        end do
    end do

	! DEBUG: Check that the matrix has been properly symmetrized
    do i = 1, nRows
        do j = 1, i
            if (M(i,j) /= 0.0) then
                if (abs(M(i,j) - M(j,i)) / M(i,j) > 0.001) then
                    write (*,*) i, j, M(i,j), M(j,i)
                end if
            end if
        end do
    end do

    ! Zero output rate constant matrix
    do i = 1, nIsom+nProd
        do j = 1, nIsom+nProd
            K(i,j) = 0.0
        end do
    end do

    ! Calculate the eigenvalues and eigenvectors
    allocate( V(1:nRows), work(1:6*nRows) )
    call DSYEV('V', 'U', nRows, M, nRows, V, work, 6*nRows, info)
    if (info > 0) then
        msg = 'DSYEV eigenvalue algorithm failed to converge.'
        return
    elseif (info < 0) then
        msg = 'Illegal argument passed to DSYEV.'
        return
    end if

    ! Assemble post-internal-energy-relaxation initial condition vector
    ! using each isomer/product channel as a starting point
    ! For isomers this is a normalized Boltzmann distribution
    allocate(p0(1:nRows,1:nIsom+nProd))
    do r = 1, nGrains
        do i = 1, nIsom
            index = indices(r,i)
            if (index > 0) then
                do j = 1, nIsom+nProd
                    p0(index,j) = 0.0
                end do
                p0(index,i) = densStates(i,r) * exp(-E(r) / 8.314472 / T)
            end if
        end do
    end do
    ! For products this is simply unity
    do i = 1, nProd
        index = nRows - nProd + i
        do j = 1, nIsom+nProd
            p0(index,j) = 0.0
        end do
        p0(index,i+nIsom) = 1.0
    end do
    
    ! Iterate over each isomer/product channel as the starting point
    allocate( Cj(1:nIsom+nProd), dXij(1:nIsom+nProd,1:nIsom+nProd) )
    do start = 1, nIsom+nProd
        ! Iterate over the chemically-significant eigenvalue-eigenvector pairs
        do j = 1, nIsom+nProd
            index = nRows - (nIsom+nProd) + j
            ! Calculate initial constant for each chemically-significant 
            ! eigenvector j (using isomer `start` as the starting position)
            Cj(j) = sum(M(:,index) * p0(:,start))
            ! Calculate dXij for isomer/product i and chemically-significant
            ! eigenvector j (using isomer `start` as the starting position)
            do i = 1, nIsom+nProd
                dXij(i,j) = - Cj(j) * sum(M(:,index) * p0(:,i))
            end do
            ! Add to rate constants
            do i = 1, nIsom+nProd
                K(i,start) = K(i,start) - V(index) * dXij(i,j)
            end do
        end do
    end do

    ! Clean up
    deallocate( Mcoll, Mrxn, M, S, Sinv, V, work, p0, Cj, dXij )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
