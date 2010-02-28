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
    real(8), dimension(:,:), allocatable :: Mcoll, Mrxn, M, X, Xinv
    real(8), dimension(:), allocatable :: S, Sinv, V0, V

    ! BLAS and LAPACK temporary variables
    real(8), dimension(:), allocatable :: work, iPiv
    integer info

    ! Indices
    integer i, j, r, index, count

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
    !do i = 1, nRows
    !    do j = 1, i
    !        if (M(i,j) /= 0.0) then
    !            if (abs(M(i,j) - M(j,i)) / M(i,j) > 0.01) then
    !                write (*,*) i, j, M(i,j), M(j,i)
    !            end if
    !        end if
    !    end do
    !end do

    ! Zero phenomenological rate constant matrix
    do i = 1, nIsom + nProd
        do j = 1, nIsom + nProd
            K(i,j) = 0.0
        end do
    end do

    ! Calculate the eigenvalues and eigenvectors
    allocate( V0(1:nRows), work(1:6*nRows) )
    call DSYEV('V', 'U', nRows, M, nRows, V0, work, 6*nRows, info)
    if (info > 0) then
        msg = 'DSYEV eigenvalue algorithm failed to converge.'
        return
    elseif (info < 0) then
        msg = 'Illegal argument passed to DSYEV.'
        return
    end if

    ! Check for proper separation of chemical timescales
    if (abs(V0(nRows-nIsom-nProd) / V0(nRows-nIsom-nProd+1)) < 2.0) then
        write (*,fmt='(A)') 'Error: Chemical eigenvalues are not distinct from internal energy eigenvalues.'
        write (*,fmt='(A,ES12.4E2,A,ES12.4E2,A,F7.4)') &
            'last IERE =', V0(nRows-nIsom-nProd), &
            ', first CSE =',V0(nRows-nIsom-nProd+1), &
            ', ratio =',abs(V0(nRows-nIsom-nProd) / V0(nRows-nIsom-nProd+1))
        msg = 'Chemical eigenvalues not distinct from internal energy eigenvalues.'
        return
    end if
    
    ! Check for zero eigenvalue
    if (abs(V0(nRows) / V0(nRows-1)) > 0.001) then
        write (*,fmt='(A)') 'Error: Zero eigenvalue not found.'
        write (*,fmt='(A,ES12.4E2,A,ES12.4E2,A,F9.6)') &
            'zero CSE =', V0(nRows), &
            ', last CSE =',V0(nRows-1), &
            ', ratio =',abs(V0(nRows) / V0(nRows-1))
        msg = 'Zero eigenvalue not found.'
        return
    end if

    ! Calculate unsymmetrized eigenvectors
    do i = 1, nRows
        do j = 1, nRows
            M(i,j) = Sinv(i) * M(i,j)
        end do
    end do
    ! Also normalize the eigenvectors
    do j = 1, nRows
        M(:,j) = M(:,j) / sum(M(:,j))
    end do

    ! Extract the chemically-significant eigenvalues and eigenvectors
    ! These are the last nIsom+nProd columns of M
    allocate( X(1:nIsom+nProd, 1:nIsom+nProd), V(1:nIsom+nProd), Xinv(1:nIsom+nProd, 1:nIsom+nProd))
    do i = 1, nIsom + nProd
        do j = 1, nIsom + nProd
            X(i,j) = 0.0
        end do
    end do
    do j = 1, nIsom + nProd
        count = nRows - (nIsom + nProd) + j
        V(j) = V0(count)
        do r = 1, nGrains
            do i = 1, nIsom
                index = indices(r,i)
                if (index > 0) then
                    X(i,j) = X(i,j) + M(index,count)
                end if
            end do
        end do
        do i = 1, nProd
            index = nRows - nProd + i
            X(i+nIsom,j) = M(index,count)
        end do
    end do
    
    ! Check: The eigenvector corresponding to the zero eigenvalue should
    ! be identical to eqRatios
    !do j = 1, nIsom+nProd
    !    if (abs(X(j,nIsom+nProd) - eqRatios(j)) / eqRatios(j) > 0.001) then
    !        return
    !    end if
    !end do

    ! Invert the chemical eigenvector matrix
    Xinv = X
    allocate(iPiv(1:nIsom+nProd))
    call DGETRF(nIsom+nProd, nIsom+nProd, Xinv, nIsom+nProd, ipiv, info)
    if (info > 0) then
        msg = 'Chemical eigenvector matrix is singular.'
        return
    elseif (info < 0) then
        msg = 'Illegal argument passed to DGETRF.'
        return
    end if
    call DGETRI(nIsom+nProd, Xinv, nIsom+nProd, ipiv, work, 6*nRows, info)
    if (info > 0) then
        msg = 'Chemical eigenvector matrix is singular.'
        return
    elseif (info < 0) then
        msg = 'Illegal argument passed to DGETRI.'
        return
    end if

    ! Construct the rate constant matrix
    do i = 1, nIsom+nProd
        do j = 1, nIsom+nProd
            K(i,j) = 0.0
            do r = 1, nIsom+nProd
                K(i,j) = K(i,j) + X(i,r) * V(r) * Xinv(r,j)
            end do
        end do
    end do

    ! Clean up
    deallocate( Mcoll, Mrxn, M, S, Sinv, V0, work, iPiv, X, V, Xinv )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
