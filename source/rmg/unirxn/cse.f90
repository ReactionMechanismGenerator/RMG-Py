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
    eqRatios, Kij, Fim, Gnj, nIsom, nReac, nProd, nGrains, K, msg)
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
    ! `nReac`    in     The number of reactant channels in the network (both A + B <=> C)
    ! `nProd`    in     The number of product channels in the network (A -> B + C only)
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
    integer, intent(in) :: nReac
    integer, intent(in) :: nProd
    integer, intent(in) :: nGrains
    real(8), dimension(1:nGrains), intent(in) :: E
    real(8), dimension(1:nIsom,1:nGrains,1:nGrains), intent(in) :: Mcoll0
    real(8), dimension(1:nIsom+nReac+nProd), intent(in) :: E0
    real(8), dimension(1:nIsom,1:nGrains), intent(in) :: densStates
    real(8), dimension(1:nIsom+nReac+nProd), intent(in) :: eqRatios
    real(8), dimension(1:nIsom,1:nIsom,1:nGrains), intent(in) :: Kij
    real(8), dimension(1:nIsom,1:nReac,1:nGrains), intent(in) :: Fim
    real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains), intent(in) :: Gnj
    real(8), dimension(1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K
    character(len=128), intent(out) :: msg

    ! The size of the master equation matrix
    integer nRows

    ! The number of chemically-significant eigenvalues
    integer nCSE

    ! The indices of the master equation matrix
    integer, dimension(1:nGrains,1:nIsom) :: indices

    ! Intermediate matrices
    real(8), dimension(:,:), allocatable :: Mcoll, Mrxn, M, X, Xinv, eqDist, dXij
    real(8), dimension(:), allocatable :: S, Sinv, V0, V, C

    ! BLAS and LAPACK temporary variables
    real(8), dimension(:), allocatable :: work
    integer info

    ! Indices
    integer i, j, n, r, h, index, count

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
    ! Reactant wells are appended to matrix, one row each
    ! Product wells are NOT added because we are neglecting reassociation
    nRows = nRows + nReac

    ! Allocate intermediate matrices
    allocate( Mcoll(1:nRows, 1:nRows), Mrxn(1:nRows, 1:nRows), M(1:nRows, 1:nRows) )
    allocate( S(1:nRows), Sinv(1:nRows) )

	! Generate full unsymmetrized master equation matrix (dimension nRows x nRows)
    call fullMEMatrix(E, E0, Mcoll0, Kij, Gnj, Fim, indices, &
        nRows, nGrains, nIsom, nReac, nProd, Mcoll, Mrxn, msg)

    ! Generate symmetrization matrix and its inverse
    do r = 1, nGrains
        do i = 1, nIsom
            index = indices(r,i)
            if (index > 0) then
                S(index) = sqrt(densStates(i,r) * exp(-E(r) / 8.314472 / T) * eqRatios(i))
                Sinv(index) = 1.0 / S(index)
            end if
        end do
    end do
    do i = 1, nReac
        index = nRows - nReac + i
        S(index) = sqrt(1.0 * eqRatios(nIsom+i))
        Sinv(index) = 1.0 / S(index)
    end do

    ! Symmetrize master equation matrix: M = S * Msymm * Sinv
    ! Since S and Sinv are diagonal we can do this very efficiently
    do i = 1, nRows
        do j = 1, nRows
            M(i,j) = Sinv(i) * (Mcoll(i,j) + Mrxn(i,j)) * S(j)
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
    do i = 1, nIsom + nReac + nProd
        do j = 1, nIsom + nReac + nProd
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

!    ! Calculate the eigenvalues and eigenvectors
!    ! Only calculate the ones we need, which are a small subset of the total
!    ! This doesn't seem to work properly; it only gives one eigenvector
!    allocate( V0(1:nRows), X0(1:nRows, 1:nRows) )
!    allocate( work(1:6*nRows), iwork(1:8*nRows), ifail(1:5*nRows) )
!    call DSYEVX('V', 'I', 'U', nRows, M, nRows, 0.0, 0.0, &
!        nRows-nIsom-nProd, nRows, 1.0e-8, nCSE, V0, X0, nRows, &
!        work, 8*nRows, iwork, ifail, info)
!    if (info > 0) then
!        msg = 'DSYEVX eigenvalue algorithm failed to converge.'
!        return
!    elseif (info < 0) then
!        msg = 'Illegal argument passed to DSYEVX.'
!        return
!    end if

    ! Count the number of chemically distinct eigenvalues
    ! This will be ideally be nIsom+nProd, but might be lower if one or more
    ! chemical eigenmodes is blended with the internal energy eigenmodes
    nCSE = 1
    do i = 1, nIsom+nReac-1
        if (abs(V0(nRows-nIsom-nReac) / V0(nRows-i)) > 5.0) then
            nCSE = nCSE + 1
        end if
    end do

    ! Check for proper separation of chemical timescales
    if (nCSE /= nIsom + nReac) then
        write (*,fmt='(A)') 'Error: Chemical eigenvalues are not distinct from internal energy eigenvalues.'
        write (*,fmt='(A,ES12.4E2,A,ES12.4E2,A,F7.4)') &
            'last IERE =', V0(nRows-nIsom-nReac), &
            ', first CSE =',V0(nRows-nIsom-nReac+1), &
            ', ratio =',abs(V0(nRows-nIsom-nReac) / V0(nRows-nIsom-nReac+1))
        msg = 'Chemical eigenvalues not distinct from internal energy eigenvalues.'
        return
    end if

    ! Check for zero eigenvalue (only if there are no product channels (sinks))
    if (nProd == 0) then
        if (abs(V0(nRows) / V0(nRows-1)) > 0.001) then
            write (*,fmt='(A)') 'Error: Zero eigenvalue not found.'
            write (*,fmt='(A,ES12.4E2,A,ES12.4E2,A,F9.6)') &
                'zero CSE =', V0(nRows), &
                ', last CSE =',V0(nRows-1), &
                ', ratio =',abs(V0(nRows) / V0(nRows-1))
            msg = 'Zero eigenvalue not found.'
            return
        end if
        ! Force largest eigenvalue to be zero (might be small numerical artifact)
        !V0(nRows) = 0.0
    end if

    allocate( X(1:nRows, 1:nCSE), V(1:nCSE), Xinv(1:nCSE, 1:nRows) )
    do j = 1, nCSE
        ! Find index in full matrix
        count = nRows - (nIsom + nReac) + j
        ! Eigenvalue matrix V
        V(j) = V0(count)
        ! Eigenvector matrix and its inverse/transpose
        do r = 1, nGrains
            do i = 1, nIsom
                index = indices(r,i)
                if (index > 0) then
                    ! Unsymmetrize as we go
                    X(index,j) = S(index) * M(index,count)
                    Xinv(j,index) = M(index,count) * Sinv(index)
                end if
            end do
        end do
        do i = 1, nReac
            index = nRows - nReac + i
            ! Unsymmetrize as we go
            X(index,j) = S(index) * M(index,count)
            Xinv(j,index) = M(index,count) * Sinv(index)
        end do
    end do

    ! Generate set of initial condition vectors, one per isomer and reactant
    ! Each initial condition vector corresponds to a case where a Boltzmann
    ! distribution of only one isomer or reactant is present
    ! The Boltzmann distribution occurs as a result of the fast relaxation of
    ! internal energy modes, which is separable from the chemical modes
    allocate( eqDist(1:nRows, 1:nIsom+nReac) )
    do j = 1, nIsom + nReac
        do r = 1, nRows
            eqDist(r,j) = 0.0
        end do
    end do
    do r = 1, nGrains
        do i = 1, nIsom
            index = indices(r,i)
            if (index > 0) then
                eqDist(index,i) = densStates(i,r) * exp(-E(r) / 8.314472 / T)
            end if
        end do
    end do
    do i = 1, nReac
        index = nRows - nReac + i
        eqDist(index,i+nIsom) = 1.0
    end do

    ! Calculate the phenomenological rate constants
    ! This version follows the notation of Miller and Klippenstein
    allocate( dXij(1:nIsom+nReac+nProd, 1:nCSE), C(1:nRows) )
    ! Iterate over isomers and reactants to determine k(T,P) with each as the
    ! reactant
    do n = 1, nIsom+nReac
        do j = 1, nCSE

            ! dXij contains the change in isomer/reactant/product i as a result 
            ! of the jth eigenmode, using isomer/reactant n as the starting 
            ! point
            do i = 1, nIsom+nReac+nProd
                dXij(i,j) = 0.0
            end do

            do r = 1, nRows
                C(r) = X(r,j) * sum(Xinv(j,:) * eqDist(:,n))
            end do

            do r = 1, nGrains
                do i = 1, nIsom
                    index = indices(r,i)
                    if (index > 0) then
                        ! Isomers
                        dXij(i,j) = dXij(i,j) + C(index)
                        ! Products
                        do h = nReac+1, nReac+nProd
                            dXij(nIsom+h,j) = dXij(nIsom+h,j) + C(index) * Gnj(h,i,r) / V(j)
                        end do
                    end if
                end do
            end do
            do i = 1, nReac
                index = nRows - nReac + i
                ! Reactants
                dXij(nIsom+i,j) = dXij(nIsom+i,j) + C(index)
            end do

            ! Convert the dXij information for eigenmode j into phenomenological
            ! rate coefficients
            do i = 1, nIsom+nReac+nProd
                K(i,n) = K(i,n) + V(j) * dXij(i,j)
            end do

        end do

    end do

    ! Clean up
    deallocate( Mcoll, Mrxn, M, S, Sinv, V0, work, eqDist, X, V, Xinv, dXij, C )

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
