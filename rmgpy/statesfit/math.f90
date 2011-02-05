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

recursive function factorial(x) result(res)
    implicit none
    integer, intent(in) :: x
    integer :: res
    if (x <= 1) then
        res = 1
    else
        res = x * factorial(x - 1)
    end if
end function

function besseli(v, z) result(res)
    implicit none
    integer, intent(in) :: v
    real(8), intent(in) :: z
    real(8) :: res
    integer k, factorial
    real(8) term
    res = 0.0
    do k = 0, 25 - v
        term = (0.25 * z * z)**k / factorial(k) / factorial(k + v)
        res = res + term
        if (abs(term / res) < 1.0e-8) then
            res = res * (0.5 * z)**v
            return
        end if
    end do
    ! If we're here, then we didn't reach that tolerance
    ! Let's try a looser tolerance
    res = 0.0
    do k = 0, 25 - v
        term = (0.25 * z * z)**k / factorial(k) / factorial(k + v)
        res = res + term
        if (abs(term / res) < 1.0e-6) then
            res = res * (0.5 * z)**v
            return
        end if
    end do
    ! If we're here, then we didn't reach that tolerance either
    ! Let's try one more looser tolerance
    res = 0.0
    do k = 0, 25 - v
        term = (0.25 * z * z)**k / factorial(k) / factorial(k + v)
        res = res + term
        if (abs(term / res) < 1.0e-4) then
            res = res * (0.5 * z)**v
            return
        end if
    end do
    ! If we're here, then something's really gone bad
    ! Let's give up
    write (*,*) 'Warning: Unable to determine value of modified Bessel function for v =', v, 'z =', z
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function besseli0(x)
    implicit none

    real(8) x, besseli0

    real(8) Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX

    data P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
    0.2659732D0,0.360768D-1,0.45813D-2/

    data Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
    0.2635537D-1,-0.1647633D-1,0.392377D-2/

    if (abs(x) < 3.75) then
        Y = x * x / 3.75 / 3.75
        besseli0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
    else
        AX = abs(x)
        Y = 3.75 / AX
        BX = exp(AX)/sqrt(AX)
        AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        besseli0 = AX*BX
    end if

end function

function besseli1(x)
    implicit none

    real(8) x, besseli1

    real(8) Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX

    data P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/

    data Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
    -0.2895312D-1,0.1787654D-1,-0.420059D-2/

    if (abs(x) < 3.75) then
        Y = x * x / 3.75 / 3.75
        besseli1 = x*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    else
        AX = abs(x)
        Y = 3.75 / AX
        BX = exp(AX)/sqrt(AX)
        AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        besseli1 = AX*BX
    end if

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function besselratioperron(x, n) result(res)
    ! Recursively evaluates Perron's continued fraction for the modified
    ! Bessel function ratio I_1(x) / I_0(x) to 25 terms. A relevant citation is
    !
    ! W. Gautschi and J. Slavic. "On the Computation of Modified Bessel
    ! Function Ratios." Mathematics of Computation 32 (143), p. 865-875 (1978).
    
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: n
    real(8) :: res

    if (n >= 25) then
        res = 0
    else
        res = (2 * n + 1) * x / (n + 2 + 2 * x - besselratioperron(x, n+1))
    end if
    
end function

function besselratio(x) result(res)
    ! Computes the ratio of modified bessel functions I_1(x) / I_0(x) using
    ! Perron's continued fraction.

    implicit none
    real(8), intent(in) :: x
    real(8) :: res
    
    real(8) :: besselratioperron

    res = x / (2 + x - besselratioperron(x, 1))

end function

function cellipk(x)
    implicit none
    real(8), intent(in) :: x
    real(8) :: cellipk

    real(8) A0, B0, A, B
    integer n

    if (x < 0 .or. x > 1) then
        cellipk = 0.0
        return
    end if

    A = 1.0
    B = sqrt(1.0 - x)
    do n = 0, 100
        A0 = A
        B0 = B
        A = (A0 + B0) / 2
        B = sqrt(A0 * B0)
        if (abs(A - B) < 1.0e-12) then
            cellipk = 3.141592654 / 2.0 / A
            return
        end if
    end do

    cellipk = 3.141592654 / 2.0 / A

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
