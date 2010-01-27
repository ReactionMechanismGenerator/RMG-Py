!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	RMG - Reaction Mechanism Generator
!
!	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
!	RMG Team (rmg_dev@mit.edu)
!
! 	Copyright (c) 2009 by Josh Allen (chemejosh@gmail.com)
!
! 	Permission is hereby granted, free of charge, to any person obtaining a
! 	copy of this software and associated documentation files (the 'Software'),
! 	to deal in the Software without restriction, including without limitation
! 	the rights to use, copy, modify, merge, publish, distribute, sublicense,
! 	and/or sell copies of the Software, and to permit persons to whom the
! 	Software is furnished to do so, subject to the following conditions:
!
! 	The above copyright notice and this permission notice shall be included in
! 	all copies or substantial portions of the Software.
!
! 	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! 	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! 	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! 	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! 	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! 	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! 	DEALINGS IN THE SOFTWARE.
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
	do k = 0, 16 - v
		term = (0.5 * z)**(v+2*k) / (factorial(k) * factorial(k+v))
		res = res + term
		if (abs(term / res) < 1.0e-8) return
	end do
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

	A = 1.0 + x
	B = 1.0 - x
	do n = 0, 100
		A0 = A
		B0 = B
		A = (A0 + B0) / 2
		B = sqrt(A0 * B0)
		if (abs(A - B) < 1.0e-8) then
			cellipk = 3.141592654 / 2.0 / A
			return
		end if
	end do

	cellipk = 3.141592654 / 2.0 / A

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
