**************************
Frequently Asked Questions
**************************

**Are there other software packages for investigating pressure-dependent
reaction networks?**

Yes. The following is an illustrative list of such packages:

=============== =============== =============== ================================
Name            Method(s)       Language        Author(s)                       
=============== =============== =============== ================================
`MultiWell`_    stochastic      Fortran         J.\  R. Barker *et al*
`UNIMOL`_       CSE             Fortran         R.\  G. Gilbert, S. C. Smith      
`ChemRate`_     CSE             C++ [#f1]_      V.\  Mokrushin, W. Tsang          
`Variflex`_     CSE             Fortran         S.\  J. Klippenstein *et al*
`MESMER`_       CSE (+ RS)      C++             S.\  H. Robertson *et al*
CHEMDIS [#f2]_  MSC             Fortran         A.\  Y. Chang, J. W. Bozzelli, A. M. Dean
=============== =============== =============== ================================

(MSC = modified strong collision, RS = reservoir state, CSE = chemically-significant eigenvalues)

Many of the above packages also provide additional functionality beyond the
approximate solving of the master equation. For example, Variflex can be used
for variational transition state theory calculations, while ChemRate provides a
(Windows) graphical user interface for exploring a database of experimental
data and physical quantities.

.. [#f1] Uses MFC for Windows graphical user interface

.. [#f2] No longer distributed

.. _MultiWell: http://aoss-research.engin.umich.edu/multiwell/
.. _UNIMOL: http://www.ccl.net/cca/software/SOURCES/FORTRAN/unimol/index.shtml
.. _ChemRate: http://www.mokrushin.com/ChemRate/chemrate.html
.. _Variflex: http://ftp.tcg.anl.gov/pub/variflex/Summary.vrfx
.. _MESMER: http://sourceforge.net/projects/mesmer


