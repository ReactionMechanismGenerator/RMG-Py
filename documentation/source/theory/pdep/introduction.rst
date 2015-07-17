************
Introduction
************

Unimolecular Reactions
======================

Unimolecular reactions are those that involve a single reactant or product
molecule, the union of isomerization and dissociation/association reactions:
 
.. math:: 
    
    \ce{A <=> B} & \hspace{40pt} \text{isomerization} \\
    \ce{A <=> B + C} & \hspace{40pt} \text{dissociation/association}
    
Gas-phase chemical reactions occur as the result of bimolecular collisions
between two reactant molecules. This presents a problem when there is only one
participating reactant molecule! The conclusion is that the above reactions
cannot be elementary as written; another step must be involved.

For a unimolecular reaction to proceed, the reactant molecule :math:`\ce{A}`
must first be excited to an energy that exceeds the barrier for reaction. A
molecule that is sufficiently excited to react is called an *activated 
species* and often labeled with an asterisk :math:`\ce{A}^\ast`. If we
replace the stable species with the activated species in the reactions above,
the reactions become elementary again:

.. math:: 
    
    \ce{A}^\ast & \ce{<=> B}^\ast \\
    \ce{A}^\ast & \ce{<=> B + C}


There are a number of ways that an activated species :math:`\ce{A}^\ast` can 
be produced:

* **Chemical activation**. :math:`\ce{A}^\ast` is produced as the adduct of
  an association reaction:
  
  .. math:: \ce{B + C <=> A}^\ast

* **Thermal activation**. :math:`\ce{A}^\ast` is produced via transfer of
  energy from an otherwise inert species :math:`\ce{M}` via bimolecular
  collision:
  
  .. math:: \ce{A + M <=> A}^\ast \ce{\mbox{} + M}

* **Photoactivation**. :math:`\ce{A}^\ast` is produced as a result of 
  absorption of a photon:
  
  .. math:: \ce{A} + h \nu \ce{-> A}^\ast

Once an activated molecule has been produced, multiple isomerization and
dissociation reactions may become competitive with one another and with
collisional stabilization (thermal deactivation); these combine to form a
network of unimolecular reactions. The major pathway will depend on the
relative rates of collision and reaction, which in turn is a function of
both temperature and pressure. At high pressure the collision rate will be
fast, and activated molecules will tend to be collisionally stabilized before
reactive events can occur; this is called the *high-pressure limit*. At low
pressures the collision rate will be slow, and activated molecules will
tend to isomerize and dissociate, often traversing multiple reactive events
before collisional stabilization can occur.

The onset of the pressure-dependent regime varies with both temperature and
molecular size. The figure below shows the approximate pressure at which 
pressure-dependence becomes important as a function of temperature and
molecular size. The parameter :math:`m \equiv N_\mathrm{vib} + \frac{1}{2} N_\mathrm{rot}`
represents a count of the internal degrees of freedom (vibrations and hindered
rotors, respectively). The ranges of the x-axis and y-axis suggest that
pressure dependence is in fact important over a wide regime of conditions of
practical interest, particularly in high-temperature processes such as
pyrolysis and combustion [Wong2003]_.

.. figure:: images/switchover_pressure.*
    :width: 67%
    
    Plot of the switchover pressure -- indicating the onset of pressure
    dependence -- as a function of temperature and molecular size. The
    value :math:`m \equiv N_\mathrm{vib} + \frac{1}{2} N_\mathrm{rot}`
    represents a count of the internal degrees of freedom. Over a wide
    variety of conditions of practical interest, even very large
    molecules exhibit significant pressure dependence. Figure adapted from
    Wong, Matheu and Green (2003).

.. [Wong2003] B. M. Wong, D. M. Matheu, and W. H. Green. *J. Phys. Chem. A* 
   **107**, p. 6206-6211 (2003).
   `doi:10.1021/jp034165g <http://dx.doi.org/10.1021/jp034165g>`_


Historical Context
==================

The importance of bimolecular collisions in unimolecular reactions was first
proposed by Lindemann in 1922 [Lindemann1922]_. It was soon recognized by
Hinshelwood and others that a rigorous treatment of these processes required
consideration of molecular energy levels [Hinshelwood1926]_. The RRKM
expression for the microcanonical rate coefficient $k(E)$ was derived in the
early 1950s [Rice1927]_ [Kassel1928]_ [Marcus1951]_. In the late 1950s master
equation models of chemical systems began appearing [Siegert1949]_
[Bartholomay1958]_ [Montroll1958]_ [Krieger1960]_ [Gans1960]_, including an
early linear integral-differential equation formulation by Widom [Widom1959]_.
Analytical solutions for a variety of simple models soon followed [Keck1965]_
[Troe1967]_ [Troe1973]_, as did the first numerical approaches [Tardy1966]_.
Numerical methods -- which are required for complex unimolecular reaction
networks -- became much more attractive in the 1970s with the appearance of
new algorithms, including Gear's method for solving stiff systems of ordinary
differential equations [Gear1971]_ and efficient algorithms for calculating
the density of states [Beyer1973]_ [Stein1973]_ [Astholz1979]_. In the 1990s
computing power had increased to the point where it was practical to solve
them numerically by discretizing the integrals over energy.

.. [Lindemann1922] F. A. Lindemann. *Trans. Faraday Soc.* **17**, 
   p. 598-606 (1922).

.. [Hinshelwood1926] C. N. Hinshelwood. *Proc. Royal Soc. A* **17**,
   p. 230-233 (1926).
   `JSTOR:94593 <http://www.jstor.org/stable/94593>`_

.. [Rice1927] O. K. Rice and H. C. Ramsperger. *J. Am. Chem. Soc.* **49**,
   p. 1617-1629 (1927).
   `doi:10.1021/ja01406a001 <http://dx.doi.org/10.1021/ja01406a001>`_

.. [Kassel1928] L. S. Kassel. *J. Phys. Chem.* **32**, 
   p. 1065-1079 (1928).
   `doi:10.1021/j150289a011 <http://dx.doi.org/10.1021/j150289a011>`_

.. [Marcus1951] R. A. Marcus and O. K. Rice. *J. Phys. Coll. Chem.* **55**,
   p. 894-908 (1951).
   `doi:10.1021/j150489a013 <http://dx.doi.org/10.1021/j150489a013>`_

.. [Siegert1949] A. J. F. Siegert. *Phys. Rev.* **76**,
   p. 1708-1714 (1949).
   `doi:10.1103/PhysRev.76.1708 <http://dx.doi.org/10.1103/PhysRev.76.1708>`_

.. [Bartholomay1958] A. F. Bartholomay. *Bull. Math. Biophys.* **20**,
   p. 175-190 (1958).
   `doi:10.1007/BF02478297 <http://dx.doi.org/10.1007/BF02478297>`_
   
.. [Montroll1958] E. W. Montroll and K. E. Shuler. *Adv. Chem. Phys.* **1**,
   p. 361-399 (1958).

.. [Krieger1960] I. M. Krieger and P. J. Gans. *J. Chem. Phys.* **32**,
   p. 247-250 (1960).
   `doi:10.1063/1.1700909 <http://dx.doi.org/10.1063/1.1700909>`_

.. [Gans1960] P. J. Gans. *J. Chem. Phys.* **33**,
   p. 691-694 (1960).
   `doi:10.1063/1.1731239 <http://dx.doi.org/10.1063/1.1731239>`_

.. [Widom1959] B. Widom. *J. Chem. Phys.* **31**, 
   p. 1387-1394 (1959).
   `doi:10.1063/1.1730604 <http://dx.doi.org/10.1063/1.1730604>`_

.. [Keck1965] J. Keck and G. Carrier. *J. Chem. Phys.* **43**, 
   p. 2284-2298 (1965).
   `doi:10.1063/1.1697125 <http://dx.doi.org/10.1063/1.1697125>`_

.. [Troe1967] J. Troe and H. Gg. Wagner. *Ber. Bunsenges. Phys. Chem.* **71**,
   p. 937 (1967).
   `doi:10.1002/bbpc.19670710904 <http://dx.doi.org/10.1002/bbpc.19670710904>`_

.. [Troe1973] J. Troe. *Ber. Bunsenges. Phys. Chem.* **77**,
   p. 665 (1973).
   `doi:10.1002/bbpc.19730770903 <http://dx.doi.org/10.1002/bbpc.19730770903>`_

.. [Tardy1966] D. C. Tardy and B. S. Rabinovitch. *J. Chem. Phys.*
   **45**, p. 3720-3730 (1966).
   `doi:10.1063/1.1727392 <http://dx.doi.org/10.1063/1.1727392>`_

.. [Gear1971] C. W. Gear. *Commun. ACM* **14**,
   p. 176-179 (1971).
   `doi:10.1145/362566.362571 <http://dx.doi.org/10.1145/362566.362571>`_

.. [Beyer1973] T. Beyer and D. F. Swinehart. *Commun. ACM* **16**,
   p. 379 (1973).
   `doi:10.1145/362248.362275 <http://dx.doi.org/10.1145/362248.362275>`_

.. [Stein1973] S. E. Stein and B. S. Rabinovitch. *J. Chem. Phys.* **58**,
   p. 2438-2444 (1973).
   `doi:10.1063/1.1679522 <http://dx.doi.org/10.1063/1.1679522>`_

.. [Astholz1979] D. C. Astholz, J. Troe, and W. Wieters. *J. Chem. Phys.* 
   **70**, p. 5107-5116 (1979).
   `doi:10.1063/1.437352 <http://dx.doi.org/10.1063/1.437352>`_

