#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module provides the :class:`PressureDependenceJob` class, which represents
a job for computing the pressure-dependent rate coefficients of a unimolecular
reaction network.
"""

import os.path
import math
import numpy
import logging
import rmgpy.quantity as quantity
from rmgpy.kinetics import Chebyshev, PDepArrhenius
from rmgpy.reaction import Reaction
from rmgpy.kinetics.tunneling import Wigner, Eckart
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.cantherm.output import prettify
from rmgpy.chemkin import writeKineticsEntry
from sensitivity import PDepSensitivity as sa
from rmgpy.exceptions import InvalidMicrocanonicalRateError, ModifiedStrongCollisionError

################################################################################

class PressureDependenceJob(object):
    """
    A representation of a pressure dependence job. The attributes are:
    
    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `Tmin`                  The minimum temperature at which to compute :math:`k(T,P)` values
    `Tmax`                  The maximum temperature at which to compute :math:`k(T,P)` values
    `Tcount`                The number of temperatures at which to compute :math:`k(T,P)` values
    `Pmin`                  The minimum pressure at which to compute :math:`k(T,P)` values
    `Pmax`                  The maximum pressure at which to compute :math:`k(T,P)` values
    `Pcount`                The number of pressures at which to compute :math:`k(T,P)` values
    `Emin`                  The minimum energy to use to compute :math:`k(T,P)` values
    `Emax`                  The maximum energy to use to compute :math:`k(T,P)` values
    `maximumGrainSize`      The maximum energy grain size to use to compute :math:`k(T,P)` values
    `minimumGrainCount`     The minimum number of energy grains to use to compute :math:`k(T,P)` values
    `method`                The method to use to reduce the master equation to :math:`k(T,P)` values
    `interpolationModel`    The interpolation model to fit to the computed :math:`k(T,P)` values
    `maximumAtoms`          The maximum number of atoms to apply pressure dependence to (in RMG jobs)
    `activeKRotor`          A flag indicating whether to treat the K-rotor as active or adiabatic
    `activeJRotor`          A flag indicating whether to treat the J-rotor as active or adiabatic
    `rmgmode`               A flag that toggles "RMG mode", described below
    ----------------------- ----------------------------------------------------
    `network`               The unimolecular reaction network
    `Tlist`                 An array of temperatures at which to compute :math:`k(T,P)` values
    `Plist`                 An array of pressures at which to compute :math:`k(T,P)` values
    `Elist`                 An array of energies to use to compute :math:`k(T,P)` values
    ======================= ====================================================
    
    In RMG mode, several alterations to the k(T,P) algorithm are made both for
    speed and due to the nature of the approximations used:
    
    * Densities of states are not computed for product channels
    
    * Arbitrary rigid rotor moments of inertia are included in the active modes;
      these cancel in the ILT and equilibrium expressions
    
    * k(E) for each path reaction is computed in the direction A -> products,
      where A is always an explored isomer; the high-P kinetics are reversed
      if necessary for this purpose
    
    * Thermodynamic parameters are always used to compute the reverse k(E)
      from the forward k(E) for each path reaction
    
    RMG mode should be turned off by default except in RMG jobs.    
    """
    
    def __init__(self, network, 
        Tmin=None, Tmax=None, Tcount=0, Tlist=None,
        Pmin=None, Pmax=None, Pcount=0, Plist=None,
        maximumGrainSize=None, minimumGrainCount=0,
        method=None, interpolationModel=None, maximumAtoms=None,
        activeKRotor=True, activeJRotor=True, rmgmode=False, sensitivity_conditions=None):
        self.network = network
        
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Tcount = Tcount
        if Tlist is not None:
            self.Tlist = Tlist
            self.Tmin = (numpy.min(self.Tlist.value_si),"K")
            self.Tmax = (numpy.max(self.Tlist.value_si),"K")
            self.Tcount = len(self.Tlist.value_si)
        else:
            self.Tlist = None

        self.Pmin = Pmin
        self.Pmax = Pmax
        self.Pcount = Pcount
        if Plist is not None:
            self.Plist = Plist
            self.Pmin = (numpy.min(self.Plist.value_si)*1e-5,"bar")
            self.Pmax = (numpy.max(self.Plist.value_si)*1e-5,"bar")
            self.Pcount = len(self.Plist.value_si)
        else:
            self.Plist = None

        self.maximumGrainSize = maximumGrainSize
        self.minimumGrainCount = minimumGrainCount
        self.Emin = None
        self.Emax = None
        self.Elist = None
        
        self.method = method
        self.interpolationModel = interpolationModel
        self.maximumAtoms = maximumAtoms
        
        self.activeKRotor = activeKRotor
        self.activeJRotor = activeJRotor
        self.rmgmode = rmgmode

        if sensitivity_conditions is not None:
            if not isinstance(sensitivity_conditions[0], list):
                sensitivity_conditions = [sensitivity_conditions]  # allow `[T, P]` as conditions input
            self.sensitivity_conditions = [[quantity.Quantity(condition[0]), quantity.Quantity(condition[1])]
                                           for condition in sensitivity_conditions]
        else:
            self.sensitivity_conditions = None
        
        if self.Tlist is None and self.Tmin is not None and self.Tmax is not None and self.Tcount is not None:
            self.generateTemperatureList()
        if self.Plist is None and self.Pmin is not None and self.Pmax is not None and self.Pcount is not None:
            self.generatePressureList()
        
    @property
    def Tmin(self):
        """The minimum temperature at which the computed k(T,P) values are valid, or ``None`` if not defined."""
        return self._Tmin
    @Tmin.setter
    def Tmin(self, value):
        self._Tmin = quantity.Temperature(value)
    
    @property
    def Tmax(self):
        """The maximum temperature at which the computed k(T,P) values are valid, or ``None`` if not defined."""
        return self._Tmax
    @Tmax.setter
    def Tmax(self, value):
        self._Tmax = quantity.Temperature(value)
    
    @property
    def Tlist(self):
        """The temperatures at which the k(T,P) values are computed."""
        return self._Tlist
    @Tlist.setter
    def Tlist(self, value):
        self._Tlist = quantity.Temperature(value)
    
    @property
    def Pmin(self):
        """The minimum pressure at which the computed k(T,P) values are valid, or ``None`` if not defined."""
        return self._Pmin
    @Pmin.setter
    def Pmin(self, value):
        self._Pmin = quantity.Pressure(value)
    
    @property
    def Pmax(self):
        """The maximum pressure at which the computed k(T,P) values are valid, or ``None`` if not defined."""
        return self._Pmax
    @Pmax.setter
    def Pmax(self, value):
        self._Pmax = quantity.Pressure(value)
    
    @property
    def Plist(self):
        """The pressures at which the k(T,P) values are computed."""
        return self._Plist
    @Plist.setter
    def Plist(self, value):
        self._Plist = quantity.Pressure(value)

    @property
    def maximumGrainSize(self):
        """The maximum allowed energy grain size, or ``None`` if not defined."""
        return self._maximumGrainSize
    @maximumGrainSize.setter
    def maximumGrainSize(self, value):
        self._maximumGrainSize = quantity.Energy(value)

    def copy(self):
        """
        Return a copy of the pressure dependence job.
        """
        return PressureDependenceJob(
            network = self.network, 
            Tmin = self.Tmax, 
            Tmax = self.Tmax, 
            Tcount = self.Tcount, 
            Tlist = self.Tlist,
            Pmin = self.Pmin, 
            Pmax = self.Pmax, 
            Pcount = self.Pcount, 
            Plist = self.Plist,
            maximumGrainSize = self.maximumGrainSize, 
            minimumGrainCount = self.minimumGrainCount,
            method = self.method, 
            interpolationModel = self.interpolationModel,
            activeKRotor = self.activeKRotor, 
            activeJRotor = self.activeJRotor,
            rmgmode = self.rmgmode,
        )

    def execute(self, outputFile, plot, format='pdf', print_summary=True):
        
        for config in self.network.isomers + self.network.reactants + self.network.products:
            for spec in config.species:
                if spec.conformer.E0 is None:
                    raise AttributeError('species {0} is missing energy for its conformer'.format(spec.label))

        # set transition state Energy if not set previously using same method as RMG pdep
        for reaction in self.network.pathReactions:
            transitionState = reaction.transitionState
            if transitionState.conformer and transitionState.conformer.E0 is None:
                transitionState.conformer.E0 = (sum([spec.conformer.E0.value_si for spec in reaction.reactants]) + reaction.kinetics.Ea.value_si,"J/mol")
                logging.info('Approximated transitions state E0 for reaction {3} from kinetics '
                             'A={0}, n={1}, Ea={2} J/mol'.format(reaction.kinetics.A.value_si,reaction.kinetics.n.value_si,reaction.kinetics.Ea.value_si,reaction.label))
        if print_summary:
            self.network.printSummary()
        
        if outputFile is not None:
            self.draw(os.path.dirname(outputFile), format)
        
        self.initialize()
        
        self.K = self.network.calculateRateCoefficients(self.Tlist.value_si, self.Plist.value_si, self.method)

        self.fitInterpolationModels()

        if outputFile is not None:
            self.save(outputFile)
            if plot:
                self.plot(os.path.dirname(outputFile))
            if self.sensitivity_conditions is not None:
                perturbation = 0.1  # kcal/mol
                logging.info('\n\nRunning sensitivity analysis...')
                for i in xrange(3):
                    try:
                        sa(self, os.path.dirname(outputFile), perturbation=perturbation)
                    except (InvalidMicrocanonicalRateError, ModifiedStrongCollisionError) as exept:
                        logging.warn("Could not complete the sensitivity analysis with a perturbation of {0}"
                                     " kcal/mol, trying {1} kcal/mol instead.".format(
                                        perturbation, perturbation / 2.0))
                        perturbation /= 2.0
                    else:
                        break
                else:
                    logging.error("Could not complete the sensitivity analysis even with a perturbation of {0}"
                                  " kcal/mol".format(perturbation))
                    raise exept
                logging.info("Completed the sensitivity analysis using a perturbation of {0} kcal/mol".format(
                    perturbation))
        logging.debug('Finished pdep job for reaction {0}.'.format(self.network.label))
        logging.debug(repr(self.network))

    def generateTemperatureList(self):
        """
        Returns an array of temperatures based on the interpolation `model`,
        minimum and maximum temperatures `Tmin` and `Tmax` in K, and the number of
        temperatures `Tcount`. For Chebyshev polynomials a Gauss-Chebyshev
        distribution is used; for all others a linear distribution on an inverse
        temperature domain is used. Note that the Gauss-Chebyshev grid does *not*
        place `Tmin` and `Tmax` at the endpoints, yet the interpolation is still
        valid up to these values.
        """
        Tmin = self.Tmin.value_si
        Tmax = self.Tmax.value_si
        Tcount = self.Tcount
        if self.Tlist is not None:
            pass
        elif self.interpolationModel[0].lower() == 'chebyshev':
            # Distribute temperatures on a Gauss-Chebyshev grid
            Tlist = numpy.zeros(Tcount, numpy.float64)
            for i in range(Tcount):
                T = -math.cos((2*i+1) * math.pi / (2*self.Tcount))
                T = 2.0 / ((1.0/Tmax - 1.0/Tmin) * T + 1.0/Tmax + 1.0/Tmin)
                Tlist[i] = T
            self.Tlist = (Tlist,"K")
        else:
            # Distribute temperatures evenly on a T^-1 domain
            Tlist = 1.0/numpy.linspace(1.0/Tmax, 1.0/Tmin, Tcount)
            self.Tlist = (Tlist,"K")
        return self.Tlist.value_si
    
    def initialize(self):
        for reaction in self.network.pathReactions:
            tunneling = reaction.transitionState.tunneling
            # throw descriptive error if tunneling not allowed
            if tunneling and reaction.transitionState.frequency is None and reaction.kinetics is not None:
                raise ValueError("""Cannot apply tunneling for reaction {0} when inverse laplace is used.
                                 Either remove tunnelling parameter or input transitionState
                                 frequencies/quantum file""".format(reaction.label))
            # add tunneling parameters
            if isinstance(tunneling, Wigner) and tunneling.frequency is None:
                tunneling.frequency = (reaction.transitionState.frequency.value_si,"cm^-1")
            elif isinstance(tunneling, Eckart) and tunneling.frequency is None:
                tunneling.frequency = (reaction.transitionState.frequency.value_si,"cm^-1")
                tunneling.E0_reac = (sum([reactant.conformer.E0.value_si for reactant in reaction.reactants])*0.001,"kJ/mol")
                tunneling.E0_TS = (reaction.transitionState.conformer.E0.value_si*0.001,"kJ/mol")
                tunneling.E0_prod = (sum([product.conformer.E0.value_si for product in reaction.products])*0.001,"kJ/mol")
            elif tunneling is not None:
                if tunneling.frequency is not None:
                    # Frequency was given by the user
                    pass
                else:
                    raise ValueError('Unknown tunneling model {0!r} for path reaction {1}.'.format(tunneling, reaction))

        maximumGrainSize = self.maximumGrainSize.value_si if self.maximumGrainSize is not None else 0.0
        
        self.network.initialize(
            Tmin = self.Tmin.value_si, 
            Tmax = self.Tmax.value_si, 
            Pmin = self.Pmin.value_si, 
            Pmax = self.Pmax.value_si, 
            maximumGrainSize = maximumGrainSize, 
            minimumGrainCount = self.minimumGrainCount, 
            activeJRotor = self.activeJRotor, 
            activeKRotor = self.activeKRotor, 
            rmgmode = self.rmgmode,
        )

        self.generateTemperatureList()
        self.generatePressureList()
    
    def generatePressureList(self):
        """
        Returns an array of pressures based on the interpolation `model`,
        minimum and maximum pressures `Pmin` and `Pmax` in Pa, and the number of
        pressures `Pcount`. For Chebyshev polynomials a Gauss-Chebyshev
        distribution is used; for all others a linear distribution on an logarithmic
        pressure domain is used. Note that the Gauss-Chebyshev grid does *not*
        place `Pmin` and `Pmax` at the endpoints, yet the interpolation is still
        valid up to these values.
        """
        Pmin = self.Pmin.value_si
        Pmax = self.Pmax.value_si
        Pcount = self.Pcount
        if self.Plist is not None:
            pass
        elif self.interpolationModel[0].lower() == 'chebyshev':
            # Distribute pressures on a Gauss-Chebyshev grid
            Plist = numpy.zeros(Pcount, numpy.float64)
            for i in range(Pcount):
                P = -math.cos((2*i+1) * math.pi / (2*self.Pcount))
                P = 10**(0.5 * ((math.log10(Pmax) - math.log10(Pmin)) * P + math.log10(Pmax) + math.log10(Pmin)))
                Plist[i] = P
            self.Plist = (Plist*1e-5,"bar")
        else:
            # Distribute pressures evenly on a log domain
            Plist = 10.0 ** numpy.linspace(math.log10(Pmin), math.log10(Pmax), Pcount)
            self.Plist = (Plist*1e-5,"bar")
        return self.Plist.value_si

    def fitInterpolationModels(self):
            
        configurations = []
        configurations.extend(self.network.isomers)
        configurations.extend(self.network.reactants)
        configurations.extend(self.network.products)
        
        self.network.netReactions = []
        
        Nreac = self.network.Nisom + self.network.Nreac
        Nprod = Nreac + self.network.Nprod
        
        Tmin = self.Tmin.value_si
        Tmax = self.Tmax.value_si
        Tdata = self.Tlist.value_si
        Pmin = self.Pmin.value_si
        Pmax = self.Pmax.value_si
        Pdata = self.Plist.value_si
        
        for prod in range(Nprod):
            for reac in range(Nreac):
                if reac == prod: continue
                reaction = Reaction(
                    reactants = configurations[reac].species,
                    products = configurations[prod].species,
                )
                
                kdata = self.K[:,:,prod,reac].copy()
                order = len(reaction.reactants)
                kdata *= 1e6 ** (order-1)
                kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
                logging.debug('Fitting master eqn data to kinetics for reaction {}.'.format(reaction))
                reaction.kinetics = self.fitInterpolationModel(Tdata, Pdata, kdata, kunits)
                
                self.network.netReactions.append(reaction)
                
    def fitInterpolationModel(self, Tdata, Pdata, kdata, kunits):
        
        Tmin = self.Tmin.value_si
        Tmax = self.Tmax.value_si
        Pmin = self.Pmin.value_si
        Pmax = self.Pmax.value_si
        
        model = self.interpolationModel[0].lower()
        
        if model == 'chebyshev':
            kinetics = Chebyshev().fitToData(Tdata, Pdata, kdata, kunits,
                self.interpolationModel[1], self.interpolationModel[2],
                Tmin, Tmax, Pmin, Pmax,
            )             
        elif model == 'pdeparrhenius':
            kinetics = PDepArrhenius().fitToData(Tdata, Pdata, kdata, kunits)
        else:
            raise Exception('Invalid interpolation model {0!r}.'.format(self.interpolationModel[0]))
        return kinetics
    
    def save(self, outputFile):
        
        logging.info('Saving pressure dependence results for network {0}...'.format(self.network.label))
        f = open(outputFile, 'a')
        f_chemkin = open(os.path.join(os.path.dirname(outputFile), 'chem.inp'), 'a')
    
        Nreac = self.network.Nisom + self.network.Nreac
        Nprod = Nreac + self.network.Nprod
        Tlist = self.Tlist.value_si
        Plist = self.Plist.value_si
        Tcount = Tlist.shape[0]
        Pcount = Plist.shape[0]
        
        count = 0
        printed_reactions = [] # list of rxns already printed
        for prod in range(Nprod):
            for reac in range(Nreac):
                if reac == prod: continue
                reaction = self.network.netReactions[count]
                count += 1
                # make sure we aren't double counting any reactions
                if not any([reaction.isIsomorphic(other_rxn,checkOnlyLabel=True) \
                            for other_rxn in printed_reactions]):
                    duplicate = False
                    # add reaction to printed reaction
                    printed_reactions.append(reaction)
                else:
                    # comment out the extra reverse reaction
                    duplicate = True

                # write chemkin output.
                string = writeKineticsEntry(reaction, speciesList=None, verbose=False, commented = duplicate)
                f_chemkin.write('{0}\n'.format(string))

                # write to 'output.py'
                kdata = self.K[:,:,prod,reac].copy()
                order = len(reaction.reactants)
                kdata *= 1e6 ** (order-1)
                kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
                
                f.write('#   =========== ')
                f.write('=========== ' * Pcount)
                f.write('\n')
                f.write('#         T \ P ')
                f.write(' '.join(['{0:11.3e}'.format(P*1e-5) for P in Plist]))
                f.write('\n')
                f.write('#   =========== ')
                f.write('=========== ' * Pcount)
                f.write('\n')
                
                for t in range(Tcount):
                    f.write('#   {0:11g}'.format(Tlist[t]))
                    for p in range(Pcount):
                        f.write(' {0:11.3e}'.format(kdata[t,p]))
                    f.write('\n')
                
                f.write('#   =========== ')
                f.write('=========== ' * Pcount)
                f.write('\n')
                
                string = 'pdepreaction(reactants={0!r}, products={1!r}, kinetics={2!r})'.format(
                    [reactant.label for reactant in reaction.reactants],
                    [product.label for product in reaction.products],
                    reaction.kinetics,
                )
                pdep_function = '{0}\n\n'.format(prettify(string))
                if duplicate:
                    # add comments to the start of the string
                    pdep_function = '#   ' + pdep_function.replace('\n','\n#   ')
                f.write(pdep_function)

        f.close()
        f_chemkin.close()

    def plot(self, outputDirectory):

        # Skip this step if matplotlib is not installed
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        import matplotlib.cm
        cm = matplotlib.cm.jet

        Nreac = self.network.Nisom + self.network.Nreac
        Nprod = Nreac + self.network.Nprod
        Tlist = self.Tlist.value_si
        Plist = self.Plist.value_si
        Tcount = Tlist.shape[0]
        Pcount = Plist.shape[0]
        
        K = self.K
        
        count = 0
        for prod in range(Nprod):
            for reac in range(Nreac):
                if reac == prod: continue
                reaction = self.network.netReactions[count]
                count += 1
                
                reaction_str = '{0} {1} {2}'.format(
                    ' + '.join([reactant.label for reactant in reaction.reactants]),
                    '<=>' if prod < Nreac else '-->',
                    ' + '.join([product.label for product in reaction.products]),
                )
                
                fig = plt.figure(figsize=(10,6))
                
                K2 = numpy.zeros((Tcount, Pcount))
                if reaction.kinetics is not None:
                    for t in range(Tcount):
                        for p in range(Pcount):
                            K2[t,p] = reaction.kinetics.getRateCoefficient(Tlist[t], Plist[p])
                
                K = self.K[:,:,prod,reac].copy()
                order = len(reaction.reactants)
                K *= 1e6 ** (order-1)
                K2 *= 1e6 ** (order-1)
                kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]

                plt.subplot(1,2,1)
                for p in xrange(Pcount):
                    plt.semilogy(1000.0 / Tlist, K[:,p], color=cm(1.*p/(Pcount-1)), marker='o', linestyle='',
                                   label=str('%.2e' % (Plist[p]/1e+5)) + ' bar')
                    if reaction.kinetics is not None:
                        plt.semilogy(1000.0 / Tlist, K2[:,p], color=cm(1.*p/(Pcount-1)), marker='', linestyle='-')
                plt.xlabel('1000 / Temperature (1000/K)')
                plt.ylabel('Rate coefficient ({0})'.format(kunits))
                plt.title(reaction_str)
                plt.legend()
                
                plt.subplot(1,2,2)
                for t in xrange(Tcount):
                    plt.loglog(Plist*1e-5, K[t,:], color=cm(1.*t/(Tcount-1)), marker='o', linestyle='',
                                   label=str('%.0d' % Tlist[t]) + ' K')
                    plt.loglog(Plist*1e-5, K2[t,:], color=cm(1.*t/(Tcount-1)), marker='', linestyle='-')
                plt.xlabel('Pressure (bar)')
                plt.ylabel('Rate coefficient ({0})'.format(kunits))
                plt.title(reaction_str)
                plt.legend()
                
                fig.subplots_adjust(left=0.10, bottom=0.13, right=0.95, top=0.92, wspace=0.3, hspace=0.3)
                if not os.path.exists(os.path.join(outputDirectory, 'plots', '')):
                    os.mkdir(os.path.join(outputDirectory, 'plots', ''))
                plt.savefig(os.path.join(outputDirectory, 'plots', 'kinetics_{0:d}.pdf'.format(count)))
                plt.close()

    def draw(self, outputDirectory, format='pdf'):
        """
        Generate a PDF drawing of the pressure-dependent reaction network.
        This requires that Cairo and its Python wrapper be available; if not,
        the drawing is not generated.
        
        You may also generate different formats of drawings, by changing format to 
        one of the following: `pdf`, `svg`, `png`.
        """
        
        # Skip this step if cairo is not installed
        try:
            import cairocffi as cairo
        except ImportError:
            try:
                import cairo
            except ImportError:
                return
        
        from rmgpy.pdep.draw import NetworkDrawer
        
        path = os.path.join(outputDirectory, 'network.' + format)
        
        NetworkDrawer().draw(self.network, format=format, path=path)

    def saveInputFile(self, path):
        """
        Save a CanTherm input file for the pressure dependence job to `path`
        on disk.
        """
        speciesList = self.network.getAllSpecies()
        
        # Add labels for species, reactions, transition states that don't have them
        for i, spec in enumerate(speciesList):
            if not spec.label:
                spec.label = 'species{0:d}'.format(i+1)
        for i, rxn in enumerate(self.network.pathReactions):
            if not rxn.label:
                rxn.label = 'reaction{0:d}'.format(i+1)
            if not rxn.transitionState.label:
                rxn.transitionState.label = 'TS{0:d}'.format(i+1)
        if not self.network.label:
            self.network.label = 'network'
        
        with open(path, 'w') as f:
            # Write species
            for spec in speciesList:
                f.write('species(\n')
                f.write('    label = {0!r},\n'.format(str(spec)))
                if len(spec.molecule) > 0:
                    f.write('    structure = SMILES({0!r}),\n'.format(spec.molecule[0].toSMILES()))
                if spec.conformer is not None:
                    if spec.conformer.E0 is not None:
                        f.write('    E0 = {0!r},\n'.format(spec.conformer.E0))
                    if len(spec.conformer.modes) > 0:
                        f.write('    modes = [\n')
                        for mode in spec.conformer.modes:
                            f.write('        {0!r},\n'.format(mode))
                        f.write('    ],\n')
                    f.write('    spinMultiplicity = {0:d},\n'.format(spec.conformer.spinMultiplicity))
                    f.write('    opticalIsomers = {0:d},\n'.format(spec.conformer.opticalIsomers))
                if spec.molecularWeight is not None:
                    f.write('    molecularWeight = {0!r},\n'.format(spec.molecularWeight))
                if spec.transportData is not None:
                    f.write('    collisionModel = {0!r},\n'.format(spec.transportData))
                if spec.energyTransferModel is not None:
                    f.write('    energyTransferModel = {0!r},\n'.format(spec.energyTransferModel))                    
                if spec.thermo is not None:
                    f.write('    thermo = {0!r},\n'.format(spec.thermo))                    
                f.write(')\n\n')
            
            # Write transition states
            for rxn in self.network.pathReactions:
                ts = rxn.transitionState
                f.write('transitionState(\n')
                f.write('    label = {0!r},\n'.format(ts.label))
                if ts.conformer is not None:
                    if ts.conformer.E0 is not None:
                        f.write('    E0 = {0!r},\n'.format(ts.conformer.E0))
                    if len(ts.conformer.modes) > 0:
                        f.write('    modes = [\n')
                        for mode in ts.conformer.modes:
                            f.write('        {0!r},\n'.format(mode))
                        f.write('    ],\n')
                    f.write('    spinMultiplicity = {0:d},\n'.format(ts.conformer.spinMultiplicity))
                    f.write('    opticalIsomers = {0:d},\n'.format(ts.conformer.opticalIsomers))
                if ts.frequency is not None:
                    f.write('    frequency = {0!r},\n'.format(ts.frequency))                    
                f.write(')\n\n')
                
            # Write reactions
            for rxn in self.network.pathReactions:
                ts = rxn.transitionState
                f.write('reaction(\n')
                f.write('    label = {0!r},\n'.format(rxn.label))
                f.write('    reactants = [{0}],\n'.format(', '.join([repr(str(spec)) for spec in rxn.reactants])))
                f.write('    products = [{0}],\n'.format(', '.join([repr(str(spec)) for spec in rxn.products])))
                f.write('    transitionState = {0!r},\n'.format(rxn.transitionState.label))
                if rxn.kinetics is not None:
                    if isinstance(rxn, LibraryReaction) and 'Reaction library:' not in rxn.kinetics.comment:
                        rxn.kinetics.comment += 'Reaction library: {0!r}'.format(rxn.library)
                    if rxn.network_kinetics is not None:
                        f.write('    kinetics = {0!r},\n'.format(rxn.network_kinetics))
                    else:
                        f.write('    kinetics = {0!r},\n'.format(rxn.kinetics))
                if ts.tunneling is not None:
                    f.write('    tunneling = {0!r},\n'.format(ts.tunneling.__class__.__name__))
                f.write(')\n\n')
            
            # Write network
            f.write('network(\n')
            f.write('    label = {0!r},\n'.format(self.network.label))
            f.write('    isomers = [\n')
            for isomer in self.network.isomers:
                f.write('        {0!r},\n'.format(str(isomer.species[0])))
            f.write('    ],\n')
            f.write('    reactants = [\n')
            for reactants in self.network.reactants:
                f.write('        ({0}),\n'.format(', '.join([repr(str(spec)) for spec in reactants.species])))
            f.write('    ],\n')
            f.write('    bathGas = {\n')
            for spec, frac in self.network.bathGas.items():
                f.write('        {0!r}: {1:g},\n'.format(str(spec), frac))
            f.write('    },\n')
            f.write(')\n\n')
            
            # Write pressure dependence
            f.write('pressureDependence(\n')
            f.write('    label = {0!r},\n'.format(self.network.label))
            f.write('    Tmin = {0!r},\n'.format(self.Tmin))
            f.write('    Tmax = {0!r},\n'.format(self.Tmax))
            f.write('    Tcount = {0:d},\n'.format(self.Tcount))
            f.write('    Tlist = {0!r},\n'.format(self.Tlist))
            f.write('    Pmin = {0!r},\n'.format(self.Pmin))
            f.write('    Pmax = {0!r},\n'.format(self.Pmax))
            f.write('    Pcount = {0:d},\n'.format(self.Pcount))
            f.write('    Plist = {0!r},\n'.format(self.Plist))
            if self.maximumGrainSize is not None:
                f.write('    maximumGrainSize = {0!r},\n'.format(self.maximumGrainSize))
            if self.minimumGrainCount != 0:
                f.write('    minimumGrainCount = {0:d},\n'.format(self.minimumGrainCount))
            f.write('    method = {0!r},\n'.format(self.method))
            if self.interpolationModel is not None:
                f.write('    interpolationModel = {0!r},\n'.format(self.interpolationModel))
            f.write('    activeKRotor = {0!r},\n'.format(self.activeKRotor))
            f.write('    activeJRotor = {0!r},\n'.format(self.activeJRotor))
            if self.rmgmode:
                f.write('    rmgmode = {0!r},\n'.format(self.rmgmode))
            f.write(')\n\n')
