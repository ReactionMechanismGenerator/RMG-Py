import rmgpy.quantity as quantity
import logging
from rmgpy.species import Species
from rmgpy.data.solvation import SolventData, SoluteData, SoluteGroups, SolvationDatabase
from rmgpy.reaction import Reaction


class DiffusionLimited():

    def __init__(self):
    # default is false, enabled if there is a solvent
        self.enabled = False

    def enable(self, solventData, solvationDatabase, comment=''):
    # diffusionLimiter is enabled if a solvent has been added to the RMG object.
        logging.info("Enabling diffusion-limited kinetics...")
        diffusionLimiter.enabled = True
        diffusionLimiter.database = solvationDatabase
        diffusionLimiter.solventData = solventData

    def getSolventViscosity(self, T):
        return self.solventData.getSolventViscosity(T)
              
    def getEffectiveRate(self, reaction, T):
        """
        Return the ratio of k_eff to k_intrinsic, which is between 0 and 1.
        
        It is 1.0 if diffusion has no effect.
        
        For 1<=>2 reactions, the reverse rate is limited.
        For 2<=>2 reactions, the faster direction is limited.
        For 2<=>1 or 2<=>3 reactions, the forward rate is limited.
        """
        intrinsicKinetics = reaction.kinetics
        reactants = len(reaction.reactants)
        products = len(reaction.products)
        k_forward = intrinsicKinetics.getRateCoefficient(T,P=100e5)
        Keq = reaction.getEquilibriumConstant(T) # Kc
        k_reverse = k_forward / Keq
        k_eff = k_forward
        
        if reactants == 1:
            if products == 1:
                k_eff = k_forward
            else: # two products; reverse rate is limited
                k_diff = self.getDiffusionLimit(T, reaction, forward=False)
                k_eff_reverse = k_reverse*k_diff/(k_reverse+k_diff)
                k_eff = k_eff_reverse * Keq
        else: # 2 reactants
            if products == 1 or products == 3:
                k_diff = self.getDiffusionLimit(T, reaction, forward=True)
                k_eff = k_forward*k_diff/(k_forward+k_diff)
            else: # 2 products
                if Keq > 1.0: # forward rate is faster and thus limited
                    k_diff = self.getDiffusionLimit(T, reaction, forward=True)
                    k_eff = k_forward*k_diff/(k_forward+k_diff)
                else: # reverse rate is faster and thus limited
                    k_diff = self.getDiffusionLimit(T, reaction, forward=False)
                    k_eff_reverse = k_reverse*k_diff/(k_reverse+k_diff)
                    k_eff = k_eff_reverse * Keq
        return k_eff        
    
    def getDiffusionFactor(self, reaction, T):
        """
        Return the diffusion factor of the specified reaction.
        """
        return self.getEffectiveRate(reaction, T)/reaction.kinetics.getRateCoefficient(T,P=0)

    
    def getDiffusionLimit(self, T, reaction, forward=True):
        """
        Return the diffusive limit on the rate coefficient, k_diff.
        
        This is the upper limit on the rate, in the specified direction.
        (ie. forward direction if forward=True [default] or reverse if forward=False)
        Returns the rate coefficient k_diff in m3/mol/s.
        """
        if forward:
            reacting = reaction.reactants
        else:
            reacting = reaction.products
        assert len(reacting)==2, "Can only calculate diffusion limit in a bimolecular direction"
        N_a = 6.022e23 # Avogadro's Number
        radii = 0.0
        diffusivities = 0.0
        for spec in reacting:
            soluteData = self.database.getSoluteData(spec)
            # calculate radius with the McGowan volume and assuming sphere
            radius = ((75 * soluteData.V / 3.14159 / N_a) ** (1. / 3)) / 100  # m
            diff = soluteData.getStokesDiffusivity(T, self.getSolventViscosity(T))
            radii += radius  # meters
            diffusivities += diff #m^2/s
        
        k_diff = 4 * 3.14159 * radii * diffusivities * N_a  # m3/mol/s
        return k_diff


# module level variable. There should only ever be one. It starts off disabled
diffusionLimiter = DiffusionLimited()
