import logging
import rmgpy.constants as constants


class DiffusionLimited():

    def __init__(self):
    # default is false, enabled if there is a solvent
        self.enabled = False

    def enable(self, solvent, solvationDatabase, comment=''):
    # diffusionLimiter is enabled if a solvent has been added to the RMG object.
        logging.info("Enabling diffusion-limited kinetics...")
        diffusionLimiter.enabled = True
        diffusionLimiter.database = solvationDatabase
        diffusionLimiter.solvent = solvent

    def getSolventViscosity(self, T):
        return self.solvent.solventData.getSolventViscosity(T)
              
    def getEffectiveRate(self, reaction, T):
        """
        Return the ratio of keff to kIntrinsic, which is between 0 and 1.
        
        It is 1.0 if diffusion has no effect.
        
        For 1<=>2 reactions, the reverse rate is limited.
        For 2<=>2 reactions, the faster direction is limited.
        For 2<=>1 or 2<=>3 reactions, the forward rate is limited.

        If the solvent species is involved in the reaction, the direction in which the solvent reacts
        is not diffusion-limited.
        """
        intrinsicKinetics = reaction.kinetics
        reactants = len(reaction.reactants)
        products = len(reaction.products)
        isSolventReactant = self.isSolventInvolved(reaction.reactants)
        isSolventProduct = self.isSolventInvolved(reaction.products)
        kf = intrinsicKinetics.getRateCoefficient(T,P=100e5) # forward rate
        Keq = reaction.getEquilibriumConstant(T) # Kc
        kr = kf / Keq #reverse rate
        
        if reactants == 1 or isSolventReactant: # 1 reactant or the solvent is one of the reactants; forward rate not limited
            if products == 1 or isSolventProduct: # 1 product or the solvent is one of the products; reverse rate not limited
                keff = kf
            else: # two products and the solvent is not the product; reverse rate is limited
                kdiff = self.getDiffusionLimit(T, reaction, forward=False) # diffusion rate
                keffReverse = kr*kdiff/(kr+kdiff)
                keff = keffReverse * Keq
        else: # 2 reactants and the solvent is not the reactant
            if products == 1 or isSolventProduct: # 1 product or the solvent is one of the products; forward rate limited
                kdiff = self.getDiffusionLimit(T, reaction, forward=True)
                keff = kf*kdiff/(kf+kdiff)
            else: # 2 products and the solvent is not the product
                if Keq > 1.0: # forward rate is faster and thus limited
                    kdiff = self.getDiffusionLimit(T, reaction, forward=True)
                    keff = kf*kdiff/(kf+kdiff)
                else: # reverse rate is faster and thus limited
                    kdiff = self.getDiffusionLimit(T, reaction, forward=False)
                    keffReverse = kr*kdiff/(kr+kdiff)
                    keff = keffReverse * Keq
        return keff        
    
    def getDiffusionFactor(self, reaction, T):
        """
        Return the diffusion factor of the specified reaction.
        """
        return self.getEffectiveRate(reaction, T)/reaction.kinetics.getRateCoefficient(T,P=0)

    
    def getDiffusionLimit(self, T, reaction, forward=True):
        """
        Return the diffusive limit on the rate coefficient, kdiff.
        
        This is the upper limit on the rate, in the specified direction.
        (ie. forward direction if forward=True [default] or reverse if forward=False)
        Returns the rate coefficient kdiff in m3/mol/s.
        """
        if forward:
            reacting = reaction.reactants
        else:
            reacting = reaction.products
        assert len(reacting)==2, "Can only calculate diffusion limit in a bimolecular direction"
        radii = 0.0
        diffusivities = 0.0
        for spec in reacting:
            soluteData = self.database.getSoluteData(spec)
            # calculate radius with the McGowan volume and assuming sphere
            radius = ((75 * soluteData.V / constants.pi / constants.Na) ** (1. / 3)) / 100  # m
            diff = soluteData.getStokesDiffusivity(T, self.getSolventViscosity(T))
            radii += radius  # meters
            diffusivities += diff #m^2/s
        
        kdiff = 4 * constants.pi * radii * diffusivities * constants.Na  # m3/mol/s
        return kdiff

    def isSolventInvolved(self, spcList):
        """
        Given the list of reactants or products, it checks whether
        the solvent species is included in the list or not.
        It returns "True" if the solvent is in the list and returns "False" if not.
        """
        solventIndex = self.solvent.solventSpecies.index
        return any([(spc.index is solventIndex) for spc in spcList])

# module level variable. There should only ever be one. It starts off disabled
diffusionLimiter = DiffusionLimited()
