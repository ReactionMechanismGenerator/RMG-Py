import os
import numpy as np

from rmgpy.thermo import ThermoData
from rmgpy.cnn_framework.predictor import Predictor

class MCNNEstimator():

    """
    A MCNN based estimator for thermochemistry prediction
    
    The attributes are:

    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `Hf298_estimator`   :class:`Predictor`      Hf298 Estimator 
    `S298_estimator`    :class:`Predictor`      S298 Estimator 
    `Cp_estimator`      :class:`Predictor`      Cp Estimator 
    =================== ======================= ====================================

    """

    def __init__(self, Hf298_path, S298_path, Cp_path):
                 
        self.Hf298_estimator = load_pretrained_estimator(Hf298_path)
        self.S298_estimator = load_pretrained_estimator(S298_path)
        self.Cp_estimator = load_pretrained_estimator(Cp_path)

    def get_thermo_data(self, molecule):
        """
        Return thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using mcnn.
        
        Returns: ThermoData
        """
        Tdata = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]
        Cp = self.Cp_estimator.predict(molecule)
        Cp = [np.float64(Cp_i) for Cp_i in Cp]
        Hf298 = self.Hf298_estimator.predict(molecule)
        S298 = self.S298_estimator.predict(molecule)
        comment = "MCNN Estimation."

        Cp0 = molecule.calculateCp0()
        CpInf = molecule.calculateCpInf()
        thermo = ThermoData( 
                           Tdata = (Tdata,"K"),
                           Cpdata = (Cp,"cal/(mol*K)"),
                           H298 = (Hf298,"kcal/mol"),
                           S298 = (S298,"cal/(mol*K)"),
                           Cp0 = (Cp0,"J/(mol*K)"),
                           CpInf = (CpInf,"J/(mol*K)"),
                           Tmin = (300.0,"K"),
                           Tmax = (2000.0,"K"),
                           comment = comment
                          )
        return thermo

def load_pretrained_estimator(model_path):

    estimator = Predictor()
    predictor_input = os.path.join(model_path, 'predictor_input.py')
    param_path = os.path.join(model_path, 'full_train.h5')
    estimator.load_input(predictor_input)
    estimator.load_parameters(param_path)
    return estimator