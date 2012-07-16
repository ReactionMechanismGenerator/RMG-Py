

import base 

class InputWriter(base.InputWriter):

    script_attempts = 5
    max_attempts = 2 * script_attempts

    "Keywords for the multiplicity"
    multiplicityKeywords = {}
    multiplicityKeywords[1] = ''
    multiplicityKeywords[2] = 'uhf doublet'
    multiplicityKeywords[3] = 'uhf triplet'
    multiplicityKeywords[4] = 'uhf quartet'
    multiplicityKeywords[5] = 'uhf quintet'
    multiplicityKeywords[6] = 'uhf sextet'
    multiplicityKeywords[7] = 'uhf septet'
    multiplicityKeywords[8] = 'uhf octet'
    multiplicityKeywords[9] = 'uhf nonet'

    "Keywords that will be added at the top of the qm input file"
    keywordsTop = {}
    keywordsTop[1] = "precise nosym"
    keywordsTop[2] = "precise nosym gnorm=0.0 nonr"
    keywordsTop[3] = "precise nosym gnorm=0.0"
    keywordsTop[4] = "precise nosym gnorm=0.0 bfgs"
    keywordsTop[5] = "precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000"

    "Keywords that will be added at the bottom of the qm input file"
    keywordsBottom = {}
    keywordsBottom[1] = "oldgeo thermo nosym precise "
    keywordsBottom[2] = "oldgeo thermo nosym precise "
    keywordsBottom[3] = "oldgeo thermo nosym precise "
    keywordsBottom[4] = "oldgeo thermo nosym precise "
    keywordsBottom[5] = "oldgeo thermo nosym precise "


    def writeInputFile(self, geometry, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        multiplicity_keyword = this.multiplicityKeywords[geometry.multiplicity]

        top_keywords = "pm3 {0} {1}".format(
                multiplicity_keyword,
                self.keywordsTop[attempt],
                )
        bottom_keywords = "{0} pm3 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keyword,
                )
        polar_keywords = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if geometry.multiplicity == 1 else 'static',
                multiplicity_keyword,
                )

        inputFilePath = os.path.join( where_stuff_lives , geometry.uniqueID + self.inputExtension)

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()

        if attempt <= self.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, geometry.getRefinedMolFilePath() )
        else:
            obConversion.ReadFile(mol, geometry.getCrudeMolFilePath() )

        mol.SetTitle( geometry.uniqueID ) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)

        input_string = obConversion.WriteString(mol)

        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keywords)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keywords)
            if qmtp.QMTP.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keywords)

        return self.molfile.name + '.mop'


class Verifier(base.Verifier):
    def checkJob():
        # Check the job was completed

class Job(base.Job):
    def runJob(self):
        # Run the job. Returns 'True' if the job successfully ran, 'False' if not.
        job = qm.jobs.MOPACJob(self)

        return self.checkJob()

    def parseResult(self):
        # cclib, Quixote?
        qmData = qm.parsers.CCLibParser.parse

        return qmData

    def processResult(self):
        # package thermo, calculator

        return self.thermo

