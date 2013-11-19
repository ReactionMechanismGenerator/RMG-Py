"""
webel - A Cinfony module that runs entirely on web services

webel can be used from all of CPython, Jython and IronPython. 

Global variables:
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  fps - a list of supported fingerprint types
"""

import re
import os
import urllib2
import StringIO

try:
    import Tkinter as tk
    import Image as PIL
    import ImageTk as piltk
except ImportError:
    tk = None

informats = {"smi":"SMILES", "inchikey":"InChIKey", "inchi":"InChI",
             "name":"Common name"}
"""A dictionary of supported input formats"""
outformats = {"smi":"SMILES", "cdxml":"ChemDraw XML", "inchi":"InChI",
              "sdf":"Symyx SDF", "names":"Common names", "inchikey":"InChIKey",
              "alc":"Alchemy", "cerius":"MSI Cerius II", "charmm":"CHARMM",
              "cif":"Crystallographic Information File",
              "cml":"Chemical Markup Language", "ctx":"Gasteiger Clear Text",
              "gjf":"Gaussian job file", "gromacs":"GROMACS",
              "hyperchem":"HyperChem", "jme":"Java Molecule Editor",
              "maestro":"Schrodinger MacroModel",
              "mol":"Symyx mol", "mol2":"Tripos Sybyl MOL2",
              "mrv":"ChemAxon MRV", "pdb":"Protein Data Bank",
              "sdf3000":"Symyx SDF3000", "sln":"Sybl line notation",
              "xyz":"XYZ", "iupac":"IUPAC name"}
"""A dictionary of supported output formats"""

fps = ["std", "maccs", "estate"]
"""A list of supported fingerprint types"""

# The following function is taken from urllib.py in the IronPython dist
def _quo(text, safe="/"):
    always_safe = ('ABCDEFGHIJKLMNOPQRSTUVWXYZ'
               'abcdefghijklmnopqrstuvwxyz'
               '0123456789' '_.-')
    _safemaps = {}
    cachekey = (safe, always_safe)
    try:
        safe_map = _safemaps[cachekey]
    except KeyError:
        safe += always_safe
        safe_map = {}
        for i in range(256):
            c = chr(i)
            safe_map[c] = (c in safe) and c or ('%%%02X' % i)
        _safemaps[cachekey] = safe_map
    res = map(safe_map.__getitem__, text)
    return ''.join(res)

def _makeserver(serverurl):
    """Curry the name of the server"""
    def server(*urlcomponents):       
        url = "%s/" % serverurl + "/".join(urlcomponents)
        resp = urllib2.urlopen(url)
        return resp.read()
    return server

rajweb = _makeserver("http://ws1.bmc.uu.se:8182/cdk")
nci = _makeserver("http://cactus.nci.nih.gov/chemical/structure")

_descs = None # Cache the list of descriptors
def getdescs():
    """Return a list of supported descriptor types"""
    global _descs
    if not _descs:
        response = rajweb("descriptors").rstrip()
        _descs = [x.split(".")[-1] for x in response.split("\n")]
    return _descs

def readstring(format, string):
    """Read in a molecule from a string.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       string

    Note: For InChIKeys a list of molecules is returned.

    Example:
    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi", input)   
    """
    format = format.lower()
    if not format in informats:
        raise ValueError("%s is not a recognised Webel format" % format)
    
    if format != "smi":
        smiles = nci(_quo(string), "smiles").rstrip()
    else:
        smiles = string
    if format == "inchikey":
        return [Molecule(smile) for smile in smiles.split("\n")]
    else:
        mol = Molecule(smiles)
        if format == "name":
            mol.title = string
        return mol

class Outputfile(object):
    """Represent a file to which *output* is to be sent.
   
    Although it's possible to write a single molecule to a file by
    calling the write() method of a molecule, if multiple molecules
    are to be written to the same file you should use the Outputfile
    class.
    
    Required parameters:
       format - see the outformats variable for a list of available
                output formats
       filename

    Optional parameters:
       overwrite -- if the output file already exists, should it
                   be overwritten? (default is False)
                   
    Methods:
       write(molecule)
       close()
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format.lower()
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % self.filename)
        if not format in outformats:
            raise ValueError("%s is not a recognised Webel format" % format)
        self.file = open(filename, "w")
    
    def write(self, molecule):
        """Write a molecule to the output file.
        
        Required parameters:
           molecule
        """
        if self.file.closed:
            raise IOError("Outputfile instance is closed.")
        output = molecule.write(self.format)
        print >> self.file, output

    def close(self):
        """Close the Outputfile to further writing."""
        self.file.close()

class Molecule(object):
    """Represent a Webel Molecule.

    Required parameter:
       smiles -- a SMILES string or any type of cinfony Molecule
 
    Attributes:
       formula, molwt, title
    
    Methods:
       calcfp(), calcdesc(), draw(), write() 
      
    The underlying SMILES string can be accessed using the attribute:
       smiles
    """
    _cinfony = True

    def __init__(self, smiles):
        
        if hasattr(smiles, "_cinfony"):
            a, b = smiles._exchange
            if a == 0:
                smiles = b
            else:
                # Must convert to SMILES
                smiles = smiles.write("smi").split()[0]            

        self.smiles = smiles
        self.title = ""
 
    @property
    def formula(self): return rajweb("mf", _quo(self.smiles))
    @property
    def molwt(self): return float(rajweb("mw", _quo(self.smiles)))
    @property
    def _exchange(self):
        return (0, self.smiles)

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, all available descriptors are
        calculated. See the descs variable for a list of available
        descriptors.
        """
        if not descnames:
            descnames = getdescs()
        else:
            for descname in descnames:
                if descname not in getdescs():
                    raise ValueError("%s is not a recognised Webel descriptor type" % descname)
        ans = {}
        p = re.compile("""Descriptor parent="(\w*)" name="([\w\-\+\d]*)" value="([\d\.]*)""")
        for descname in descnames:
            longname = "org.openscience.cdk.qsar.descriptors.molecular." + descname
            response = rajweb("descriptor", longname, _quo(self.smiles))
            for match in p.findall(response):
                if match[2]:
                    ans["%s_%s" % (match[0], match[1])] = float(match[2])
        return ans
    
    def calcfp(self, fptype="std"):
        """Calculate a molecular fingerprint.
        
        Optional parameters:
           fptype -- the fingerprint type (default is "std"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        fptype = fptype.lower()
        if fptype not in fps:
            raise ValueError("%s is not a recognised Webel Fingerprint type" % fptype)
        fp = rajweb("fingerprint/%s/%s" % (fptype, _quo(self.smiles))).rstrip()
        return Fingerprint(fp)

    def write(self, format="smi", filename=None, overwrite=False):
        """Write the molecule to a file or return a string.
        
        Optional parameters:
           format -- see the informats variable for a list of available
                     output formats (default is "smi")
           filename -- default is None
           overwite -- if the output file already exists, should it
                       be overwritten? (default is False)

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.

        To write multiple molecules to the same file you should use
        the Outputfile class.
        """
        format = format.lower()
        if not format in outformats:
            raise ValueError("%s is not a recognised Webel format" % format)
        if format == "smi":
            output = self.smiles
        elif format == "names":
            try:
                output = nci(_quo(self.smiles), "%s" % format).rstrip().split("\n")
            except urllib2.URLError, e:
                if e.code == 404:
                    output = []
        elif format in ['inchi', 'inchikey']:
            format = "std" + format
            output = nci(_quo(self.smiles), "%s" % format).rstrip()
        elif format == 'iupac':
            format = format + "_name"
            try:
                output = nci(_quo(self.smiles), "%s" % format).rstrip()
            except urllib2.URLError, e:
                if e.code == 404:
                    output = ""
        else:
            output = nci(_quo(self.smiles), "file?format=%s" % format).rstrip()

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % filename)
            outputfile = open(filename, "w")
            print >> outputfile, output
            outputfile.close()
        else:
            return output

    def __str__(self):
        return self.write()

    def draw(self, show=True, filename=None):
        """Create a 2D depiction of the molecule.

        Optional parameters:
          show -- display on screen (default is True)
          filename -- write to file (default is None)

        Tkinter and Python Imaging Library are required for
        image display.
        """
        imagedata = nci(_quo(self.smiles), "image")
        if filename:
            print >> open(filename, "wb"), imagedata
        if show:
            if not tk:
                errormessage = ("Tkinter or Python Imaging "
                                "Library not found, but is required for image "
                                "display. See installation instructions for "
                                "more information.")
                raise ImportError, errormessage                 
            root = tk.Tk()
            root.title(self.smiles)
            frame = tk.Frame(root, colormap="new", visual='truecolor').pack()
            image = PIL.open(StringIO.StringIO(imagedata))
            imagedata = piltk.PhotoImage(image)
            label = tk.Label(frame, image=imagedata).pack()
            quitbutton = tk.Button(root, text="Close", command=root.destroy).pack(fill=tk.X)
            root.mainloop()

class Fingerprint(object):
    """A Molecular Fingerprint.
    
    Required parameters:
       fingerprint -- a string of 0's and 1's representing a binary fingerprint

    Attributes:
       fp -- the underlying fingerprint object
       bits -- a list of bits set in the Fingerprint

    Methods:
       The "|" operator can be used to calculate the Tanimoto coeff. For example,
       given two Fingerprints 'a', and 'b', the Tanimoto coefficient is given by:
          tanimoto = a | b
    """
    def __init__(self, fingerprint):
        self.fp = fingerprint
    def __or__(self, other):
        mybits = set(self.bits)
        otherbits = set(other.bits)
        return len(mybits&otherbits) / float(len(mybits|otherbits))
    @property
    def bits(self):
        return [i for i,x in enumerate(self.fp) if x=="1"]
    def __str__(self):
        return self.fp

class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern
    
    Methods:
       match(molecule)
    
    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> smarts.match(mol) 
    True
    """
    def __init__(self, smartspattern):
        """Initialise with a SMARTS pattern."""
        self.pat = smartspattern
    def match(self, molecule):
        """Does a SMARTS pattern match a particular molecule?
        
        Required parameters:
           molecule
        """
        resp = rajweb("substruct", _quo(molecule.smiles), _quo(self.pat)).rstrip()
        return resp == "true"
 
if __name__=="__main__": #pragma: no cover
    import doctest
    doctest.run_docstring_examples(rajweb, globals())
