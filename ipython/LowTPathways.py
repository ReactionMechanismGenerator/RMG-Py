#!/usr/bin/env python
# coding: utf-8

# # Generate Species for Low Temperature Pathways
# 
# The goal of this notebook is to allow you to generate all the species that would occur by systematically applying known low-temperature pathways to a starting fuel molecule(s).
# 
# You might do this to generate species that you then put into an RMG input file, to help RMG find the low temperature pathways. An RMG input file species block is generated at the end, to facilitate this.

# These are the pathways, from https://doi.org/10.1016/j.combustflame.2015.07.005
# ![image2.png](attachment:image2.png)

# In[1]:


import sys, os
sys.path.insert(0,os.path.expandvars("$RMGpy"))


# In[2]:


from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy import settings
from IPython.display import display
from arkane.output import prettify


# Declare database variables here by changing the thermo and reaction libraries. These are not actually used.

# In[3]:


database = """
database(
    thermoLibraries = ['BurkeH2O2','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = [
        'H_Abstraction',
        'R_Recombination',
        'R_Addition_MultipleBond',
        'intra_H_migration',
        'HO2_Elimination_from_PeroxyRadical',
        'intra_OH_migration',
    ],
    kineticsEstimator = 'rate rules',
)

options(
    verboseComments=True,  # Set to True for detailed kinetics comments
)
"""


# Just put O2 in the model for now

# In[4]:


speciesList = """
species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]")
)
"""


# In[5]:


# Write input file to disk
os.makedirs('temp', exist_ok=True)
inputFile = open('temp/input.py','w')
inputFile.write(database)
inputFile.write(speciesList)
inputFile.close()


# In[6]:


# initialize RMG instance
# This does all the hard work of loading the databases etc.
# which can take a while
from rmgpy.tools.generate_reactions import RMG
kwargs = {
            'restart': '',
            'walltime': '00:00:00:00',
            'kineticsdatastore': True
    }
rmg = RMG(input_file='temp/input.py', output_directory='temp')

rmg.initialize(**kwargs)


# In[7]:


from rmgpy.molecule import Molecule


# In[8]:


from collections import defaultdict
molecules = defaultdict(set)

def union(*args):
    out = set()
    for a in args:
        out.update(molecules[a])
    return out


# # Put your fuel molecule(s) here:

# In[9]:


# you may have more than one if you wish, just repeat this line
molecules['fuel'].add(Molecule(smiles='CCCCCCCCCC'))


# In[10]:


for m in molecules['fuel']:
    display(m)


# # ☝️☝️☝️☝️☝️☝️☝️☝️

# In[11]:


molecules['H'].add(Molecule(smiles='[H]'))
molecules


# In[12]:


union('fuel','H')


# In[13]:


# React fuel with H via H_Abstraction to get the radicals R

h = list(molecules['H'])[0]
for s in molecules['fuel']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s, h), only_families='H_Abstraction')
    for r in reactions:
        print(r)
        m = r.products[1].molecule[0]
        display(m)
        molecules['R'].add(m)
molecules


# In[14]:


molecules['O2'].add(Molecule(smiles='[O][O]'))
molecules


# In[15]:


# React R with O2 to get the ROO
o2 = list(molecules['O2'])[0]
for s in molecules['R']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s, o2), only_families='R_Recombination')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        display(m)
        molecules['ROO'].add(m)
molecules


# In[16]:


# Isomerize ROO to get QOOH
for s in molecules['ROO']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s, ), only_families='intra_H_migration')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        display(m)
        molecules['QOOH'].add(m)
molecules


# In[17]:


# React QOOH with O2 to get the O2QOOH
o2 = list(molecules['O2'])[0]
for s in molecules['QOOH']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s, o2), only_families='R_Recombination')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        display(m)
        molecules['O2QOOH'].add(m)
molecules


# In[18]:


# What next? OH + keto-hydroperoxide   or   HO2 + alkenyl hydroperoxide 


# In[19]:


# Next we want O2QOOH on its way towards OH + keto-hydroperoxide
# but in elementary steps it goes via HOOQjOOH made by intra_H_migration
# we will need a template to filter out just these species
# from all the others that can be made by intra_H_migration
from rmgpy.molecule import Group
rjooh = Group().from_adjacency_list("""
1 C u1 {2,S}
2 O u0 {1,S} {3,S}
3 O u0 {2,S} {4,S}
4 H u0 {3,S}
""")
rjooh


# In[20]:


# O2QOOH on its way towards OH + keto-hydroperoxide
print("Applying intra_H_migration but finding just the HOOQjOOH radicals")
for s in molecules['O2QOOH']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s,), only_families='intra_H_migration')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        if m.is_subgraph_isomorphic(rjooh):
            display(m)
            molecules['HOOQjOOH'].add(m)
print('HOOQjOOH')
molecules['HOOQjOOH']


# In[21]:


# HOOQjOOH now beta scissions to OH + ketohydroperoxide
print("Applying R_Addition_MultipleBond to HOOQjOOH to find its Beta Scission products")
for s in molecules['HOOQjOOH']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s,), only_families='R_Addition_MultipleBond')
    for r in reactions:
        print(r)
        r1 = r.reactants[0].molecule[0]
        r2 = r.reactants[1].molecule[0]
        if r2.get_formula() == 'HO':
            display(r1)
            molecules['ketohydroperoxide'].add(r1)

print('ketohydroperoxide')
molecules['ketohydroperoxide']


# In[22]:


# The ketohydroperoxide can now break its O-OH and form an alkoxy radical
for s in molecules['ketohydroperoxide']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s,), only_families='R_Recombination')
    for r in reactions:
        print(r)
        m1 = r.reactants[0].molecule[0]
        m2 = r.reactants[1].molecule[0]
        if m1.get_formula() == 'HO':
            display(m2)
            molecules['alkoxy_radical'].add(m2)
            molecules['OH'].add(m1)
        if m2.get_formula() == 'HO':
            display(m1)
            molecules['alkoxy_radical'].add(m1)
            molecules['OH'].add(m2)
molecules['alkoxy_radical']


# In[23]:


# HO2_Elimination_from_PeroxyRadical applied to QOOH makes an Alkene-OOH
print("Applying HO2_Elimination_from_PeroxyRadical")
for s in molecules['O2QOOH']:
    reactions = rmg.database.kinetics.generate_reactions_from_families(
        (s,), only_families='HO2_Elimination_from_PeroxyRadical')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        display(m)
        molecules['QeneOOH'].add(m)
print('QeneOOH')
molecules['QeneOOH']


# In[24]:


# Another chain-brainching route.
# Add H to RO2 to make ROOH (which in the next cell becomes RO and OH)
# This would happen by abstracting an H from some fuel, but 
# here we can make the same species by just adding an H atom.
h = list(molecules['H'])[0]
for s in molecules['ROO']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s,h), only_families='R_Recombination')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        display(m)
        molecules['ROOH'].add(m)
print('ROOH')
molecules['ROOH']


# In[25]:


# The ROOH then becomes RO and OH which is chain-brainching
for s in molecules['ROOH']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s,), only_families='R_Recombination')
    for r in reactions:
        print(r)
        m1 = r.reactants[0].molecule[0]
        m2 = r.reactants[1].molecule[0]
        
        if m1.get_formula() == 'HO':
            display(m2)
            molecules['RO'].add(m2)
            molecules['OH'].add(m1)
        if m2.get_formula() == 'HO':
            display(m1)
            molecules['RO'].add(m1)
            molecules['OH'].add(m2)

print("\n\nRO")
molecules['RO']


# In[26]:


# intra_OH_migration can happen to QOOH
# Not sure it's helpful for low-T chemistry
# So we won't actually save them, but this is what they'd look like
print("Applying intra_OH_migration to QOOH")
for s in molecules['QOOH']:
    reactions = rmg.database.kinetics.generate_reactions_from_families((s,), only_families='intra_OH_migration')
    for r in reactions:
        print(r)
        m = r.products[0].molecule[0]
        display(m)
        # molecules['weird'].add(m)


# In[27]:


# Here is what we've made
molecules


# In[28]:


_=[ display(m) for s in molecules.values() for m in s]


# In[29]:


molecules.keys()


# In[30]:


import datetime
datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# In[31]:



# Print them out for an RMG input file

print(f"""
######################################
# RMG input file species block
# Generated on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# to help discover the low temperature combustion of:
#   {', '.join([m.to_smiles() for m in molecules['fuel']])}
######################################
""")

special = 'OH H O2'.split()

print('#'*30)
print(f'# Small molecules')
for name in special:
    mols = molecules[name]
    for i, m in enumerate(mols):
        print(f"species(label='{name}', reactive=True, structure=SMILES('{m.to_smiles()}'))")
        
for name,mols in molecules.items():
    
    if name in special:
        continue
    
    print('#'*30)
    print(f'# {name}')
    for i,m in enumerate(mols):
        print(f"species(label='{name}_{i+1}', reactive=True, structure=SMILES('{m.to_smiles()}'))")
    


# In[ ]:




