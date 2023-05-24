
options(
    title='RMS_constantVIdealGasReactor_superminimal',
    tolerance=0.1
)

observable(
    label='H2',
    structure=SMILES("[H][H]"),
)

observable(
   label='O2',
   structure=SMILES("[O][O]"),
)

reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([0.01],'s'),
    initialMoleFractionsList=[{
        'H2':.67,
        'O2':.33,
    }],
    temperatures=([1000],'K'),
    pressures=([1.0],'bar'),
)