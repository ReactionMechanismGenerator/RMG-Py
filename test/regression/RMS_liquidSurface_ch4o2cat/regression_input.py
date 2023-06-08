
options(
    title='RMS_liquidSurface_ch4o2cat',
    tolerance=0.1
)


observable(
    label='CH4',
    structure=SMILES("[CH4]"),
)

observable(
   label='O2',
   structure=adjacencyList(
       """
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)


reactorSetups(
    reactorTypes=['IdealGasReactor'],
    terminationTimes=([0.01],'s'),
    initialMoleFractionsList=[{
        "CH4": 0.1,
        "O2": 0.2,
        "vacantX": 1.0,
    }],
    temperatures=([600],'K'),
    pressures=([1.0],'bar'),
)