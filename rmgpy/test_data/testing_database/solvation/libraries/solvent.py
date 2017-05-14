#!/usr/bin/env python
# encoding: utf-8

name = "Solvent Descriptors"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = 1,
    label = "water",
    molecule = "O",
    solvent = SolventData(
        s_g = 2.743,
        b_g = 4.814,
        e_g = 0.822,
        l_g = -0.213,
        a_g = 3.904,
        c_g = -1.271,
        s_h = 2.836,
        b_h = -41.816,
        e_h = 9.91,
        l_h = -6.354,
        a_h = -32.01,
        c_h = -13.31,
        A = -52.843,
        B = 3703.6,
        C = 5.866,
        D = -5.88e-29,
        E = 10,
        alpha = 0.353,
        beta = 0.38,
        eps = 80.4,
    ),
    shortDesc = u""" """,
    longDesc =
u"""

""",
)

entry(
    index = 2,
    label = "1-octanol",
    molecule = "CCCCCCCCO",
    solvent = SolventData(
        s_g = 0.56,
        b_g = 0.702,
        e_g = -0.203,
        l_g = 0.939,
        a_g = 3.56,
        c_g = -0.12,
        s_h = 5.89,
        b_h = -8.99,
        e_h = -1.04,
        l_h = -9.18,
        a_h = -53.99,
        c_h = -6.49,
        A = -0.022128,
        B = 3018.4,
        C = -2.8054,
        D = 1.3141e-05,
        E = 2,
        alpha = 0.328,
        beta = 0.45,
        eps = 10.3,
    ),
    shortDesc = u""" """,
    longDesc =
u"""
alpha = 0.328, #primary alcohols
beta = 0.45, #primary alcohols,
""",
)

entry(
    index = 3,
    label = "benzene",
    molecule = "C1=CC=CC=C1",
    solvent = SolventData(
        s_g = 1.053,
        b_g = 0.169,
        e_g = -0.313,
        l_g = 1.02,
        a_g = 0.457,
        c_g = 0.107,
        s_h = -12.599,
        b_h = -4.023,
        e_h = 4.446,
        l_h = -8.488,
        a_h = -9.775,
        c_h = -4.637,
        A = 7.5117,
        B = 294.68,
        C = -2.794,
        D = 0,
        E = 0,
        alpha = 0,
        beta = 0.14,
        eps = 2.3,
    ),
    shortDesc = u""" """,
    longDesc =
u"""

""",
)

entry(
    index = 4,
    label = "cyclohexane",
    molecule = "C1CCCCC1",
    solvent = SolventData(
        s_g = 0,
        b_g = 0,
        e_g = -0.11,
        l_g = 1.013,
        a_g = 0,
        c_g = 0.163,
        s_h = 0.0,
        b_h = 0.0,
        e_h = 3.375,
        l_h = -9.078,
        a_h = 0.0,
        c_h = -6.507,
        A = -33.763,
        B = 2497.2,
        C = 3.2236,
        D = 0,
        E = 0,
        alpha = 0,
        beta = 0,
        eps = 2.0,
    ),
    shortDesc = u""" """,
    longDesc =
u"""

""",
)
