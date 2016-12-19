****************
Reaction recipes
****************

.. _reaction-recipe-actions:

A reaction recipe is a procedure for applying a reaction to a set of chemical
species. Each reaction recipe is made up of a set of actions that, when applied 
sequentially, a set of chemical reactants to chemical products via that
reaction's characteristic chemical process. Each action requires a small set of
parameters in order to be fully defined.

We define the following reaction recipe actions:

    ============= ============================= ================================
    Action name   Arguments                     Action
    ============= ============================= ================================
    CHANGE_BOND   `center1`, `order`, `center2` change the bond order of the bond between `center1` and `center2` by `order`; do not break or form bonds
    FORM_BOND     `center1`, `order`, `center2` form a new bond between `center1` and `center2` of type `order`
    BREAK_BOND    `center1`, `order`, `center2` break the bond between `center1` and `center2`, which should be of type `order`
    GAIN_RADICAL  `center`, `radical`           increase the number of free electrons on `center` by `radical`
    LOSE_RADICAL  `center`, `radical`           decrease the number of free electrons on `center` by `radical`
    ============= ============================= ================================
