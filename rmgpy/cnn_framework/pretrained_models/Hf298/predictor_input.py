# added atom_in_rings feature
# added bond_in_rings feature
# fixed a bug in GraphFP updating attribute matrix
# change depth to 3
# add fixed shuffle_seed for data split
# standard data splitting scheme: 10% testing and 5-folds for 90%
# swith to MoleculeConv, RMG-Py sha=33c3534a90e27, rebase to b9d2dec4a0
# added padding, RMG-Py sha=fdfab62242ebef
# in-house train_model allows batch training, sha=b9c862d7e924

predictor_model(prediction_task="Hf298(kcal/mol)",
				embedding_size=512, depth=3, scale_output=0.05, 
				padding=True, padding_final_size=25,
				add_extra_atom_attribute=True, add_extra_bond_attribute=True,
				differentiate_atom_type=False,
				differentiate_bond_type=False,
				mol_conv_inner_activation='tanh',
				mol_conv_outer_activation='softmax',
				hidden=50, hidden_activation='tanh',
				output_activation='linear', output_size=1, 
				lr=0.01, optimizer='adam', loss='mse')

