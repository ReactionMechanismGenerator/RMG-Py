predictor_model(embedding_size=512, depth=2, 
				add_extra_atom_attribute=True, add_extra_bond_attribute=False,
				scale_output=0.05, padding=False, 
				hidden=50, hidden_activation='tanh',
				output_activation='linear', output_size=1, 
				lr=0.01, optimizer='adam', loss='mse')