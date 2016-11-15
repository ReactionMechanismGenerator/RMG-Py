predictor_model(embedding_size=512, attribute_vector_size=33, depth=2, scale_output=0.05, padding=False, 
				hidden=50, hidden_activation='tanh',
				output_activation='linear', output_size=1, 
				lr=0.01, optimizer='adam', loss='mse')