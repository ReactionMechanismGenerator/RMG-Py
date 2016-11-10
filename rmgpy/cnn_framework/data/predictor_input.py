predictor_model(embedding_size=512, attribute_vector_size=9, depth=2, scale_output=0.05, padding=True, 
				hidden=0, hidden_activation='tanh',
				output_activation='linear', output_size=1, 
				lr=0.01, optimizer='adam', loss='mse')