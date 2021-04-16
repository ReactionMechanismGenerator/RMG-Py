#make file for running formatting 
SRC = .

isort = isort $(SRC)
black = black $(SRC)
autoflake = autoflake -ir --remove-all-unused-imports --ignore-init-module-imports --remove-unused-variables   $(SRC)

.PHONY: format
format:
	$(autoflake)
	$(isort)
	$(black)