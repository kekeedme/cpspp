black:
	black .

lint:
	pylint test.py

blt: black lint

.PHONY: black lint blt