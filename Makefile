# 2 liy ki swiv yo vle di ke, le mwen tape `make black` li a egzekute `black .`
black:
	black src

# 2 liy ki swiv yo vle di ke, le mwen tape `make lint` li a egzekute `pylint test.py`
lint:
	pylint src

# liy ki swiv lan vle di ke, le mwen tape `make blt` lap fe sa ki nesese pou ni `black` ni `lint` 
blt: black lint

# liy ki swiv lan vle di ke, black lint ak blt pa produi ankenn fichye
.PHONY: black lint blt
