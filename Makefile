CMD   = python
FLAGS = install
DOC   = setup.py

all : setup

setup : $(DOC)
	$(CMD) $(DOC) install

script :
	$(CMD) -m cahnhilliard_2d.tests.driver cahnhilliard_2d/tests/input_driver.dat

cpp :
	g++ -std=c++11 ./cahnhilliard_2d/src/cahn_hilliard.cpp

clean :
	rm -rf build/ dist/ *egg-info

