CMD   = python
FLAGS = install
DOC   = setup.py

all : setup

setup : $(DOC)
	$(CMD) $(DOC) install

script :
	$(CMD) -m cahnhilliard_2d.tests.driver cahnhilliard_2d/tests/input_driver.dat

clean :
	rm -rf build/ dist/ *egg-info

