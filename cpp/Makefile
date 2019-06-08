# this makefile is used to wrap the cpp solver in a python3.6 module
CPP   = cahnhilliard
DIR   = swig

all : python_module

python_module :
	swig -threads -c++ -python src/$(CPP).i
	g++ -O3 -fopenmp -std=c++11 -c src/$(CPP).cpp src/$(CPP)_wrap.cxx -I/usr/include/python3.6 -fPIC
	g++ -shared $(CPP).o $(CPP)_wrap.o -o _$(CPP).so
	mv $(CPP).o $(CPP)_wrap.o _$(CPP).so src/$(CPP).py src/$(CPP)_wrap.cxx $(DIR)/

clean :
	rm -rf $(DIR)/*.o $(DIR)/*~ $(DIR)/*.so $(DIR)/$(CPP).py $(DIR)/$(CPP)_wrap.cxx $(DIR)/__pycache__

