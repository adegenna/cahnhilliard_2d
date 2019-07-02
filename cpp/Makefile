# this makefile is used to wrap the cpp solver in a python3.6 module
CPP_NONTHERMAL   = cahnhilliard
CPP_THERMAL      = cahnhilliard_thermal
CPP_PARAMS       = chparams
CPP_CHBASE       = cahnhilliard_nonlocal
DIR              = swig

all : ch_nonthermal

ch_nonthermal :
	swig -threads -c++ -python src/$(CPP_NONTHERMAL).i
	g++ -O3 -fopenmp -std=c++11 -c src/$(CPP_NONTHERMAL).cpp src/$(CPP_THERMAL).cpp src/$(CPP_PARAMS).cpp src/$(CPP_CHBASE).cpp src/$(CPP_NONTHERMAL)_wrap.cxx -I/usr/include/python3.6 -fPIC
	g++ -shared $(CPP_NONTHERMAL).o $(CPP_THERMAL).o $(CPP_PARAMS).o $(CPP_CHBASE).o $(CPP_NONTHERMAL)_wrap.o -o _$(CPP_NONTHERMAL).so
	mv $(CPP_NONTHERMAL).o $(CPP_THERMAL).o $(CPP_PARAMS).o $(CPP_CHBASE).o $(CPP_NONTHERMAL)_wrap.o _$(CPP_NONTHERMAL).so src/$(CPP_NONTHERMAL).py src/$(CPP_NONTHERMAL)_wrap.cxx $(DIR)/

clean :
	rm -rf $(DIR)/*.o $(DIR)/*~ $(DIR)/*.so $(DIR)/$(CPP_NONTHERMAL).py $(DIR)/$(CPP_NONTHERMAL)_wrap.cxx $(DIR)/__pycache__

