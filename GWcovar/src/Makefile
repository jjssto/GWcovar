all: covar.so 

covar.so:
	$(R_HOME)/bin/R CMD SHLIB covar.c wendland.c -lm -lgsl -fPIC

clean:
	rm wendland.o covar.o covar.so

