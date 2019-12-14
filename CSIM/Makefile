include config.mk
lflags=-L. `gsl-config --libs`
iflags=`gsl-config --cflags`

# Lists of files to be built
objs=common.o common_tasks.o networks.o depl_sim_nl.o
src=$(patsubst %.o,%.cc,$(objs))
execs=constant_n_v2

all:
	$(MAKE) executables

executables: $(execs)

include Makefile.dep

depend:
	$(cxx) $(cflags) $(iflags) -MM $(src) >Makefile.dep

clean:
	rm -rf $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

libdepl.a: $(objs)
	rm -f libdepl.a
	ar rs libdepl.a $^

constant_n_v2: constant_n.cc libdepl.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -ldepl -lgsl -lgslcblas

.PHONY: clean depend executables
