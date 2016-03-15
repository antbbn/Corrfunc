include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target=DDrppi_pvs

SRC1= DDrppi_pvs.c pvs_countpairs_rp_pi.c ../../utils/utils.c ../../utils/gridlink_pvs.c ../../utils/progressbar.c ../../io/io.c ../../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = pvs_countpairs_rp_pi.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_pvs.h \
		../../utils/function_precision.h ../../utils/cellarray_pvs.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h 

all: $(target) ../../theory.options ../../common.mk Makefile 

DDrppi_pvs: $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

#OBJS2= wprp.o ../../utils/utils.o 
#wprp: $(OBJS2) ../../theory.options ../../common.mk Makefile
#	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(OBJS2) $(CLINK)

libs: lib
lib: libpvs_countpairs_rp_pi.a | ../../lib

OBJS3=pvs_countpairs_rp_pi.o ../../utils/gridlink_pvs.o ../../utils/utils.o ../../utils/progressbar.o 
libcountpairs_rp_pi.a: $(OBJS3) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS3)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" pvs_countpairs_rp_pi.h > ../../include/pvs_countpairs_rp_pi.h

../../lib ../../bin ../../include:
	mkdir -p $@

install: $(target) | ../../bin 
	cp -p $(target) ../../bin/

#.PHONY: clean clena celan install lib tests
.PHONY: clean clena celan install lib

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) $(OBJS3) libpvs_countpairs_rp_pi.a

clena: clean

celan: celan

distclean:clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libpvs_countpairs_rp_pi.a ../../include/countpairs_rp_pi.h
	cd ../../bin && $(RM) $(target)

#tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
#	$(MAKE) -C ../tests xi_rp_pi
#
#../tests/tests_periodic:
#	$(MAKE) -C ../tests
#
#../tests/tests_nonperiodic:
#	$(MAKE) -C ../tests


