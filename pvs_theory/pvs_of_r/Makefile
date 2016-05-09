include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif


target := DDpvs
SRC1   := DDpvs.c pvs_countpairs.c ../../utils/gridlink_pvs.c ../../utils/utils.c ../../io/ftread.c ../../io/io.c 
INCL   := pvs_countpairs.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_pvs.h \
		../../utils/function_precision.h ../../utils/cellarray_pvs.h ../../utils/avx_calls.h \
		../../utils/defs.h 

ifeq ($(OUTPUT_PGBAR), 1)
SRC1 += ../../utils/progressbar.c
INCL += ../../utils/progressbar.h 
endif
OBJS1  := $(SRC1:.c=.o)

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = pvs_countpairs.c ../../utils/gridlink_pvs.c ../../utils/utils.c
ifeq ($(OUTPUT_PGBAR), 1)
SRC2 += ../../utils/progressbar.c
endif

OBJS2  = $(SRC2:.c=.o)
libs: lib
lib:  ../../lib/libcountpairs_pvs.a | ../../lib

libcountpairs_pvs.a: $(OBJS2) $(INCL) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)

../../lib/libcountpairs_pvs.a: libcountpairs_pvs.a
	cp -p $< $@
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_pvs.h > ../../include/countpairs_pvs.h

.PHONY: clean clena celan install lib tests distclean

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $< $@

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_pvs.a ../../include/countpairs_pvs.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_pvs.a

clena: clean
celan: clean


#tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
#	$(MAKE) -C ../tests	xi_of_r
#
#../tests/tests_periodic:
#	$(MAKE) -C ../tests
#
#../tests/tests_nonperiodic:
#	$(MAKE) -C ../tests
#
