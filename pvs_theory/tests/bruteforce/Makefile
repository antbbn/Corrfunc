ROOTDIR=../../..

include $(ROOTDIR)/theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include $(ROOTDIR)/common.mk
endif

INCLUDE=-I$(ROOTDIR)/io -I$(ROOTDIR)/utils


target=DDrppi_pvs_bruteforce DDpvs_bruteforce

SRC1= DDrppi_pvs_bruteforce.c pvs_countpairs_rp_pi_bruteforce.c $(ROOTDIR)/utils/utils.c $(ROOTDIR)/utils/progressbar.c $(ROOTDIR)/io/io.c $(ROOTDIR)/io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = pvs_countpairs_bruteforce.h pvs_countpairs_rp_pi_bruteforce.h $(ROOTDIR)/io/ftread.h $(ROOTDIR)/io/io.h $(ROOTDIR)/utils/utils.h \
		$(ROOTDIR)/utils/function_precision.h  $(ROOTDIR)/utils/avx_calls.h \
		$(ROOTDIR)/utils/defs.h $(ROOTDIR)/utils/progressbar.h 

SRC2= DDpvs_bruteforce.c pvs_countpairs_bruteforce.c $(ROOTDIR)/utils/utils.c $(ROOTDIR)/utils/progressbar.c $(ROOTDIR)/io/io.c $(ROOTDIR)/io/ftread.c 
OBJS2  = $(SRC2:.c=.o)

all: $(target) $(ROOTDIR)/theory.options $(ROOTDIR)/common.mk Makefile 

DDrppi_pvs_bruteforce: $(OBJS1) $(INCL) $(ROOTDIR)/theory.options $(ROOTDIR)/common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@

DDpvs_bruteforce: $(OBJS2) $(INCL2) $(ROOTDIR)/theory.options $(ROOTDIR)/common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) -o $@

%.o: %.c $(INCL) $(ROOTDIR)/theory.options $(ROOTDIR)/common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan install lib

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) $(OBJS3) 

clena: clean

celan: celan

