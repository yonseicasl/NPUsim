############################################
# Developed by                             #
# Intelligent Computing Systems Lab (ICSL) #
# Yonsei University, Seoul, Korea          #
############################################

HDR=$(wildcard *.h)
# main.cc includes the library headers (npu.h -> stats.h, etc.) through the -I paths,
# so the model objects must be rebuilt whenever any of those change. Without this the
# model links against a stale object whose npu_t/stats_t layout disagrees with the freshly
# built libnpusim.so (ABI skew -> heap corruption at exit).
HDR+=$(wildcard $(DIR)/scheduler/*.h)
HDR+=$(wildcard $(DIR)/utils/*.h)
HDR+=$(wildcard $(DIR)/components/*.h)
CCSRC=$(wildcard *.cc)

# Build profile isolation (functional-sim plan G8): objects and the linked executable
# live in obj/<profile>/ and the executable is installed to ./$(EXE) only when its
# content differs -- see library/Makefile for the rationale.
PROFILE?=timing
OBJDIR=obj/$(PROFILE)
CCOBJ=$(patsubst %.cc,$(OBJDIR)/%.o,$(CCSRC))
EXEARTIFACT=$(OBJDIR)/$(EXE)
FLAGSTAMP=$(OBJDIR)/.ccflags

default: install

.PHONY: install FORCE
install: $(EXEARTIFACT)
	@cmp -s $(EXEARTIFACT) $(EXE) 2>/dev/null || cp -f $(EXEARTIFACT) $(EXE)

$(EXEARTIFACT): $(CCOBJ) $(LIB)
	$(LC) $(LDFLAG) -o $@ $(CCOBJ) $(LIBFLAG)

$(FLAGSTAMP): FORCE
	@mkdir -p $(OBJDIR)
	@printf '%s\n' '$(CCFLAG) $(STD)' | cmp -s - $@ 2>/dev/null || printf '%s\n' '$(CCFLAG) $(STD)' > $@

$(OBJDIR)/%.o: %.cc $(HDR) $(FLAGSTAMP)
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAG) $(STD) -o $@ -c $<

clean:
	rm -rf obj
	rm -f $(EXE)
