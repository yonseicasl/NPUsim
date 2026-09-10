############################################
# Developed by                             #
# Intelligent Computing Systems Lab (ICSL) #
# Yonsei University, Seoul, Korea          #
############################################

HDR=$(wildcard $(DIR)/scheduler/*.h)
HDR+=$(wildcard $(DIR)/utils/*.h)
HDR+=$(wildcard $(DIR)/components/*.h)

CCSRC=$(wildcard $(DIR)/scheduler/*.cc)
CCSRC+=$(wildcard $(DIR)/utils/*.cc)
CCSRC+=$(wildcard $(DIR)/components/*.cc)

# Build profile isolation (functional-sim plan G8): -DFUNCTIONAL changes data_t
# (utils/user-def.h), so timing and functional objects must never mix. Each profile
# compiles into its own obj/<profile>/ tree and the finished library is installed to
# $(LIB) only when its content differs, so toggling profiles relinks from cached
# objects instead of silently reusing stale ones.
PROFILE?=timing
OBJDIR=$(DIR)/library/obj/$(PROFILE)
CCOBJ=$(patsubst $(DIR)/%.cc,$(OBJDIR)/%.o,$(CCSRC))
LIBARTIFACT=$(OBJDIR)/$(notdir $(LIB))
# Compile-flag stamp: refreshed only when the flag string changes, so a flag change
# WITHIN a profile (e.g. DEBUG=1) still rebuilds every object of that profile.
FLAGSTAMP=$(OBJDIR)/.ccflags

default: install

.PHONY: install FORCE
install: $(LIBARTIFACT)
	@cmp -s $(LIBARTIFACT) $(LIB) 2>/dev/null || cp -f $(LIBARTIFACT) $(LIB)

$(LIBARTIFACT): $(CCOBJ)
	$(CC) -g -shared -Wl,-soname,$(LIB) -o $@ $(CCOBJ)

$(FLAGSTAMP): FORCE
	@mkdir -p $(OBJDIR)
	@printf '%s\n' '$(CCFLAG) $(STD)' | cmp -s - $@ 2>/dev/null || printf '%s\n' '$(CCFLAG) $(STD)' > $@

$(OBJDIR)/%.o: $(DIR)/%.cc $(HDR) $(FLAGSTAMP)
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAG) $(STD) -o $@ -c $<

clean:
	rm -rf $(DIR)/library/obj
	rm -f $(LIB)
