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
CCOBJ=$(CCSRC:.cc=.o)

default: $(LIB)

$(LIB): $(CCOBJ)
	$(CC) -g -shared -Wl,-soname,$@ -o $@ $(CCOBJ)

%.o: %.cc $(HDR)
	$(CC) $(CCFLAG) $(STD) -o $@ -c $<

clean:
	rm -f $(CCOBJ) $(LIB)
