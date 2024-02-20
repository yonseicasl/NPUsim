############################################
# Developed by                             #
# Intelligent Computing Systems Lab (ICSL) #
# Yonsei University, Seoul, Korea          #
############################################

HDR=$(wildcard *.h)
CCSRC=$(wildcard *.cc)
CCOBJ=$(CCSRC:.cc=.o)

default: $(EXE)

$(EXE): $(CCOBJ) $(LIB)
	$(LC) $(LDFLAG) -o $@ $(CCOBJ) $(LIBFLAG)

%.o: %.cc $(HDR)
	$(CC) $(CCFLAG) $(STD) -o $@ -c $<

clean:
	rm -f $(CCOBJ) $(CUOBJ) $(EXE)
