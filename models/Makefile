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
CCOBJ=$(CCSRC:.cc=.o)

default: $(EXE)

$(EXE): $(CCOBJ) $(LIB)
	$(LC) $(LDFLAG) -o $@ $(CCOBJ) $(LIBFLAG)

%.o: %.cc $(HDR)
	$(CC) $(CCFLAG) $(STD) -o $@ -c $<

clean:
	rm -f $(CCOBJ) $(CUOBJ) $(EXE)
