CC        = gcc
EXEC      = CSA
CFLAGS    = -Wall -Wextra -Wunused -mpopcnt
CDEBUG    = -g -ggdb
COPTIMIZE = -Wuninitialized -O9 -fomit-frame-pointer
CLIBS     = -lm

CSRCS     = $(wildcard *.c)
CHDRS     = $(wildcard *.h)
TXTS      = $(wildcard *.txt README* LICENSE*)
SCRIPTS   = $(wildcard Makefile* *.sh *.py)

ifeq ($(MAKECMDGOALS),debug)
	CFLAGS += $(CDEBUG)
else
	CFLAGS += $(COPTIMIZE)
endif

.PHONY: all clean pack

all: clean bin

debug: all

bin:
	@echo :: Compiling \"$(EXEC)\" ...
	$(CC) $(CFLAGS) $(CSRCS) -o $(EXEC) $(CLIBS)
	@echo :: Done

clean:
	@echo :: Cleaning up ...
	@rm -f $(EXEC) $(EXEC).zip

pack:
	@echo :: Packing files ...
	zip -9 -v $(EXEC).zip $(CSRCS) $(CHDRS) $(TXTS) $(SCRIPTS)
	@echo :: Done
