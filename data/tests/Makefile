#COMMENT=@echo making "'$@'" stem "'$*'" out of "'$+'" due to "'$?'"
#COMMENT=\#@ignore
COMMENT=
BBXS=bb2s
BBXP=bb2p
MPIOPT:=${shell mpirun --version | perl -ane '$$.==1 and $$v=pop @F and $$f=substr($$v,0,1) and $$f>=3 and print("--oversubscribe");'}
MPIRUN:=export TMPDIR=/tmp; mpirun $(MPIOPT) -np

TEST=echo -n test $* :
OK=echo ' ok'
FAIL=echo ' failed'; ./admitfailure

NEGATIVE=nodev depth nexist notype reserved noname declared divzero diff0d nomem nomem-4
POSITIVE=funcs matout record reduce reduce-4 \
	fhn0 fhn1 fhn2 fhn3 fhn3side fhn3side-16 \
	bkl3cont-0 bkl3cont-4 br0 hh52mrl \
	ballf ballfcont-0 ballfcont-4 balldcont-0 balldcont-4 
ALL=$(NEGATIVE) $(POSITIVE)
all: clean closeenough admitfailure ${foreach t,$(ALL),$t.test}
clean:
	-rm -f *.log *.out **.fullout *.sorted
	-rm -f bkl3cont-*.exact
	-rm -rf ball*.exact ball*.ctl ball*BAK ball*ct~

.PHONY: %.test

##########################
# Negative tests: scripts leading to exceptions; check for error messages in the log files.
%.test: %.log %.errmsg
	$(COMMENT)
	@$(TEST)
	@if grep "`cat $*.errmsg`" $< > /dev/null; then $(OK); else $(FAIL); fi
	@rm $*.log
%.log: %.bbs
	@-$(BBXS) $< -log $@ > /dev/null 2>/dev/null

%-4.log: %.bbs
	@-$(MPIRUN) 4 $(BBXP) $< > $@ 2>/dev/null

##########################
# Positive tests: output files compared to what they should be: either exact or approximate. 
%.test: %.out %.exact
	$(COMMENT)
	@$(TEST)
	@if cmp $+; then $(OK); else $(FAIL); fi
	@rm $*.out
%.test: closeenough %.out %.approx
	$(COMMENT)
	@$(TEST)
	@if ./$+; then $(OK); else $(FAIL); fi
	@rm $*.out
# Outputs of parallel jobs may come out in arbitrary order
%.test: closeenough %.sorted %.sortedapprox
	$(COMMENT)
	@$(TEST)
	@if ./$+; then $(OK); else $(FAIL); fi
	@rm $*.sorted
	@rm $*.out
%.sorted: %.out
	$(COMMENT)
	@sort $< > $@
# .PRECIOUS: fhn0.out
%.out: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(BBXS) $< $@ -log /dev/null > /dev/null
%.fullout: %.bbs
	$(COMMENT)
	$(BBXS) $< $@ -log $@.log -verbose | tee $@.stdout
%-4.out: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 4 $(BBXP) $< $@ -log /dev/null > /dev/null
%-4.fullout: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 4 $(BBXP) $+ $@ -log $@.log -verbose > $@.stdout
%-16.out: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 16 $(BBXP) $+ $@ -decomp 4x4x1 -log /dev/null > /dev/null
%-16.fullout: %.bbs
	$(COMMENT)
	@-rm -rf $@
	@$(MPIRUN) 16 $(BBXP) $+ $@ -decomp 4x4x1 -log /dev/null -verbose > $@.stdout
fhn1.out: fhn0.out fhn.par
fhn2.out: fhn1.out fhn.par
fhn3.out: fhn1.out fhn.par

bkl3.half bkl3.full: bkl3.bbs
	$(COMMENT)
	@$(BBXS) $< bkl3.half bkl3.full  -log /dev/null > /dev/null
.PRECIOUS: bkl3cont-%.exact
bkl3cont-%.exact: bkl3.full
	$(COMMENT)
	@cp -p $< $@
bkl3cont-0.out: bkl3cont.bbs bkl3.half
	$(COMMENT)
	@$(BBXS) $+ $@ -log /dev/null > /dev/null
bkl3cont-4.out: bkl3cont.bbs bkl3.half
	$(COMMENT)
	@$(MPIRUN) 4 $(BBXP) $+ $@ -log /dev/null > /dev/null

.PRECIOUS: ballfcont-%.exact # seq and mpi ctlpoints done in one go
ballfcont-0.exact ballfcont-4.exact ballf.out ballf.ctl: ballf.bbs ball.bbg
	$(COMMENT)
	@if [ -e ballf.ctl ]; then rm -f ballf.ctl; fi # do not start from control point
#	@$(BBXS) $+ ballf.ctl ballf.out -log /dev/null > /dev/null
	@$(BBXS) $+ ballf.ctl ballf.out -verbose -log $@.log  -debug log > /dev/null
	@cp -p ballf.out ballfcont-0.exact
	@cp -p ballf.out ballfcont-4.exact
ballfcont-0.out: ballfcont.bbs ball.bbg ballf.ctl # do start from control point
	$(COMMENT)
#	@$(BBXS) $+ $@ -log /dev/null > /dev/null
	@$(BBXS) $+ $@ -verbose -log $@.log  -debug log > /dev/null
ballfcont-4.out: ballfcont.bbs ball.bbg ballf.ctl # do start from control point
	$(COMMENT)
#	@$(MPIRUN) 4 $(BBXP) $+ $@ -log /dev/null > /dev/null
	@$(MPIRUN) 4 $(BBXP) $+ $@ -verbose -debug stdout > $@.log


.PRECIOUS: balldcont-%.exact # seq and mpi ctlpoints are incompatible and done separately
balldcont-0.exact balld-0.out balld-0.ctl: balld.bbs ball.bbg
	$(COMMENT)
	@if [ -e balld-0.ctl ]; then rm -rf balld-0.ctl; fi # do not start from control point
#	@$(BBXS) $+ balld-0.ctl balld-0.out -log /dev/null > /dev/null
	@$(BBXS) $+ balld-0.ctl balld-0.out -verbose -log $@.log  -debug log > /dev/null
	@cp -p balld-0.out balldcont-0.exact
balldcont-4.exact balld-4.out balld-4.ctl: balld.bbs ball.bbg
	$(COMMENT)
	@if [ -e balld-4.ctl ]; then rm -rf balld-4.ctl; fi # do not start from control point
#	@$(MPIRUN) 4 $(BBXP) $+ balld-4.ctl balld-4.out -log /dev/null > /dev/null
	@$(MPIRUN) 4 $(BBXP) $+ balld-4.ctl balld-4.out -verbose -debug stdout > $@.log
	@cp -p balld-4.out balldcont-4.exact
balldcont-0.out: balldcont.bbs ball.bbg balld-0.ctl # do start from control point
	$(COMMENT)
#	@$(BBXS) $+ $@ -log /dev/null > /dev/null
	@$(BBXS) $+ $@ -verbose -log $@.log  -debug log > /dev/null
balldcont-4.out: balldcont.bbs ball.bbg balld-4.ctl # do start from control point
	$(COMMENT)
#	@$(MPIRUN) 4 $(BBXP) $+ $@ -log /dev/null > /dev/null
	@$(MPIRUN) 4 $(BBXP) $+ $@ -verbose -debug stdout > $@.log

##########################
# Auxiliary utilities
closeenough: closeenough.c
	$(COMMENT)
	cc $< -o $@
admitfailure: admitfailure.c
	$(COMMENT)
	cc $< -o $@
