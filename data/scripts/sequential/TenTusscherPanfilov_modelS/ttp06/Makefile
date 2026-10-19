
##################
# The calculations
BEATBOX=bb2s

# Create the geometry
%.bbg: domain.pl
	$(COMMENT)
	$< $* > $@

# All space steps
ALLHX=0.5 0.2 0.1
# All time steps
ALLHT=0.05 0.01 0.005
# All their combinations
ALLDAT=${foreach x,$(ALLHX),${foreach t,$(ALLHT),$x_$t.dat}}
.PHONY: alldat
alldat: $(ALLDAT)

# The runs
.SECONDEXPANSION:
%.dat: benchmark.bbs $${firstword $${subst _, , $$*}}.bbg
	$(BEATBOX) $+ ${subst _, , $*} $@ -verbose -log $*.log
#	$(BEATBOX) $+ ${subst _, , $*} $@ -verbose -log $*.log -nograph

####################
# Table in the Supplement
table.tex: dat2tex.pl $(wildcard *_*.dat)
	$+ | tee $@

####################
# Figures
%.pdf: %.eps
	epstopdf $< --outfile=$@

# Figure 2
RED=0.1_0.005.dat.dir/actdiag.txt
GREEN=0.2_0.005.dat.dir/actdiag.txt
BLUE=0.5_0.005.dat.dir/actdiag.txt
e:
	echo e > $@
graph.dat: $(RED) e $(GREEN) e $(BLUE) e
	cat $+ > $@
graph.eps: graph.gnu graph.dat
	cat $+ | gnuplot > $@

# Figure 3
%.pgm: %.dat.dir/activation.ppm
	ppmtorgb3 < $<
	mv noname.red $@
	rm noname.*
# Discard last row and column to make them proportional
# (this is due to points vs intervals mismatch)
# and scale up to the same size
upper.pgm: 0.5_0.05.pgm
	pamcut -right=39 -top=1 $< | pamstretch 10 > $@
lower.pgm: 0.1_0.005.pgm
	pamcut -right=199 -top=1 $< | pamstretch 2 > $@
colormap: colormap.c
	cc $< -o $@
colormap.ppm: colormap
	$+ > $@
%.eps: colormap.ppm %.pgm
	pgmtoppm -map $+ | pnmtops -noturn -imagewidth=2 >  $@


####################
# Generate the report
report.pdf: report.tex table.tex graph.pdf upper.pdf lower.pdf
	pdflatex $<

