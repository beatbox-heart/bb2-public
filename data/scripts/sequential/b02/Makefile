# Make panels (a) and (b) of fig. 5 in 
# V.N. Biktashev, Phys. Rev. Lett., 89(16): 168102, 2002
# in png format, using BeatBox and netpbm
BBX=bb2s

all: fig5a.png fig5b.png

# Script fig5.bbs takes two command line arguments: the c value
# for the initial condition, and the name of the ppm output file. 
# As the figure caption says, c=0.624 for the left panel
fig5a.pgm: fig5.bbs
	$(BBX) $+ 0.624 $@

# and c=0.625 for the right panel
fig5b.pgm: fig5.bbs
	$(BBX) $+ 0.625 $@

# Conversion from ppm to png include:
# * contrast normalization to engage the full range of
#   gray values from black to white,
# * top/bottom flip
#   to make earlies time values at the bottom,
# * vertical stretching to make the 1:1 aspect
#   ratio as in the paper
%.png: %.pgm
	pnmnorm $< | pnmflip -tb | pnmscale -xscale=1 -yscale=2.5 | pnmtopng > $@
