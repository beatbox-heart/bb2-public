# Reference tests of the elliptic device: convergence of the solution to the exact
# solution of the PDE, convergence rate of the multigrid, agreement with plain
# Gauss-Seidel and with the diff device, the singular (no-flux) problem, various
# diffusion coefficients, isolated=1, and the independence of the number of MPI processes.
#
#   make            or  make all     all the tests: serial and mpi
#   make serial     the sequential tests: dirichlet neumann coef isolated
#   make mpi        the MPI tests, each with every number of processes in NPS
#   make dirichlet | neumann | coef | isolated     one group of tests, as a table
#   make u2d32      (any case name, see the lists below) just one case
#   make list       list the groups and their cases
#   make clean      remove the work directory, so that everything is run again
#
# Each case is run once and its row of results kept in $(WORK)/<case>.row, so making
# a group again just prints its table; a case is rerun if the scripts, or the Beatbox
# executable found in PATH (or given by BBXS, BBXP), are newer than its row.
# Variables (e.g. make BBXS=../../../../src/bb2s NPS="1 4" TLIM=600 mpi):
#   BBXS, BBXP   sequential and MPI Beatbox executables
#   NPS          numbers of processes for the MPI tests
#   TLIM         time limit of one run, seconds (a deadlock busy-waits rather than fails)
#   WORK         work directory
# For the meaning of the columns, see elltest.sh, and for the arguments of the cases,
# the comments in the scripts: ell_plane.bbs, ell_neum.bbs, ell_coef.par.

ifeq (,${shell echo $(MAKEFLAGS) | perl -ape 's/\-\-.*$$//' | grep n})
COMMENT=${warning making $@ out of $+ because of $? :}
else
COMMENT=${warning } ${shell echo making "\x1B[1;31m" $@ "\x1B[0m" out of $+  because of "\033[32m" $? "\033[0m" :}
endif
BBXS=bb2s
BBXP=bb2p
MPIOPT:=${shell mpirun --version 2>/dev/null | perl -ane '$$.==1 and $$v=pop @F and $$f=substr($$v,0,1) and $$f>=3 and print("--oversubscribe");'}
MPIRUN=mpirun $(MPIOPT) -np
NPS=1 2 4
TLIM=300
WORK=work
export BBXS BBXP MPIRUN NPS TLIM WORK
HELP=./elltest.sh
# The executables, as prerequisites of the results
BBXSEXE:=${shell command -v $(BBXS)}
BBXPEXE:=${shell command -v $(BBXP)}

##########################
# Part 1, Dirichlet problem, exact solution a plane wave (ell_plane.bbs):
# name = nx ny nz | th0 th1 al | rhs (0 analytic, 1 from diff) | upper_level | gs (1: also plain GS)
PLANE_u2d16       = 16  16  0 30  0  0 0 2 1
PLANE_u2d32       = 32  32  0 30  0  0 0 3 1
PLANE_u2d64       = 64  64  0 30  0  0 0 4 1
PLANE_u2d128      = 128 128 0 30  0  0 0 5 0
PLANE_v2d32       = 32  32  0 30 90  0 0 3 1
PLANE_v2d64       = 64  64  0 30 90  0 0 4 1
PLANE_v2d128      = 128 128 0 30 90  0 0 5 0
PLANE_u2d64d      = 64  64  0 30  0  0 1 4 0
PLANE_v2d64d      = 64  64  0 30 90  0 1 4 0
PLANE_r2d48x24    = 48  24  0 30 90  0 0 3 1
PLANE_r2d40x26    = 40  26  0 30 90  0 0 3 1
PLANE_u3d16       = 16  16 16 30  0 20 0 2 1
PLANE_u3d32       = 32  32 32 30  0 20 0 3 0
PLANE_v3d32       = 32  32 32 30 90 20 0 3 0
PLANE_r3d24x16x12 = 24  16 12 30 90 20 0 2 1
PLANE_r3d16x24x32 = 16  24 32 30 90 20 0 2 1
PLANE=u2d16 u2d32 u2d64 u2d128 v2d32 v2d64 v2d128 u2d64d v2d64d r2d48x24 r2d40x26 \
	u3d16 u3d32 v3d32 r3d24x16x12 r3d16x24x32

# Part 2, no-flux (singular) problem, pinned (ell_neum.bbs):
# name = nx ny nz | shape | th0 al | upper_level | reference cycles (0: no reference solve)
NEUM_n2d16        = 16  16  0 box   0  0 2 3000
NEUM_n2d32        = 32  32  0 box   0  0 3 8000
NEUM_n2d64        = 64  64  0 box   0  0 4 0
NEUM_n2d128       = 128 128 0 box   0  0 5 0
NEUM_o2d32        = 32  32  0 box  30  0 3 8000
NEUM_b2d32        = 32  32  0 ball 30  0 3 8000
NEUM_b2d64        = 64  64  0 ball 30  0 4 0
NEUM_b2d48x30     = 48  30  0 ball 30  0 3 8000
NEUM_n3d16        = 16  16 16 box   0  0 2 3000
NEUM_b3d16        = 16  16 16 ball 30 20 2 3000
NEUM_b3d32        = 32  32 32 ball 30 20 3 0
NEUM_b3d24x16x12  = 24  16 12 ball 30 20 2 3000
NEUM=n2d16 n2d32 n2d64 n2d128 o2d32 b2d32 b2d64 b2d48x30 n3d16 b3d16 b3d32 b3d24x16x12

# Part 3, various coefficients, Dirichlet, plane wave (ell_coef*.bbs, see ell_coef.par):
# name = kind | script | coef1 coef2 | nx ny nz | rhs | upper_level
VARD='D=~D0*(1+va*sin(bx*(x-1)*hx+by*(y-1)*hx))'
VARDL='Dl=~Dl0*(1+va*sin(bx*(x-1)*hx))'
COEF_c0n32  = 0 ell_coefng.bbs  D=D0 manypoint=0     32 32  0 0 3
COEF_c1n32  = 1 ell_coefiso.bbs D=D0 manypoint=0     32 32  0 0 3
COEF_c1n64  = 1 ell_coefiso.bbs D=D0 manypoint=0     64 64  0 0 4
COEF_c2n32  = 2 ell_coefiso.bbs D=D0 manypoint=1     32 32  0 0 3
COEF_c2n64  = 2 ell_coefiso.bbs D=D0 manypoint=1     64 64  0 0 4
COEF_c3n32  = 3 ell_coefiso.bbs $(VARD) manypoint=0  32 32  0 0 3
COEF_c3n64  = 3 ell_coefiso.bbs $(VARD) manypoint=0  64 64  0 0 4
COEF_c3n64d = 3 ell_coefiso.bbs $(VARD) manypoint=0  64 64  0 1 4
COEF_c4n32  = 4 ell_coef.bbs    $(VARDL) Dt=Dt0      32 32  0 0 3
COEF_c4n64  = 4 ell_coef.bbs    $(VARDL) Dt=Dt0      64 64  0 0 4
COEF_c4n64d = 4 ell_coef.bbs    $(VARDL) Dt=Dt0      64 64  0 1 4
COEF_c1n3d  = 1 ell_coefiso.bbs D=D0 manypoint=0     16 16 16 0 2
COEF_c2n3d  = 2 ell_coefiso.bbs D=D0 manypoint=1     16 16 16 0 2
COEF_c3n3d  = 3 ell_coefiso.bbs $(VARD) manypoint=0  16 16 16 0 2
COEF_c3n3dd = 3 ell_coefiso.bbs $(VARD) manypoint=0  16 16 16 1 2
COEF_c4n3d  = 4 ell_coef.bbs    $(VARDL) Dt=Dt0      16 16 16 0 2
COEF_c4n3dd = 4 ell_coef.bbs    $(VARDL) Dt=Dt0      16 16 16 1 2
COEF=c0n32 c1n32 c1n64 c2n32 c2n64 c3n32 c3n64 c3n64d c4n32 c4n64 c4n64d \
	c1n3d c2n3d c3n3d c3n3dd c4n3d c4n3dd

# Part 4, isolated=1 on x<=xb of the tissue vs the tissue cut at xb (ell_neum.bbs):
# name = nx ny nz | shape | th0 al | upper_level | xb
ISOL_i2d32  = 32 32  0 box   0  0 3 20
ISOL_i2d32o = 32 32  0 box  30  0 3 20
ISOL_ib2d64 = 64 64  0 ball 30  0 4 40
ISOL_ib3d16 = 16 16 16 ball 30 20 2 12
ISOL=i2d32 i2d32o ib2d64 ib3d16

# Part 5, MPI: each case with every number of processes in NPS
# name = kind (see elltest.sh) | nx ny nz | th0 th1 al | shape | upper_level
MPI_m2d32       = d 32 32  0 30 90  0 box  3
MPI_m2d40x26    = d 40 26  0 30 90  0 box  3
MPI_m3d16       = d 16 16 16 30  0 20 box  2
MPI_m3d24x16x12 = d 24 16 12 30 90 20 box  2
MPI_mn2d32      = n 32 32  0  0  0  0 box  3
MPI_mb2d64a     = n 64 64  0  0  0  0 ball 4
MPI_mb3d16      = n 16 16 16 30  0 20 ball 2
MPI_mb2d64o     = n 64 64  0 30  0  0 ball 4
MPI_mc3n64      = c3 64 64  0  0  0  0 box  4
MPI_mc4n3d      = c4 16 16 16  0  0  0 box  2
MPI_mi2d32o     = i 32 32  0 30  0  0 box  3
MPI=m2d32 m2d40x26 m3d16 m3d24x16x12 mn2d32 mb2d64a mb3d16 mb2d64o mc3n64 mc4n3d mi2d32o

##########################
.PHONY: all serial mpi dirichlet neumann coef isolated list clean \
	$(PLANE) $(NEUM) $(COEF) $(ISOL) $(MPI)
all: serial mpi
serial: dirichlet neumann coef isolated

# A group: run its cases (if needed), then print the table
SHOWROWS=for f in $^; do if [ -e $$f ]; then cat $$f; else cat $$f.failed; fi; done
dirichlet: $(PLANE:%=$(WORK)/%.row)
	@echo; $(HELP) header plane; $(SHOWROWS)
neumann: $(NEUM:%=$(WORK)/%.row)
	@echo; $(HELP) header neum; $(SHOWROWS)
coef: $(COEF:%=$(WORK)/%.row)
	@echo; $(HELP) header coef; $(SHOWROWS)
isolated: $(ISOL:%=$(WORK)/%.row)
	@echo; $(HELP) header isol; $(SHOWROWS)
mpi: $(MPI:%=$(WORK)/%.mpirow)
	@echo; $(HELP) header mpi; $(SHOWROWS)

# A single case: run it (if needed), then print its row with the header
$(PLANE): %: $(WORK)/%.row
	@$(HELP) header plane; $(SHOWROWS)
$(NEUM): %: $(WORK)/%.row
	@$(HELP) header neum; $(SHOWROWS)
$(COEF): %: $(WORK)/%.row
	@$(HELP) header coef; $(SHOWROWS)
$(ISOL): %: $(WORK)/%.row
	@$(HELP) header isol; $(SHOWROWS)
$(MPI): %: $(WORK)/%.mpirow
	@$(HELP) header mpi; $(SHOWROWS)

# Running a case. On failure, the row (saying FAILED and why) is kept as .failed,
# so the case is shown in the table and rerun next time.
COMMON=elltest.sh mkgeom.pl $(BBXSEXE)
$(PLANE:%=$(WORK)/%.row): $(WORK)/%.row: ell_plane.bbs $(COMMON) | $(WORK)
	$(COMMENT)
	@$(HELP) plane $* $(PLANE_$*) > $@ || mv $@ $@.failed
$(NEUM:%=$(WORK)/%.row): $(WORK)/%.row: ell_neum.bbs $(COMMON) | $(WORK)
	$(COMMENT)
	@$(HELP) neum $* $(NEUM_$*) > $@ || mv $@ $@.failed
$(COEF:%=$(WORK)/%.row): $(WORK)/%.row: ell_coef.bbs ell_coefiso.bbs ell_coefng.bbs ell_coef.par ell_coef.inc $(COMMON) | $(WORK)
	$(COMMENT)
	@$(HELP) coef $* $(COEF_$*) > $@ || mv $@ $@.failed
$(ISOL:%=$(WORK)/%.row): $(WORK)/%.row: ell_neum.bbs $(COMMON) | $(WORK)
	$(COMMENT)
	@$(HELP) isol $* $(ISOL_$*) > $@ || mv $@ $@.failed
$(MPI:%=$(WORK)/%.mpirow): $(WORK)/%.mpirow: ell_plane.bbs ell_neum.bbs ell_coef.bbs ell_coefiso.bbs \
		ell_coef.par ell_coef.inc elltest.sh mkgeom.pl $(BBXPEXE) | $(WORK)
	$(COMMENT)
	@$(HELP) mpi $(word 1,$(MPI_$*)) $* $(wordlist 2,99,$(MPI_$*)) > $@ || mv $@ $@.failed

$(WORK):
	@mkdir -p $@

list:
	@echo "dirichlet: $(PLANE)"; echo "neumann:   $(NEUM)"; echo "coef:      $(COEF)"
	@echo "isolated:  $(ISOL)"; echo "mpi:       $(MPI)"
clean:
	-rm -rf $(WORK)
