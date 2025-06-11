# If this file is modified you need to run:
#
# autoreconf --install
#
# To propagate the changes to the Makefile.in files.
#

GSL = gsl-1.16.extract
SUBDIRS = $(GSL)

bin_PROGRAMS = bb2p bb2s bb2d
bb2p_CPPFLAGS = -D MPI -D NO_TMPFILES ${common_cppflags} -I$(builddir)/$(GSL)
bb2p_LDADD = ${LDFLAGS}
bb2p_LDADD +=  $(GSL)/.libs/libgsl.a $(GSL)/cblas/.libs/libgslcblas.a

bb2s_CPPFLAGS = ${common_cppflags} -I$(builddir)/$(GSL)/
bb2s_LDFLAGS =  ${common_ldflags}
bb2s_LDADD = ${bb2p_LDADD}

bb2d_CPPFLAGS = ${common_cppflags} 
bb2d_LDFLAGS =  ${common_ldflags} 

common_sources = \
	beatbox.c beatbox.h \
	activation.c \
	ap.c \
	b02.c \
	bgi.c \
	bgi.h \
	bikt.c \
	bikt.h \
	hh52.c hh52_other.h hh52_tgate.h hh52_par.h hh52_ftab.h\
	hh52m.c hh52m_ik.h hh52m_ina.h\
	br77.c br77_fddt.h br77_ftab.h br77_fun.h br77_ngate.h br77_other.h br77_par.h br77_tgate.h  br77_const.h\
	byte.c \
	byte.h \
	byteout.c \
	bzoregon.c \
	bzr.c \
	cgle.c \
	cglel.c \
	cglep.c \
	channel.h \
	chrecruit.c \
	chrecruitpar.h \
	chrecruitvar.h \
	chtb.c \
	chtbpar.h \
	chtbvar.h \
	clock.c \
	crn98.c crn98_fun.h  crn98_ngate.h crn98_other.h crn98_par.h crn98_tgate.h crn98_const.h crn98_ftab.h crn98_fddt.h \
	crn.c crn_par.h crn_step.h crn_var.h \
	crnach.c crnach_par.h crnach_step.h crnach_var.h \
	ctlpoint.c ctlpoint.h \
	d_dt.c \
	device.c \
	device.h \
	devlist.h \
	diff.c diff1.h diff2.h diff3.h diff1d.h diff2d.h diff3d.h diff_printweights.h \
	displayName.h \
	dump.c \
	dump_types.h \
	dynamic.c \
	dynamic.h \
	ecg.c \
	elliptic.c \
	error.h \
	euler.c \
	extern.h \
	fhncub.c \
	fhnbkl.c \
	fhncubpar.c \
	fk.c \
	fk.h \
	fkmod.c \
	fkmod.h \
	flow2d.c \
	front.c \
	geometry.c \
	geometry.h \
	glerr.h \
	gpcell.c \
	gpcell.h \
        gpc_var.h \
	grad2d.c \
	ha.c \
	ha_cur.h \
	ha_par.h \
	ha_un.h \
	ha_var.h \
	hh.c \
	hh.h \
	init.c init.h \
	ionic.h ioniclist.h \
	kernik19.c kernik19_const.h kernik19_fddt.h kernik19_ftab.h kernik19_fun.h kernik19_ngate.h kernik19_other.h kernik19_par.h kernik19_tgate.h \
	k_.c \
	k_.h \
	k_calc.h \
	k_code.h \
	k_comp.c \
	k_comp.h \
	k_def.h \
	k_errors.h \
	k_exec.c \
	k_exec.h \
	k_expr.h \
	k_free.h \
	k_func.c \
	k_poincare.c \
	k_print.c \
	linear1.c\
	linear2.c\
	load.c \
	lrd.c \
	lrd_par.h \
	lrd_step.h \
	lrd_var.h \
	mchtb.c \
	mchtb.h \
	mchtbpar.h \
	mchtbvar.h \
	mouse.c  mouse_init.h  mouse_par.h  mouse_step.h  \
	mousevent.c  mousevent_init.h  mousevent_par.h  mousevent_step.h  \
	mpi_io_choice.h \
	ms.c \
	n62.c n62.h \
	neum2d.c \
	neum3d.c \
	noble.c noble.h noblecns.h nobleft1.h nobleft2.h nobleft3.h nobleftn.c nobleftn.h nobleglb.h nobleini.h nobleloc.h noblest2.c noblest2.h noblestu.c noblestu.h \
	noble2.c noble2.h \
	noblen.c noblen.h \
	noble_gbk.c \
	nograph.h \
	ord.c ord_par.h ord_step.h ord_init.h \
	p2clib.c p2c.h \
	pause.c \
	pipe.c pipe.h \
	poincare.c \
	ppcpash18.c ppcpash18_fun.h  ppcpash18_ngate.h ppcpash18_other.h ppcpash18_par.h ppcpash18_tgate.h ppcpash18_const.h ppcpash18_ftab.h ppcpash18_fddt.h \
	paci18.c paci18_fun.h  paci18_ngate.h paci18_other.h paci18_par.h paci18_tgate.h paci18_const.h paci18_ftab.h paci18_fddt.h \
	paci20.c paci20_fun.h  paci20_ngate.h paci20_other.h paci20_par.h paci20_tgate.h paci20_const.h paci20_ftab.h paci20_fddt.h \
	ppmout.c \
	pw_mult.c \
	qpp.c qpp.h \
	reduce.c \
	rhs.h rhslist.h \
	rk4.c \
	rushlarsen.c \
	sample.c \
	scr.h \
	screen.c screen.h \
	singz.c \
	sequence.c sequence.h \
	shell.c \
	state.c state.h \
	stop.c \
	system.c system.h \
	tb.c tbpar.h tbvar.h \
	ttnp.c ttnp_par.h ttnp_step.h ttnp_init.h ttnnp04.c ttnnp04_par.h ttnnp04_var.h ttnnp04_step.h \
	ttp06.c  ttp06_par.h ttp06_step.h ttp06_var.h\
	ttp2006.c ttp2006_const.h ttp2006_fddt.h ttp2006_ftab.h ttp2006_fun.h ttp2006_ngate.h ttp2006_other.h ttp2006_par.h ttp2006_tgate.h\
	trivial.c\
	update.c \
	vtkout2.c \
	windraw.c \
	windraw.h \
	zeeman.c \
	zfk.c

bb2p_SOURCES = \
	${common_sources} \
	decomp.c \
	decomp.h

bb2s_SOURCES = \
	${common_sources} \
	adi3d.c \
	bytein.c \
	ezpaint.c ezpaint.h \
	ezstep.c ezstep.h  ezstep1.h  ezstep2.h ezstep3.h ezstep4.h ezstep0d.h ezstep2d.h ezstep3d.h \
	ezview.c ezgraph3d.c ezmarching.c ezview.h ezpar.h ezfix.h ezvar.h ezini.h ezarr.h ezgraph3d.h ezmarching.h ezgenerators.h \
	imgout.c \
	k_clock.c \
	k_draw.c \
	k_imgout.c \
	k_paint.c \
	k_paintgl.c \
	k_plot.c \
	matout.c \
	ppmin.c \
	record.c \
	skrecord.c \
	screen_dump.c \
	torx.c \
	tory.c \
	torz.c

# Obsoleted in 2024, kept for backward compatibility
common_sources += \
	diff2dv.c \
	diff3dv.c \
	diffold.c \
	diffstep.c

bb2d_SOURCES = bbxd.c
