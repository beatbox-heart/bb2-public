# If this file is modified you need to run:
#
# autoreconf --install
#
# To propagate the changes to the Makefile.in files.
#

# Do not behave necessarily like a GNU package
AUTOMAKE_OPTIONS = foreign

# Directories that must be recoursed (expect to find a Makefile.am
# in each directory named in SUBDIRS).
SUBDIRS   = src

# Define some directories for use below.
datadir   = data
scriptdir = $(datadir)/scripts
workshopdir = $(scriptdir)/BeatboxWorkshop
hectordir = $(scriptdir)/parallel_Hector
paralldir = $(scriptdir)/parallel
seqdir    = $(scriptdir)/sequential
vizdir    = $(scriptdir)/viz3d
paramsdir = $(datadir)/parameters
geomdir   = $(datadir)/geometry
testdir   = $(datadir)/tests

# This is what will actually be contained in the distribution file.
EXTRA_DIST    = README \
		configure \
                $(SVNVER) \
		mydebug.sh \
		myrebuild.sh \
                $(docs) \
                $(dist_examples) 

# Document files.
docs          = \
                doc/bbg_howto.html \
                doc/beatbox.html \
                doc/beatbox_quickstart.html \
                doc/docstyle.css \
                doc/faq.html \
                doc/README  \
                doc/ToC.js \
		doc/Ohara_Rudy_2011.c \
		doc/Ohara_Rudy_2011.cellml \
		doc/Ohara_Rudy_2011.h \
		doc/beatbox.html \
		doc/ord.c \
		doc/ord0.bbs \
		doc/ord_init.h \
		doc/ord_par.h \
		doc/ord_step.h \
                doc/imgs/0000.png \
                doc/imgs/0001.png \
                doc/imgs/0002.png \
                doc/imgs/0003.png \
                doc/imgs/0004.png \
                doc/imgs/0005.png \
                doc/imgs/0006.png \
                doc/imgs/0007.png \
                doc/imgs/0008.png \
                doc/imgs/0009.png \
                doc/imgs/0010.png \
                doc/imgs/0011.png \
                doc/imgs/0012.png \
                doc/imgs/0013.png \
                doc/imgs/0014.png \
                doc/imgs/0015.png \
                doc/imgs/0016.png \
                doc/imgs/0017.png \
                doc/imgs/0018.png \
                doc/imgs/0019.png \
                doc/imgs/0020.png \
                doc/imgs/G1.png \
                doc/imgs/G2.png \
                doc/imgs/crn0_ap.png \
                doc/imgs/device-ring.png \
                doc/imgs/fhn0_ap.png \
                doc/imgs/fk0_ap.png \
                doc/imgs/fkmod_atrial.png \
                doc/imgs/fkmod_san.png \
                doc/imgs/lrd0_ap.png \
                doc/imgs/ord0_ap.png \
                doc/imgs/ttnnp04.png \
                doc/imgs/ttnp0_ap.png \
                doc/imgs/ttp06.png \
                doc/imgs/ttp06_ap.png \
                doc/imgs/mc_mat_diag.png \
                doc/imgs/str_ionic.png

# Example files
dist_examples= $(datadir)/README \
               $(dist_scripts_BeatboxWorkshop) \
               $(dist_scripts_hector) \
               $(dist_scripts_parallel) \
               $(dist_scripts_seq) \
               $(dist_scripts_viz3d) \
               $(dist_geom) \
               $(dist_params) \
               $(dist_tests)

# Geometry files
dist_geom    = \
               README

# Parameter files
dist_params  = \
               $(paramsdir)/br.par \
               $(paramsdir)/fhn.par \
               $(paramsdir)/lrd.par

# Beatbox Workshop scripts
dist_scripts_BeatboxWorkshop = \
               $(workshopdir)/CRN_model/pd_crn2.bbs \
               $(workshopdir)/CRN_model/pd_crn0.bbs \
               $(workshopdir)/CRN_model/pd_crn1.bbs \
               $(workshopdir)/CRN_model/pd_crn1.rec \
               $(workshopdir)/FitzHughNagumo_model/fhn3.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn3_PositiveTension.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn3_NegativeTension.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn_ffr_iso_xz.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn_crossFieldStim_ffr_slice.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn2.bbs \
               $(workshopdir)/FitzHughNagumo_model/surface.pl \
               $(workshopdir)/FitzHughNagumo_model/fhn_spiral_ffr_slice_aniso.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn.par \
               $(workshopdir)/FitzHughNagumo_model/fhn_spiral_ffr_slice.bbs \
               $(workshopdir)/FitzHughNagumo_model/inside.pl \
               $(workshopdir)/FitzHughNagumo_model/fhn_ffr_xz.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn_ffr.bbs \
               $(workshopdir)/FitzHughNagumo_model/fhn0.bbs \
               $(workshopdir)/FitzHughNagumo_model/view_ffr.pl \
               $(workshopdir)/FitzHughNagumo_model/ffr.bbg \
               $(workshopdir)/FitzHughNagumo_model/fhn_ffr_iso.bbs \
               $(workshopdir)/FitzHughNagumo_model/view.pl \
               $(workshopdir)/FitzHughNagumo_model/fhn1.bbs \
               $(workshopdir)/FitzHughNagumo_model/ffr_slice.bbg \
               $(workshopdir)/parallel_Hector/ffr.bbg\
               $(workshopdir)/parallel_Hector/fhn1_NegativeTension.rec \
               $(workshopdir)/parallel_Hector/fhn1_PositiveTension.rec \
               $(workshopdir)/parallel_Hector/fhn3_NegativeTension.bbs \
               $(workshopdir)/parallel_Hector/fhn_ffr_xz.bbs \
               $(workshopdir)/parallel_Hector/fhn_ffr_xz_scaling.bbs \
               $(workshopdir)/parallel_Hector/fhn3_NegativeTension.pbs \
               $(workshopdir)/parallel_Hector/fhn_ffr_xz.pbs

# HECToR example scripts
dist_scripts_hector = \
               $(hectordir)/box3D/bigbox_crn.bbs \
               $(hectordir)/box3D/bigbox_crn.sh  \
               $(hectordir)/box3D/bigbox_fhn.bbs \
               $(hectordir)/box3D/bigbox_fhn.sh  \
               $(hectordir)/box3D/bigbox_lrd.bbs \
               $(hectordir)/box3D/bigbox_lrd.sh  \
               $(hectordir)/phasedbn_parallel/pd_crn0.rec \
               $(hectordir)/phasedbn_parallel/pd_crn1.rec \
               $(hectordir)/phasedbn_parallel/pd_crn2.bbs \
               $(hectordir)/phasedbn_parallel/submit_pdcrn2.sh

# Parallel  example scripts
dist_scripts_parallel    = \
                $(paralldir)/Barkley/Makefile \
                $(paralldir)/Barkley/bkl1.rec \
                $(paralldir)/Barkley/scroll.bbs \
                $(paralldir)/Barkley/view.pl \
                $(paralldir)/CRN_model/crn_ffr.bbs \
                $(paralldir)/CRN_model/pd_crn1.rec \
                $(paralldir)/FitzHughNagumo_model/ffr.bbg \
                $(paralldir)/FitzHughNagumo_model/fhn1_NegativeTension.rec \
                $(paralldir)/FitzHughNagumo_model/fhn3_NegativeTension.bbs \
                $(paralldir)/FitzHughNagumo_model/fhn3_NegativeTension.pbs \
                $(paralldir)/FitzHughNagumo_model/fhn_ffr_xz.bbs \
                $(paralldir)/FitzHughNagumo_model/fhn_ffr_xz.pbs \
                $(paralldir)/FitzHughNagumo_model/fhn_ffr_xz_scaling.bbs \
                $(paralldir)/FitzHughNagumo_model/fhncub.par \
                $(paralldir)/FitzHughNagumo_model/sample.bbs \
                $(paralldir)/FitzHughNagumo_model/sample.rec

# Sequential example scripts
dist_scripts_seq    = \
                $(seqdir)/Barkley_model/bkl1.rec \
                $(seqdir)/Barkley_model/ez.bbs \
                $(seqdir)/Barkley_model/ezscroll.bbs \
                $(seqdir)/Barkley_model/ezspiral.bbs \
                $(seqdir)/Barkley_model/implicit.bbs \
                $(seqdir)/Barkley_model/bkl1.rec \
                $(seqdir)/Barkley_model/rabbit.bbs \
                $(seqdir)/Barkley_model/rabbit3.bbg \
		$(seqdir)/BeelerReuter_model/otg.bbs \
		$(seqdir)/BeelerReuter_model/tgo.bbs \
		$(seqdir)/BeelerReuter_model/tog.bbs \
                $(seqdir)/CRN_model/compare1D.bbs \
                $(seqdir)/CRN_model/crn0.bbs \
                $(seqdir)/CRN_model/crn0_ap.bbs \
                $(seqdir)/CRN_model/crn0_rest.bbs \
                $(seqdir)/CRN_model/crn1.bbs \
                $(seqdir)/CRN_model/pd_crn0.bbs \
                $(seqdir)/CRN_model/pd_crn1.bbs \
                $(seqdir)/CRN_model/pd_crn2.bbs \
                $(seqdir)/CRN_model/pd_runcrn0.sh \
                $(seqdir)/CRN_model/pd_runcrn1.sh \
                $(seqdir)/CRN_model/pd_runcrn2.sh \
                $(seqdir)/CRN_model/restitutioncrn.sh \
                $(seqdir)/CRN_model/runcrn0.sh \
                $(seqdir)/CRN_model/runcrn1.sh \
                $(seqdir)/CRN_model/compare1D.bbs \
                $(seqdir)/FentonKarma_model/Makefile \
                $(seqdir)/FentonKarma_model/fk0_ap.bbs \
                $(seqdir)/FentonKarma_model/fk0_atrial.bbs \
                $(seqdir)/FentonKarma_model/fk0_san.bbs \
                $(seqdir)/FentonKarma_model/fk2_san.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn.par \
                $(seqdir)/FitzHughNagumo_model/fhn0.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn0_ap.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn1.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn1_load.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn2.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn3.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn_crossFieldStim_ffr_slice.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn_ffr.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn_spiral_ffr_slice.bbs \
                $(seqdir)/FitzHughNagumo_model/fhn_spiral_ffr_slice_aniso.bbs \
                $(seqdir)/FitzHughNagumo_model/RunTest.sh \
		$(seqdir)/FitzHughNagumo_model/singz.bbs \
                $(seqdir)/LRD_model/RunTest.sh \
                $(seqdir)/LRD_model/lrd.par \
                $(seqdir)/LRD_model/lrd0_ap.bbs \
                $(seqdir)/LRD_model/LRD_3d_box.bbs \
                $(seqdir)/LRD_model/LRD_ffr.bbs \
                $(seqdir)/LRD_model/LRD_ffr_slice.bbs \
                $(seqdir)/MitchellSchaeffer_model/RunTest.sh \
                $(seqdir)/MitchellSchaeffer_model/ms.par \
                $(seqdir)/MitchellSchaeffer_model/ms0.bbs \
                $(seqdir)/MitchellSchaeffer_model/ms1.bbs \
                $(seqdir)/MitchellSchaeffer_model/ms2.bbs \
		$(seqdir)/RK4/Makefile \
		$(seqdir)/RK4/rk4test-control.eps \
		$(seqdir)/RK4/rk4test-control.txt \
		$(seqdir)/RK4/rk4test.bbs \
		$(seqdir)/RK4/rk4test.eps \
		$(seqdir)/RK4/rk4test.gnu \
		$(seqdir)/RK4/rk4test.txt \
                $(seqdir)/TenTusscherPanfilov_modelS/ttnnp04.png \
                $(seqdir)/TenTusscherPanfilov_modelS/ttnnp04_ap.bbs \
                $(seqdir)/TenTusscherPanfilov_modelS/ttnp0.bbs \
                $(seqdir)/TenTusscherPanfilov_modelS/ttnp0_ap.bbs \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06.bbs \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06_ap.bbs \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06 \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/Makefile \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/benchmark.bbs \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/colormap.c \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/dat2tex.pl \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/domain.pl \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/graph.gnu \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/report-original.pdf \
                $(seqdir)/TenTusscherPanfilov_modelS/ttp06/report.tex \
                $(seqdir)/ZFK_model/Makefile \
                $(seqdir)/ZFK_model/square.pl \
                $(seqdir)/ZFK_model/square_20.bbg \
                $(seqdir)/ZFK_model/zfkmon2.bbs \
                $(seqdir)/ZFK_model/zfkbid2.bbs \
                $(seqdir)/ord0_ap.bbs \
                $(seqdir)/b02/Makefile \
                $(seqdir)/b02/fig5.bbs \
                $(seqdir)/b02/fig5a.pgm \
                $(seqdir)/b02/fig5b.pgm \
                $(seqdir)/b02/fig5a.png \
                $(seqdir)/b02/fig5b.png \
		$(seqdir)/hh52/hh52m.bbs \
		$(seqdir)/hh52/Makefile

dist_tests = \
		$(testdir)/Makefile \
		$(testdir)/README \
		$(testdir)/admitfailure.c \
		$(testdir)/ball.bbg \
		$(testdir)/ballf.bbs \
		$(testdir)/ballfcont.bbs \
		$(testdir)/balld.bbs \
		$(testdir)/balldcont.bbs \
		$(testdir)/ballf.approx \
		$(testdir)/bkl1.rec \
		$(testdir)/bkl3.bbs \
		$(testdir)/bkl3cont.bbs \
		$(testdir)/br0.bbs \
		$(testdir)/br0.approx \
		$(testdir)/closeenough.c \
		$(testdir)/declared.bbs \
		$(testdir)/declared.errmsg \
		$(testdir)/depth.bbs \
		$(testdir)/depth.errmsg \
		$(testdir)/diff0d.bbs \
		$(testdir)/diff0d.errmsg \
		$(testdir)/divzero.bbs \
		$(testdir)/divzero.errmsg \
		$(testdir)/fhn.par \
		$(testdir)/fhn0.approx \
		$(testdir)/fhn0.bbs \
		$(testdir)/fhn1.approx \
		$(testdir)/fhn1.bbs \
		$(testdir)/fhn2.approx \
		$(testdir)/fhn2.bbs \
		$(testdir)/fhn3.approx \
		$(testdir)/fhn3.bbs \
		$(testdir)/fhn3side-16.sortedapprox \
		$(testdir)/fhn3side.approx \
		$(testdir)/fhn3side.bbs \
		$(testdir)/funcs.approx \
		$(testdir)/funcs.bbg \
		$(testdir)/funcs.bbs \
		$(testdir)/hh52mrl.approx \
		$(testdir)/hh52mrl.bbs \
		$(testdir)/matout.bbs \
		$(testdir)/matout.exact \
		$(testdir)/nexist.bbs \
		$(testdir)/nexist.errmsg \
		$(testdir)/nodev.bbs \
		$(testdir)/nodev.errmsg \
		$(testdir)/nomem-4.errmsg \
		$(testdir)/nomem.bbs \
		$(testdir)/nomem.errmsg \
		$(testdir)/noname.bbs \
		$(testdir)/noname.errmsg \
		$(testdir)/notype.bbs \
		$(testdir)/notype.errmsg \
		$(testdir)/record.bbs \
		$(testdir)/record.exact \
		$(testdir)/reduce-4.approx \
		$(testdir)/reduce.approx \
		$(testdir)/reduce.bbs \
		$(testdir)/reserved.bbs \
		$(testdir)/reserved.errmsg 


# viz3d examples
dist_scripts_viz3d  = \
                $(vizdir)/ppmz.task12
