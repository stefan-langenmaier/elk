
#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
MAKE = make
AR = ar
F90 = f90
F90_OPTS =
F77 = f90
F77_OPTS =
LIB_SYS =
LIB_LPK = lapack.a blas.a
LIB_FFT = fftlib.a
LIB_XC =

include ../make.inc

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

# modules
SRC_libxc = libxcifc_stub.f90
include libxc.inc
SRC_modules = modmain.f90 $(SRC_libxc) modxcifc.f90 modldapu.f90 modrdm.f90 \
 modphonon.f90 modtest.f90

# main program
SRC_main = main.f90

# subroutines and functions
SRC_routines = \
 autoradmt.f90 findprim.f90 findsym.f90 findsymcrys.f90 genppts.f90 rfinp.f90 \
 rvfcross.f90 seceqn.f90 symrf.f90 symrfir.f90 wavefmt.f90 force.f90 \
 elfplot.f90 rhonorm.f90 energy.f90 gensdmat.f90 \
 zpotclmt.f90 writegeom.f90 nfftifc.f90 zfftifc.f90 allatoms.f90 gridsize.f90 \
 poteff.f90 genrmesh.f90 readfermi.f90 potcoul.f90 gensfacgp.f90 \
 checkmt.f90 zfinp.f90 match.f90 forcek.f90 writeefg.f90 \
 findsymlat.f90 genffacg.f90 mixerifc.f90 mixpack.f90 \
 genlofr.f90 atom.f90 writefermi.f90 writekpts.f90 fsmfield.f90 mossbauer.f90 \
 occupy.f90 writelinen.f90 writeinfo.f90 readinput.f90 readspecies.f90 \
 charge.f90 moment.f90 writesym.f90 genidxlo.f90 gencore.f90 addrhocr.f90 \
 gengvec.f90 genshtmat.f90 plot1d.f90 plot2d.f90 plot3d.f90 updatpos.f90 \
 writeiad.f90 symvect.f90 vecplot.f90 genylmg.f90 linengy.f90 init0.f90 \
 init1.f90 gengpvec.f90 genveffig.f90 \
 gencfun.f90 genapwfr.f90 seceqnfv.f90 seceqnsv.f90 seceqnss.f90 getngkmax.f90 \
 dos.f90 rhoinit.f90 potplot.f90 writestate.f90 potxc.f90 zpotcoul.f90 \
 gndstate.f90 rhomagk.f90 readstate.f90 bandstr.f90 writeeval.f90 rhoplot.f90 \
 gwf2val.f90 gwf2cr.f90 rfarray.f90 fermisurf.f90 symrvfir.f90 symrvf.f90 \
 symrfmt.f90 wfplot.f90 hmlrad.f90 olprad.f90 olpistl.f90 olpaa.f90 olpalo.f90 \
 olplolo.f90 hmlistl.f90 hmlaa.f90 hmlalo.f90 hmllolo.f90 init2.f90 \
 writelsj.f90 wavefcr.f90 gradrf.f90 dbxcplot.f90 effmass.f90 projsbf.f90 \
 seceqnhf.f90 writelat.f90 energynn.f90 findigp0.f90 writeengy.f90 \
 writechg.f90 writeforce.f90 genjlgpr.f90 findkpt.f90 putevalfv.f90 \
 getevalfv.f90 putevalsv.f90 getevalsv.f90 putoccsv.f90 getoccsv.f90 \
 putevecfv.f90 getevecfv.f90 putevecsv.f90 getevecsv.f90 delevec.f90 \
 genvmatk.f90 genvmat.f90 genwfsv.f90 gendmat.f90 findsymsite.f90 geomplot.f90 \
 hartfock.f90 genkinmat.f90 energykncr.f90 rfpack.f90 rfint.f90 \
 writefsm.f90 genbmatk.f90 seceqnit.f90 gdft.f90 genexpiqr.f90 writeexpiqr.f90 \
 elnes.f90 symdmat.f90 genlmirep.f90 rhomagsh.f90 testcheck.f90 mixander.f90 \
 fermisurfbxsf.f90 checkfsm.f90 genzqrmt.f90 genbeffmt.f90 writesf.f90 \
 ggamt_sp_1.f90 ggair_sp_1.f90 ggamt_1.f90 ggair_1.f90 genspecies.f90 \
 ggamt_2a.f90 ggair_2a.f90 ggamt_2b.f90 ggair_2b.f90 ggamt_sp_2a.f90 \
 ggair_sp_2a.f90 ggamt_sp_2b.f90 ggair_sp_2b.f90 writeqpts.f90 potefield.f90 \
 findswidth.f90

SRC_phonon = \
 phonon.f90 dyntask.f90 phcell.f90 phdisp.f90 readdyn.f90 dynqtor.f90 \
 dynrtoq.f90 dynsymapp.f90 dynsym.f90 dyndiag.f90 phdos.f90 sumrule.f90 \
 writephn.f90 phveff.f90 phdveff.f90 writedveff.f90 readdveff.f90 \
 phfext.f90 epcouple.f90 genepmat.f90 writegamma.f90 readgamma.f90 \
 phlwidth.f90 writelambda.f90 alpha2f.f90 phdelete.f90

SRC_optics = \
 genpmat.f90 writepmat.f90 moke.f90 dielectric.f90 dielectricq.f90

SRC_advanced = \
 vnlrho.f90 vnlrhomt.f90 genwiq2.f90 writewiq2.f90 genjlgq0r.f90 \
 exxengy.f90 exxengyk.f90 oepmain.f90 oepresk.f90 oepvnl.f90 \
 oepvnlk.f90 oepmagmt.f90 oepmag.f90 zrhogp.f90

SRC_LDAU = \
 gendmatlu.f90 genveelu.f90 genvmatlu.f90 writeldapu.f90 readalphalu.f90 \
 tensmom.f90  dm2tensmom.f90  tensmomengy.f90 dmplz.f90 fluykwa.f90 \
 fluykwa0.f90 genflufr.f90 engyflu.f90 writeengyflu.f90 zbessela.f90 \
 zbesselb.f90 readdmatlu.f90 genflu.f90 findlambdalu.f90 readflu.f90 \
 writetensmom.f90

SRC_rdmft = \
 rdmft.f90 rdmdkdc.f90 rdmminc.f90 rdmvaryc.f90 rdmdedc.f90 \
 rdmenergy.f90 rdmengyxc.f90 rdmwritededn.f90 rdmwriteengy.f90 rdmminn.f90 \
 rdmvaryn.f90 rdmdedn.f90 rdmdexcdn.f90 rdmgenvnl_ijji.f90 rdmgenvnl_ijjk.f90 \
 rdmgetvnl_ijji.f90 rdmgetvnl_ijjk.f90 rdmputvnl_ijji.f90 rdmputvnl_ijjk.f90 \
 rdmdexcdc.f90 rdmdtsdn.f90 rdmentropy.f90 rdmeval.f90

SRC_pw = \
 writewfpw.f90 genwfpw.f90 getwfpw.f90 yambo.f90

SRC_lib = \
 euler.f90 wigner3j.f90 gaunt.f90 gauntyry.f90 r3mm.f90 r3mtm.f90 r3mmt.f90 \
 r3mv.f90 r3mtv.f90 r3cross.f90 r3dist.f90 r3taxi.f90 r3minv.f90 \
 r3mdet.f90 r3frac.f90 i3mdet.f90 i3mtv.f90 factnm.f90 factr.f90 hermite.f90 \
 brzint.f90 sphcrd.f90 sphcover.f90 erf.f90 clebgor.f90 wigner6j.f90 \
 sbessel.f90 sbesseldm.f90 genylm.f90 genrlm.f90 zmatinp.f90 \
 lopzflm.f90 sortidx.f90 gcd.f90 zfmtinp.f90 rfmtinp.f90 findband.f90 \
 gradzfmt.f90 gradrfmt.f90 ztorflm.f90 rtozflm.f90 zflmconj.f90 rotzflm.f90 \
 polynom.f90 sdelta.f90 stheta.f90 sdelta_mp.f90 stheta_mp.f90 sdelta_fd.f90 \
 stheta_fd.f90 sdelta_sq.f90 stheta_sq.f90 rdiracint.f90 rdiracdme.f90 \
 rdirac.f90 rschrodint.f90 rschroddme.f90 rschrodapp.f90 reciplat.f90 \
 connect.f90 flushifc.f90 spline.f90 rfinterp.f90 rfmtctof.f90 \
 fderiv.f90 fsmooth.f90 rotaxang.f90 i3minv.f90 axangsu2.f90 z2mm.f90 \
 z2mctm.f90 z2mmct.f90 potnucl.f90 timesec.f90 rfirvec.f90 vecfbz.f90 \
 mixadapt.f90 mixpulay.f90 ylmrot.f90 ylmroty.f90 rlmrot.f90 rotrflm.f90 \
 sdelta_lr.f90 stheta_lr.f90 grad2rfmt.f90

SRC_xc = \
  xc_pzca.f90 xc_pwca.f90 xc_pbe.f90 xc_am05.f90 xc_xalpha.f90 xc_wc06.f90 \
  x_wc06.f90 x_pbe.f90 c_pbe.f90 c_pbe_gcor.f90 xc_vbh.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines) $(SRC_lib) $(SRC_advanced) \
 $(SRC_phonon) $(SRC_optics) $(SRC_xc) $(SRC_LDAU) $(SRC_rdmft) $(SRC_pw)

OBJ = $(SRC:.f90=.o)
EXE = elk

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
elk:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_LPK) $(LIB_FFT) $(LIB_XC) 
         $(LIB_SYS)

blas:
	cd BLAS; $(MAKE); cp blas.a ..

lapack:
	cd LAPACK; $(MAKE); cp lapack.a ..

fft:
	cd fftlib; $(MAKE); cp fftlib.a ..

all:	blas lapack fft elk

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out *.aux *.dvi *.log \
         elk.pdf *.tex *.toc $(EXE)

cleanall:
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	$(MAKE) clean

doc:
	rm -f elk.aux elk.bbl elk.blg elk.log elk.pdf elk.tex elk.toc elk.lst
	ls $(SRC_modules) $(SRC_routines) $(SRC_lib) $(SRC_xc) $(SRC_advanced) $(SRC_phonon) $(SRC_optics) $(SRC_LDAU) $(SRC_rdmft) > elk.lst
	./protex -F -s $(SRC_main) $$(cat elk.lst) > elk.tex
	pdflatex elk;pdflatex elk;pdflatex elk

backup:
	tar -czf elk.tgz $(SRC) fftlib eos spacegroup junk BLAS LAPACK \
         Makefile notes.txt docs ../examples ../tests protex COPYING \
         README elk_silhouette.pdf depend ../release gentest genalltests \
         checklist.txt rmspaces

VERSION = $$(awk -F"/" '/data version/ {print $$2}' modmain.f90 | sed 's/ //g;s/,/./g')$
RELEASE = ../release/elk

release:
	rm -rf $(RELEASE)
	mkdir $(RELEASE)
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	cd eos; $(MAKE) clean
	cd spacegroup; $(MAKE) clean
	mkdir $(RELEASE)/src
	cp -r BLAS $(RELEASE)/src
	cp -r LAPACK $(RELEASE)/src
	cp -r fftlib $(RELEASE)/src
	cp -r eos $(RELEASE)/src
	cp -r spacegroup $(RELEASE)/src
	cp $(SRC) Makefile protex elk_silhouette.pdf $(RELEASE)/src
	cp libxc.inc libxc_funcs.f90 libxc.f90 libxcifc.f90 $(RELEASE)/src
	cd ../examples; find . -type f -name *~ -exec rm -f {} \;
	cp -r ../examples $(RELEASE)
	cd ../tests; find . -type f -name *~ -exec rm -f {} \;
	cp -r ../tests $(RELEASE)
	mkdir $(RELEASE)/species
	cp ../species/*.in $(RELEASE)/species
	cp -r ../utilities $(RELEASE)
	cp ../COPYING $(RELEASE)
	cp ../README $(RELEASE)
	cp release_notes.txt $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../setup $(RELEASE)
	$(MAKE) doc
	cp elk.pdf ../docs
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../../docs;$(MAKE) clean
	cp -r ../docs $(RELEASE)
	tar -C ../release -czf ../release/elk-$(VERSION).tgz elk
	cp ../docs/elk.pdf ../release
	cp ../docs/spacegroup.pdf ../release

lines:
	cat $(SRC) | wc -l

spaces:
	./rmspaces $(SRC)

