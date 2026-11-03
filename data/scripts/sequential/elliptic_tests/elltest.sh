#!/bin/bash
# Run one case of the elliptic device tests and print its row of results;
# used by the Makefile in this directory (see there for the cases).
#
# Usage: elltest.sh header plane|neum|coef|isol|mpi   print the table header
#        elltest.sh plane NAME nx ny nz th0 th1 al rhs ul gs
#        elltest.sh neum  NAME nx ny nz shape th0 al ul refcyc
#        elltest.sh coef  NAME kind script coef1 coef2 nx ny nz rhs ul
#        elltest.sh isol  NAME nx ny nz shape th0 al ul xb
#        elltest.sh mpi   KIND NAME nx ny nz th0 th1 al shape ul   (one row per np in $NPS)
# Environment: BBXS, BBXP (sequential and MPI Beatbox), MPIRUN (e.g. "mpirun -np"),
#   NPS (numbers of processes), TLIM (time limit of a run, seconds), WORK (work directory).
#
# Columns:
#   Eu    max|u - u_exact| (the device's solution vs the exact solution of the PDE)
#   Eg    the same for the plain GS solution (upper_level=-1), if run
#   Epg   max|u_MG - u_GS|, should be ~tolerance
#   Er    max|u_MG - u_ref|, reference: upper_level=0, many sweeps, if run
#   rres  final residual of the reference solve (Er is meaningful only if small)
#   Pp    u_MG(pin) - u_exact(pin), should be 0
#   dEu   |Eu - Eu(np=1)| (mpi), or |Eu(sub) - Eu(cut)| (isol), should be 0
#   ncyc  number of V-cycles to reach the tolerance
#   rho   mean V-cycle residual reduction factor, from the 2nd cycle on
set -u
set -f				# coefficients contain * and ()
here=$(cd "$(dirname "$0")" && pwd)
BBXS=${BBXS:-bb2s}
BBXP=${BBXP:-bb2p}
MPIRUN=${MPIRUN:-mpirun -np}
NPS=${NPS:-1 2 4}
TLIM=${TLIM:-300}
WORK=${WORK:-work}
mkdir -p "$WORK"; cd "$WORK" || exit 1

# Run a command with a time limit of TLIM seconds (macOS has no timeout(1));
# the watchdog's output goes to /dev/null so that $(...) does not wait for it
runlim() {
  "$@" & local pid=$!
  ( sleep $TLIM; kill -TERM $pid; sleep 5; kill -KILL $pid ) >/dev/null 2>&1 &
  local w=$!
  wait $pid; local rc=$?
  kill $w 2>/dev/null; wait $w 2>/dev/null
  return $rc
}
# V-cycle count and mean reduction factor (from the 2nd cycle on) from a debug file
vstats() {
  grep -E 'V-cycle [0-9]+ done, rsn\[0\]' $1 2>/dev/null | \
    sed -E 's/.*rsn\[0\]=([^ ]*).*/\1/' | \
    awk '{r[NR]=$1} END{n=NR; if(n<3){print n, "-"; exit}
         printf "%d %.3f\n", n, exp((log(r[n])-log(r[2]))/(n-2))}'
}
getval() { echo "$2" | sed -E "s/.*[ ]$1=([^ ]*).*/\1/"; }
absdiff() { awk -v a=$1 -v b=$2 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.3g", d}'; }
# The last result line of a run, or FAILED with the reason
result() {
  local out="$1" res
  res=$(echo "$out" | grep -E 'Eu=' | tail -1)
  if [ -z "$res" ]; then
    echo "FAILED: $(echo "$out" | grep -m1 -i -E 'error|abort|fault' || echo "timeout or crash (TLIM=$TLIM)")"
  else
    echo "$res"
  fi
}

case $1 in
header)
  case $2 in
  plane) printf "%-12s %4s %4s %4s %5s %5s %3s %3s %3s  %-11s %-11s %-11s %4s %6s\n" \
	   case nx ny nz th0 th1 al rhs ul Eu Eg Epg ncyc rho ;;
  neum)  printf "%-12s %4s %4s %4s %-4s %4s %3s %3s  %-11s %-11s %-11s %-11s %4s %6s\n" \
	   case nx ny nz shp th0 al ul Eu Er rres Pp ncyc rho ;;
  coef)  printf "%-8s %4s %4s %4s %4s %3s %3s  %-11s %-11s %-11s %4s %6s\n" \
	   case kind nx ny nz rhs ul Eu Eg Epg ncyc rho ;;
  isol)  printf "%-8s %4s %4s %4s %-4s %4s %3s %3s %3s  %-11s %-11s %-11s %4s %4s\n" \
	   case nx ny nz shp th0 al ul xb Eu_sub Eu_cut dEu ncs ncc ;;
  mpi)   printf "%-14s %3s %4s %4s %4s %-4s %3s  %-11s %-11s %4s %6s\n" \
	   case np nx ny nz shp ul Eu dEu ncyc rho ;;
  esac ;;

plane)	# Dirichlet, plane wave, ell_plane.bbs
  name=$2 nx=$3 ny=$4 nz=$5 th0=$6 th1=$7 al=$8 rhs=$9 ul=${10} gs=${11}
  "$here/mkgeom.pl" $nx $ny $nz box $th0 $th1 $al > $name.bbg
  ln -sf "$here/ell_plane.bbs" $name.bbs
  res=$(result "$(runlim $BBXS -nograph -- $name.bbs $name.bbg $nx $ny $nz $th0 $th1 $al $rhs $ul $gs 0 200 2>&1)")
  case $res in FAILED*) printf "%-12s %s\n" $name "$res"; exit 1;; esac
  Eg=$(getval Eg "$res"); Epg=$(getval Epg "$res")
  [ $gs = 0 ] && Eg=- && Epg=-
  read ncyc rho < <(vstats ${name}_p.dbg)
  printf "%-12s %4s %4s %4s %5s %5s %3s %3s %3s  %-11s %-11s %-11s %4s %6s\n" \
    $name $nx $ny $nz $th0 $th1 $al $rhs $ul $(getval Eu "$res") $Eg $Epg $ncyc $rho ;;

neum)	# no-flux (singular), pinned, ell_neum.bbs
  name=$2 nx=$3 ny=$4 nz=$5 shape=$6 th0=$7 al=$8 ul=$9 ref=${10}
  "$here/mkgeom.pl" $nx $ny $nz $shape $th0 0 $al > $name.bbg
  ln -sf "$here/ell_neum.bbs" $name.bbs
  res=$(result "$(runlim $BBXS -nograph -- $name.bbs $name.bbg $nx $ny $nz $th0 $al $ul $ref 0 400 0 0 2>&1)")
  case $res in FAILED*) printf "%-12s %s\n" $name "$res"; exit 1;; esac
  Er=$(getval Er "$res"); rres=-
  if [ $ref = 0 ]; then Er=-; else
    rres=$(grep -E 'V-cycle [0-9]+ done, rsn\[0\]' ${name}_g.dbg | tail -1 | sed -E 's/.*rsn\[0\]=([^ ]*).*/\1/')
  fi
  read ncyc rho < <(vstats ${name}_p.dbg)
  printf "%-12s %4s %4s %4s %-4s %4s %3s %3s  %-11s %-11s %-11s %-11s %4s %6s\n" \
    $name $nx $ny $nz $shape $th0 $al $ul $(getval Eu "$res") $Er $rres $(getval Pp "$res") $ncyc $rho ;;

coef)	# Dirichlet, plane wave, various coefficients, ell_coef*.bbs
  name=$2 kind=$3 script=$4 c1=$5 c2=$6 nx=$7 ny=$8 nz=$9 rhs=${10} ul=${11}
  for f in ell_coef.par ell_coef.inc; do ln -sf "$here/$f" .; done
  ln -sf "$here/$script" $name.bbs
  res=$(result "$(runlim $BBXS -nograph -- $name.bbs $kind $nx $ny $nz $c1 $c2 $rhs $ul 1 200 2>&1)")
  case $res in FAILED*) printf "%-8s %s\n" $name "$res"; exit 1;; esac
  read ncyc rho < <(vstats ${name}_p.dbg)
  printf "%-8s %4s %4s %4s %4s %3s %3s  %-11s %-11s %-11s %4s %6s\n" \
    $name $kind $nx $ny $nz $rhs $ul $(getval Eu "$res") $(getval Eg "$res") $(getval Epg "$res") $ncyc $rho ;;

isol)	# isolated=1 on x<=xb of the tissue ("sub") vs the tissue cut at xb ("cut")
  name=$2 nx=$3 ny=$4 nz=$5 shape=$6 th0=$7 al=$8 ul=$9 xb=${10}
  "$here/mkgeom.pl" $nx $ny $nz $shape $th0 0 $al > ${name}s.bbg
  "$here/mkgeom.pl" $nx $ny $nz $shape $th0 0 $al $xb > ${name}c.bbg
  ln -sf "$here/ell_neum.bbs" ${name}s.bbs
  ln -sf "$here/ell_neum.bbs" ${name}c.bbs
  ress=$(result "$(runlim $BBXS -nograph -- ${name}s.bbs ${name}s.bbg $nx $ny $nz $th0 $al $ul 0 0 400 $xb 1 2>&1)")
  resc=$(result "$(runlim $BBXS -nograph -- ${name}c.bbs ${name}c.bbg $nx $ny $nz $th0 $al $ul 0 0 400 $xb 1 2>&1)")
  for r in "$ress" "$resc"; do case $r in FAILED*) printf "%-8s %s\n" $name "$r"; exit 1;; esac; done
  Eus=$(getval Eu "$ress"); Euc=$(getval Eu "$resc")
  printf "%-8s %4s %4s %4s %-4s %4s %3s %3s %3s  %-11s %-11s %-11s %4s %4s\n" \
    $name $nx $ny $nz $shape $th0 $al $ul $xb $Eus $Euc $(absdiff $Eus $Euc) \
    $(vstats ${name}s_p.dbg | cut -d' ' -f1) $(vstats ${name}c_p.dbg | cut -d' ' -f1) ;;

mpi)	# a case run with each number of processes in NPS
  # kind d: Dirichlet (ell_plane.bbs); n: no-flux (ell_neum.bbs); c3: isotropic D=~,
  # c4: anisotropic Dl=~ (ell_coef*.bbs, with tensor layers d0..d1); i: no-flux on
  # x<=nx/2+4 with isolated=1 (ell_neum.bbs)
  kind=$2 name=$3 nx=$4 ny=$5 nz=$6 th0=$7 th1=$8 al=$9 shape=${10} ul=${11}
  "$here/mkgeom.pl" $nx $ny $nz $shape $th0 $th1 $al > $name.bbg
  case $kind in
  d)  ln -sf "$here/ell_plane.bbs" $name.bbs
      args="$name.bbg $nx $ny $nz $th0 $th1 $al 0 $ul 0 0 300" ;;
  n)  ln -sf "$here/ell_neum.bbs" $name.bbs
      args="$name.bbg $nx $ny $nz $th0 $al $ul 0 0 300 0 0" ;;
  c3) for f in ell_coef.par ell_coef.inc; do ln -sf "$here/$f" .; done
      ln -sf "$here/ell_coefiso.bbs" $name.bbs
      args="3 $nx $ny $nz D=~D0*(1+va*sin(bx*(x-1)*hx+by*(y-1)*hx)) manypoint=0 0 $ul 0 300" ;;
  c4) for f in ell_coef.par ell_coef.inc; do ln -sf "$here/$f" .; done
      ln -sf "$here/ell_coef.bbs" $name.bbs
      args="4 $nx $ny $nz Dl=~Dl0*(1+va*sin(bx*(x-1)*hx)) Dt=Dt0 0 $ul 0 300" ;;
  i)  ln -sf "$here/ell_neum.bbs" $name.bbs
      args="$name.bbg $nx $ny $nz $th0 $al $ul 0 0 300 $((nx/2+4)) 1" ;;
  *)  echo "unknown mpi kind $kind" >&2; exit 2 ;;
  esac
  Eu1= ; rc=0
  for np in $NPS; do
    res=$(result "$(runlim $MPIRUN $np $BBXP -nograph -- $name.bbs $args 2>&1)")
    case $res in FAILED*) printf "%-14s %3s %s\n" $name $np "$res"; rc=1; continue;; esac
    Eu=$(getval Eu "$res"); [ -z "$Eu1" ] && Eu1=$Eu
    read ncyc rho < <(vstats ${name}_p.dbg)
    printf "%-14s %3s %4s %4s %4s %-4s %3s  %-11s %-11s %4s %6s\n" \
      $name $np $nx $ny $nz $shape $ul $Eu $(absdiff $Eu $Eu1) $ncyc $rho
  done
  exit $rc ;;

*) echo "Usage: see the comment in $0" >&2; exit 2 ;;
esac
