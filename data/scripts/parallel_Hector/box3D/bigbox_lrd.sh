#!/bin/bash
# Scaling simulation for 3D LRD complete box.
# this script, and bigbox_crn3d.bbs needs to be put in your work on Hector.

start=1
finish=3

for (( i = $start; i < $finish; i++  ))
do

if [ "$i" -eq 1 ]; then
    j=1
    k=1
elif [ "$i" -eq 2 ]; then
    j=712
    k=28
fi

    dir=$HOME/work/bigboxcrn_run${i}
    mkdir -p ${dir}
    mkdir -p ${dir}/ppm
    mkdir -p ${dir}/out
#   Copy the bbs files, the par files, the geometry, the Beatbox binary to the directory.
    cp $HOME/bin/Beatbox ${dir}
    cp $HOME/beatbox/data/parameters/lrd.par ${dir}
    cp bigbox_lrd.bbs ${dir}

cat >Beatbox${i}.sh <<EOF
#!/bin/bash --login
#PBS -N Beatbox${i}
#PBS -l mppwidth=${j}
#PBS -l mppnppn=${k}

#PBS -l walltime=1:00:00
#PBS -A e203

cd ${dir}

ulimit -s unlimited

# n is the total number of processes
# N is the number of processes per node
# Launch the parallel job
aprun -n ${j} -N ${k} ./Beatbox bigbox_lrd.bbs -verbose -profile
EOF
qsub Beatbox${i}.sh
done
