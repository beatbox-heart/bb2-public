#!/bin/bash
# Scaling simulation for 3D CRN complete box.
# this script, and bigbox_crn3d.bbs needs to be put in your work on Hector.

start=1
finish=3

for (( i = $start; i < $finish; i++  ))
do

if [ "$i" -eq 1 ]; then
    j=1
    k=1
elif [ "$i" -eq 2 ]; then
    j=128
    k=28
fi

    dir=$HOME/work/bigboxcrn_run${i}
    mkdir -p ${dir}
    mkdir -p ${dir}/ppm
    mkdir -p ${dir}/out
#   Copy the bbs files, the par files, the geometry, the Beatbox binary to the directory.
    cp $HOME/bin/Beatbox ${dir} 
    cp bigbox_crn.bbs ${dir}

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
aprun -n ${j} -N ${k} ./Beatbox bigbox_crn.bbs -verbose -profile
EOF
qsub Beatbox${i}.sh
done
