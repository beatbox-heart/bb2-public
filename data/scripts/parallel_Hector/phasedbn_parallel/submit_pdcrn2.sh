#!/bin/bash
# runs.sh
# Hector/Clusters multiple submits for qsub.

j=8
k=8
    dir=run${i}
    mkdir -p ${dir}
    mkdir -p ${dir}/ppm
    mkdir -p ${dir}/out
#   Copy the bbs files, the par files, the geometry, the Beatbox binary to the directory.
    cp $HOME/bin/Beatbox ${dir} 
    cp pd_crn2.bbs ${dir}
    cp pd_crn1.rec ${dir}

cat >Beatbox${dir}.sh <<EOF
#!/bin/bash --login
#PBS -N Beatbox${i}
#PBS -l mppwidth=${j}
#PBS -l mppnppn=${k}

#PBS -l walltime=00:20:00
#PBS -A e203

cd $HOME/work/crn2d/${dir}

ulimit -s unlimited

# n is the total number of processes
# N is the number of processes per node
# Launch the parallel job
aprun -n ${j} -N ${k} ./Beatbox crn2.bbs ${gk1}
EOF
qsub Beatbox${dir}.sh
done
