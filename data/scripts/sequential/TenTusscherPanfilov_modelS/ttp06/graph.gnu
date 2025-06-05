set term post eps color 'Times-Roman' 18 size 2,1.5
unset key;

set xlabel 'distance (mm)'
set xrange [0:21.4]
set xtics 10.7;
set mxtics 2;
set grid x;

set ylabel 'activation time (ms)' offset 1,0
set yrange [0:150];
set ytics 50;
set grid y;


set style data lines;
set style line 1 lt 1 lc 1 lw 3;
set style line 2 lt 1 lc 2 lw 3;
set style line 3 lt 1 lc 3 lw 3;
plot '-' ls 1, '-' ls 2, '-' ls 3



