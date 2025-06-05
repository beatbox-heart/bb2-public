set term post eps color enh 'Helvetica' 18 size 3.5in,3in;
set out 'rk4test.eps';

set key top left reverse Left spacing 1.5;
set logscale xy;

set xlabel 'h, time step';
xmin=5e-5; xmax=3e-1;
set xrange [xmin:xmax];
set format x '10^{%T}';
xmid=sqrt(xmin*xmax);

set ylabel 'error norm';
ymin=1e-18; ymax=2e-4;
set yrange [ymin:ymax];
set format y '10^{%T}';
ymid=sqrt(ymin*ymax);

plot \
  'rk4test.txt' u 2:3 tit 'rk4' w lp lc rgbc 'red', \
  ymid*(x/xmid)**4 tit 'h^4' w l dt 2 lc rgbc 'blue', \
  1/0 tit ''

