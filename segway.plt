reset
set xlabel "Time [s]" 
set grid 

set term postscript eps color "Times-Roman" 20  
set output "figure1.eps" 
set ylabel "displacements" 
plot 'segway.res' using 1:2 title 'q_0' with line , 'segway.res' using 1:5 title 'q_1' with line , 'segway.res' using 1:8 title 'q_2' with line , 'segway.res' using 1:11 title 'q_3' with line , 'segway.res' using 1:14 title 'q_4' with line , 'segway.res' using 1:17 title 'q_5' with line , 'segway.res' using 1:20 title 'q_6' with line , 'segway.res' using 1:23 title 'q_7' with line 
set term pop 
replot 
pause -1 'Next plot (velocity level)?' 

set term postscript eps color "Times-Roman" 20  
set output "figure2.eps" 
set ylabel "velocities" 
plot 'segway.res' using 1:3 title 'qd_0' with line , 'segway.res' using 1:6 title 'qd_1' with line , 'segway.res' using 1:9 title 'qd_2' with line , 'segway.res' using 1:12 title 'qd_3' with line , 'segway.res' using 1:15 title 'qd_4' with line , 'segway.res' using 1:18 title 'qd_5' with line , 'segway.res' using 1:21 title 'qd_6' with line , 'segway.res' using 1:24 title 'qd_7' with line 
set term pop 
replot 
pause -1 'Next plot (acceleration level)?' 

set term postscript eps color "Times-Roman" 20  
set output "figure3.eps" 
set ylabel "accelerations" 
plot 'segway.res' using 1:4 title 'qdd_0' with line , 'segway.res' using 1:7 title 'qdd_1' with line , 'segway.res' using 1:10 title 'qdd_2' with line , 'segway.res' using 1:13 title 'qdd_3' with line , 'segway.res' using 1:16 title 'qdd_4' with line , 'segway.res' using 1:19 title 'qdd_5' with line , 'segway.res' using 1:22 title 'qdd_6' with line , 'segway.res' using 1:25 title 'qdd_7' with line 
set term pop 
replot 
pause -1 
