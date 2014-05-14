timingf = [0.011760,0.077914, 1.422031,11.798825,60.582009]
timingc = [0.009619,0.175789,2.113036,16.767863]
timinga = [0.002079,0.005414,0.016068,0.047425]
ftcs = [5,10,20,30,50]
cranktime = [3,4,5,6]
aditime = [3,4,5,6]

figure;
subplot(3,1,1);
plot(ftcs, timingf)

title('FTCS');


subplot(3, 1, 2);
plot(cranktime, timingc)

title('Crank-Nicholson');



subplot(3, 1, 3);
plot(aditime,timinga)
xlabel('Size of cube NxNxN');
ylabel('Time');
title('ADI Method');





