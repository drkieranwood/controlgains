disp('Creating control...');



close all;
clear all;
globalRate = 0.1;
globalDelay = 0.0;
tuneXYControl;
writeSS(reg1,'../controlgains/controlX');
writeSS(reg1,'../controlgains/controlY');

close all;
clear all;
globalRate = 0.1;
globalDelay = 0.0;
tuneZControl;
writeSS(reg1,'../controlgains/controlZ');

close all;
clear all;
globalRate = 0.1;
globalDelay = 0.0;
tuneWControl;
writeSS(reg1,'../controlgains/controlW');

close all;
clear all;
disp('...done');