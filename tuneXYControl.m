close all;

%=============================================
% Create system model
%=============================================

%Setup for control system
delayHd  = 0.12;
sampleTs = 0.2;
ptam_on = 0;

%Best parameters from the step responses.
%wn and zeta are the natural frequency and damping ratio of a second order
%system. kc is the gain and delay is the loop delay when flying with Vicon.
%ptamDelay is the additional feedback delay when using ptam.
wn   = 7.27;
zeta = 0.542;
kc   = 0.532*9.81;
delay = 0.08 + delayHd;
if ptam_on
    delay = delay + 0.14;
end


%Find gd and n
n=0;
while ((delay-(n*sampleTs))>sampleTs)
    n=n+1;
end
delayGd = delay-(n*sampleTs);


%First create a continuous version of the system. This is used later to
%check the discretisation process.
tfcd = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1 0 0],'inputdelay',delay);
tfc  = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1 0 0]);
ssc = ss(tfc);

tfcvel  = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1 0]);
sscvel = ss(tfcvel);

tfcaccel  = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1]);
sscaccel = ss(tfcaccel);

h1 = figure('name','Cont. vs Disc Model');
hold on;
step(0.1*tfcd,3);

%Now discretise the system with delay. Need the full time step Ad and Bd
%and Cd matrices first.
ssd = c2d(ssc,sampleTs);
Ad = ssd.a;
Bd = ssd.b;
Cd = ssd.c;

TT = eye(length(Ad));
TT(1,1) = Cd(4);
TT(2,2) = Cd(4);
TT(3,3) = Cd(4);
TT(4,4) = Cd(4);
Bd = TT*Bd;

%Need Ad and Bd for Ts-gd
ssd = c2d(ssc,(sampleTs-delayGd));
Adtemp = ssd.a;
Bd2 = ssd.b;

%Need Bd for gd
ssd = c2d(ssc,delayGd);
Bd1 = Adtemp*(ssd.b);

%Now have Ad, Bd1, Bd2, Bd, Cd. Construct into the augmented form.
if (n==0)
    disp('n==0');
    %If the delay is less than a single sample period.
    AD = [Ad Bd1; zeros(1,kw(Ad)) 0];
    BD = [Bd2 ; 1];
    CD = [Cd 0];
    DD = [0];    
else
    disp('n>0');
    %If the delay spans multiple sample periods.
    %Need n+1 extra states adding.
    AD = [Ad Bd1 Bd2 zeros(kl(Ad),(n-1)) ; zeros(n+1,kw(Ad)+2) zeros(n+1,(n-1))];
    BD = [zeros(kl(AD)-1,1);1];
    CD = [Cd zeros(1,n+1)];
    DD = [0];
    AD(kl(Ad)+1:kl(Ad)+n,kl(Ad)+2:kl(Ad)+1+n) = eye(n);  
end

TT = eye(kl(AD));
TT(1,1) = CD(4);
TT(2,2) = CD(4);
TT(3,3) = CD(4);
TT(4,4) = CD(4);
AD = TT*AD*inv(TT);
BD = TT*BD;
CD = CD*inv(TT);


SSD = ss(AD,BD,CD,DD,sampleTs);
figure(h1);
step(0.1*SSD,3);

%At this point all the matrices for an LQG design with input noise have
%been created. Hence create the control.

%=============================================
% Cost and noise (tuning matrices)
%=============================================
%Noise standard distributions
inputNoise = 0.05;
if ptam_on
    measuNoise = sqrt(0.0039);
else
    measuNoise = 0.001;   %i.e. using Vicon so very accurate
end
%Convert to covariance matrices
Bdnoise = [Bd ; zeros(n+1,1)];
QE = Bdnoise*inputNoise*inputNoise*Bdnoise';
RE = measuNoise*measuNoise;

%State and input weightings (2nd ord lag, accel, velocity, and position) 1/maxdev^2
QR = zeros(kl(AD),kl(AD));
QR(1,1) = 1/(0.1^2);
QR(2,2) = 1/(0.1^2);
QR(3,3) = 1/(0.1^2);
QR(4,4) = 1/(0.01^2);
RR = 1/(0.37^2);

QXU = blkdiag(QR,RR);
QWV = blkdiag(QE,RE);

reg1 = lqg(SSD,QXU,QWV);


%=============================================
% Create h2 system
%=============================================
h2AA = AD;
h2BB = cat(2,Bdnoise*inputNoise,zeros(kl(SSD.a),1));
h2BB = cat(2,h2BB,BD);
h2CC = cat(1,zeros(1,kl(SSD.a)),sqrt(QR));
h2CC = cat(1,h2CC,CD);
h2DD = zeros(kl(h2CC),kw(h2BB));
h2DD(1,3) = sqrt(RR);
h2DD(end,2) = measuNoise;

ssh2 = ss(h2AA,h2BB,h2CC,h2DD,sampleTs);
[reg2,CL,GAM,INFO]=h2syn(ssh2,1,1);


%=============================================
% Test closed loop step response
%=============================================
closedLoop1 = feedback(series(reg1,SSD),-1);
closedLoop2 = feedback(series(reg2,SSD),-1);
h2 = figure('name','Closed Loop');
hold on;
step(closedLoop1,10,'-xr');
step(closedLoop2,10,'-g');
legend('LQG','H2');


ctrlcomp = sum(sum(reg1.a-reg2.a)) + sum(sum(reg1.b-reg2.b)) + sum(sum(reg1.c-reg2.c)) + sum(sum(reg1.d-reg2.d))