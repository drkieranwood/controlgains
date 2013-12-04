close all;

%=============================================
% Create system model
%=============================================

%Setup for control system
delayHd  = 0.5;
sampleTs = 0.1;
ptam_on = 0;       %changing this to 1 adds extra delay and slightly more measurement noise.

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
%check the discretisation process. Three continuous-time systems are
%created. One for position, one for velocity and one for acceleration.
%These are used to check the states of the discrete system actually
%represent these values, and not some arbitrary state selected by Matlab's
%c2d process.
tfcd = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1 0 0],'inputdelay',delay);
Ac=[-2.0*zeta*wn -1.0*wn*wn 0 0;1 0 0 0;0 1 0 0;0 0 1 0];
Bc=[kc*wn*wn;0;0;0];
Cc=[0 0 0 1];
Dc=[0];
ssc = ss(Ac,Bc,Cc,Dc);

tfc  = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1]);
tfcvel  = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1 0]);
sscvel = ss(tfcvel);
tfcaccel  = kc*tf(1,[1/(wn*wn) (2.0*zeta)/wn 1])*tf(1,[1]);
sscaccel = ss(tfcaccel);

h1 = figure('name','Cont. vs Disc Model');
hold on;
step(0.1*tfcd,3);


%====================
% Ad, Bd, Cd
%====================
%Now discretise the system with delay. Need the full time step Ad and Bd
%and Cd matrices first.
[Ad Bd Cd Dd] = discmat(ssc.a,ssc.b,ssc.c,sampleTs,0);

%====================
% Bd2
%====================
%Need Ad and Bd for Ts-gd
[Adt Bdt Cdt Ddt] = discmat(ssc.a,ssc.b,ssc.c,sampleTs-delayGd,0);
Adtemp = Adt;
Bd2 = Bdt;

%====================
% Bd1
%====================
%Need Bd for sampling interval gd
[Adt Bdt Cdt Ddt] = discmat(ssc.a,ssc.b,ssc.c,delayGd,0);
Bd1 = Adtemp*(Bdt);


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

SSD = ss(AD,BD,CD,DD,sampleTs);
figure(h1);
step(0.1*SSD,3);

sim('testDisc',3);
h1 = figure('name','Compare States');
ax1 = subplot(3,1,1);hold on;
plot(posCompare.time,posCompare.signals.values(:,1),'-r');
stairs(posCompare.time,posCompare.signals.values(:,2),'-g');
ax2 = subplot(3,1,2);hold on;
plot(velCompare.time,velCompare.signals.values(:,1),'-r');
stairs(velCompare.time,velCompare.signals.values(:,2),'-g');
ax3 = subplot(3,1,3);hold on;
plot(accelCompare.time,accelCompare.signals.values(:,1),'-r');
stairs(accelCompare.time,accelCompare.signals.values(:,2),'-g');
linkaxes([ax1 ax2 ax3],'x');

%At this point all the matrices for an LQG design with input noise have
%been created. Hence create the control.



%=============================================
% Cost and noise (tuning matrices)
%=============================================
%Noise standard distributions
inputNoise = 0.1;
if ptam_on
    measuNoise = sqrt(0.0039);
else
    measuNoise = 0.001;   %i.e. using Vicon so very accurate
end
%Convert to covariance matrices
Bdnoise = [Bd ; zeros(n+1,1)];
QE = Bdnoise*inputNoise*inputNoise*Bdnoise';
RE = measuNoise*measuNoise;

%State and input weightings (jerk, accel, velocity, position) 1/(maxdev^2)
QR = zeros(kl(AD),kl(AD));
QR(1,1) = 1/(1.0^2);
QR(2,2) = 1/(1.0^2);
QR(3,3) = 1/(1.0^2);
QR(4,4) = 1/(0.02^2);
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
if 1
    closedLoop1 = feedback(series(reg1,SSD),1,+1);
    closedLoop2 = feedback(series(reg2,SSD),1,+1);
    h2 = figure('name','Closed Loop');
    hold on;grid on;
    TT = 0:sampleTs:10;
    UU = 0.1*ones(1,length(TT));
    YY1=lsim(closedLoop1,UU,TT);
    YY2=lsim(closedLoop2,UU,TT);
    plot(TT,YY1,'-r');
    plot(TT,YY2,'-g');
    legend('LQG','H2');
end


%Create the LQR and LQE gains
Kc = dlqr(SSD.a,SSD.b,QR,RR);
Lc1 = dlqe(Ad,eye(length(Ad)),Cd,Bd*inputNoise*inputNoise*Bd',RE);
Lc2 = dlqe(Ad,Bd,Cd,inputNoise*inputNoise,RE);


%Compare the two state-space systems. Should be a small number.
ctrlcomp = sum(sum(reg1.a-reg2.a)) + sum(sum(reg1.b-reg2.b)) + sum(sum(reg1.c-reg2.c)) + sum(sum(reg1.d-reg2.d))


%==================================
%==================================
%Simulate in Simulink
matErrA=1.0;
matErrB=1.0;
matErrC=1.0;
stepDist = 0.01;
sim('testXY',12);
h363 = figure('name','Step Closed Loop Sim.');
ax1 = subplot(3,1,1);hold on;
plot(closedLoopSim.time,closedLoopSim.signals.values(:,1),'-r');
ax2 = subplot(3,1,2);hold on;
plot(closedLoopSimInp.time,closedLoopSimInp.signals.values(:,1),'-r');


%==================================
%==================================
%Simulate in code here.
%Need a system model to represent the ground truth, and the input reference
%positions.
clear XX YY UU RR TT corEst predEst delayStore;
ssTruth = SSD;
TT=0:sampleTs:12;
RR=zeros(length(AD),length(TT));
RR(4,:) = stepDist;
XX(:,1:1) = zeros(kl(AD),1);
YY(:,1:1) = zeros(1,1);
UU(:,1:1) = zeros(1,1);

%Need estimator states.
corEst(:,1)  = zeros(kl(Ad),1);
predEst(:,1) = zeros(kl(Ad),1); 
delayStore(:,1) = zeros(n+2,1)'

%Run a simulation for all time
for ii=1:1:(length(TT))
    %At ii=1 -> TT=0
   
    %Correct the state prediction to form the current state estimate used
    %for control.
    corEst(:,ii) = predEst(:,ii) + Lc2*( YY(:,ii) - Cd*predEst(:,ii) );
    
    %Create new control. First augment states, then find error.
    %The augmented state consists of the estimate and the previous n+1
    %inputs. The newly created control does not feature in order to avoid 
    %a loop.
    xAug = [corEst(:,ii) ; delayStore(2:n+2)];
    UU(:,ii) = Kc*(RR(:,ii) - xAug);
    if(UU(:,ii) > 0.37)
        UU(:,ii) = 0.37;
    end
    if (UU(:,ii)< -0.37)
        UU(:,ii) = -0.37;
    end
    
    %Now the control is found use it to move the real system state along.
    %Also create the measurement output for the future.
    XX(:,ii+1) = AD*XX(:,ii) + BD*UU(:,ii);
    YY(:,ii+1) = CD*XX(:,ii+1);
    
    %Update the delay storage. The new control is used and added to the
    %bottom.    
    delayStore(1:n+1) = delayStore(2:n+2);
    delayStore(n+2) = UU(:,ii);
    
    %Predict the future state using the stored inputs
    predEst(:,ii+1) = Ad*corEst(:,ii) + Bd1*delayStore(1) + Bd2*delayStore(2);
    

    
end

figure(h363);
subplot(3,1,1);hold on;
plot(TT(1:ii),XX(4,1:end-1),'-g');
subplot(3,1,2);hold on;
plot(TT(1:ii),UU,'-g');
subplot(3,1,3);hold on;
plot(TT(1:ii),corEst(4,:),'-g');



%Write gains to file
writeRedord(Ad,Bd1,Bd2,Cd,n,sampleTs,Kc,Lc1,'redordX');
writeSS(reg1,'controlX');

Lc1
Kc


