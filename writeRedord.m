function [ output_args ] = writeRedord( Ad,Bd1,Bd2,Cd,nd,Ts,Kc,Lc,filename )
%WRITESS Write a state-space system to file of teh reduced order SISO
%estimator
%   This function writes the size properties and parameters of a
%   discrete-time state-space system to file. The system to be written is
%   the first argument and the file name (as a string) is the second. This
%   script automatically overwrites previous files.
%
%   TODO: Make the file include text descriptions of the entries, rather
%   than only numeric values.

%First check if all the arguments are given
if ((nargin<9) || (nargin>9))
    disp('writeRedord: Incorrect number of arguments');
    return;
end

%Check if the specified system is discrete-time and state-space
if (Ts==0)
    disp('writeRedord: Sample-time==0. Is this discrete-time?');
    return;
end

%Check nd is zero or greater
if (nd<0)
    disp('writeRedord: Number of delays must be positive');
    return;
end

%Now find the state-space size
[a b] = size(Ad);
numSta = a;
[a b] = size(Bd1);
numInp = b;
[a b] = size(Cd);
numOut = a;

%Open a file for the output. If exists then overwrite.
fd = fopen(filename,'w');

%Check if file opened
if (fd==-1)
    disp('writeRedord: Could not open file for writing.');
    return;
end

%First write the matrix sizes onto the first three lines.
fprintf(fd,'%d\n',numSta);
fprintf(fd,'%d\n',numInp);
fprintf(fd,'%d\n',numOut);
fprintf(fd,'%d\n',nd);

%Write the sample rate.
fprintf(fd,'%20.10f\n',Ts);

%Work through the Ad matrix adding the values to file. This is in row order
%(i.e. first row, second row, third row, ...)
for ii=1:1:numSta
    for jj=1:1:numSta
        fprintf(fd,'%10.5f\n',Ad(ii,jj));
    end
end

%Write Bd1 matrix similarly
for ii=1:1:numSta
    for jj=1:1:numInp
        fprintf(fd,'%10.5f\n',Bd1(ii,jj));
    end
end
%Write Bd2 matrix similarly
for ii=1:1:numSta
    for jj=1:1:numInp
        fprintf(fd,'%10.5f\n',Bd2(ii,jj));
    end
end

%Write Cd matrix similarly
for ii=1:1:numOut
    for jj=1:1:numSta
        fprintf(fd,'%10.5f\n',Cd(ii,jj));
    end
end

%Write Kc matrix similarly
for ii=1:1:numInp
    for jj=1:1:(numSta+nd+1)
        fprintf(fd,'%10.5f\n',Kc(ii,jj));
    end
end

%Write Lc matrix similarly
for ii=1:1:numSta
    for jj=1:1:numOut
        fprintf(fd,'%10.5f\n',Lc(ii,jj));
    end
end

%Close the file
fclose(fd);

%Return success
output_args = 1;
end

