function [ output_args ] = writeSS( SSSYS , filename )
%WRITESS Write a state-space system to file
%   This function writes the size properties and parameters of a
%   discrete-time state-space system to file. The system to be written is
%   the first argument and the file name (as a string) is the second. This
%   script automatically overwrites previous files.
%
%   TODO: Make the file include text descriptions of the entries, rather
%   than only numeric values.

%First check if all the arguments are given
if ((nargin<2) || (nargin>2))
    disp('writeSS: Incorrect number of arguments');
    return;
end

%Check if the specified system is discrete-time and state-space
SSSYS.a;   %This breaks if there is no A matrix. A bit clumsy but it works.
if (SSSYS.Ts==0)
    disp('writeSS: Sample-time==0. Is this discrete-time?');
    return;
end

%Now find the state-space size
[a b] = size(SSSYS.a);
numSta = a;
[a b] = size(SSSYS.b);
numInp = b;
[a b] = size(SSSYS.c);
numOut = a;

%Open a file for the output. If exists then overwrite.
fd = fopen(filename,'w');

%Check if file opened
if (fd==-1)
    disp('writeSS: Could not open file for writing.');
    return;
end

%First write the matrix sizes onto the first three lines.
fprintf(fd,'%d\n',numSta);
fprintf(fd,'%d\n',numInp);
fprintf(fd,'%d\n',numOut);

%Write the sample rate.
fprintf(fd,'%20.10f\n',SSSYS.Ts);

%Work through the A matrix adding the values to file. This is in row order
%(i.e. first row, second row, third row, ...)
for ii=1:1:numSta
    for jj=1:1:numSta
        fprintf(fd,'%10.5f\n',SSSYS.a(ii,jj));
    end
end

%Write B matrix similarly
for ii=1:1:numSta
    for jj=1:1:numInp
        fprintf(fd,'%10.5f\n',SSSYS.b(ii,jj));
    end
end

%Write C matrix similarly
for ii=1:1:numOut
    for jj=1:1:numSta
        fprintf(fd,'%10.5f\n',SSSYS.c(ii,jj));
    end
end

%Write D matrix similarly
for ii=1:1:numOut
    for jj=1:1:numInp
        fprintf(fd,'%10.5f\n',SSSYS.d(ii,jj));
    end
end

%Close the file
fclose(fd);

%Return success
output_args = 1;
end

