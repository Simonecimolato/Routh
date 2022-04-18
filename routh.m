% ROUTH TABLE AND STABILITY CALCULATOR
% Simone Cimolato, 10/03/2022

close all
clear
clc

%%  INPUT DATA
coeffs = [1 3 4 2]; % in order, include 0s, [1 2 2 4 1 2] says it's stable but it's not, it has repeated poles on the imaginary axis

%%  ROUTH TABLE CONSTRUCTION
R = zeros(length(coeffs)-1+1,ceil((length(coeffs)-1+1)/2)+1);

R(:,1) = length(coeffs)-1:-1:0;
R(1,2:end) = coeffs(1:2:end);

%checking if length(coeffs) is even or odd
if rem(length(coeffs),2)
    R(2,2:end-1) = coeffs(2:2:end);
    R(2,end) = 0;
else
    R(2,2:end) = coeffs(2:2:end);
end

%adding one column of zeros to avoid errors in the routh formula
R(:,size(R,2)+1) = zeros(size(R,1),1);

for i = 3:size(R,1)
    for j = 2:size(R,2)-1
        R(i,j) = (-1/R(i-1,2)) * det([R(i-2,2) R(i-2,j+1); R(i-1,2) R(i-1,j+1)]); %routh table formula
        if any(R(i,2:end))         %correcting first row zeros, if the row is all zeros skip correction
            pos = find(R(i,2:end));
            R(i,2) = R(i,pos(1)+1)*((-1)^pos(1));
        end
    end
end

%removing first and last column 
R = R(1:end,2:end-1);

%checking if a row is all zeros or first element is NaN and only saving previous rows in Routh table
for i = 1:size(R,1)
    if isnan(R(i,1)) || not(any(R(i,:)))
        R = R(1:(i-1),:);
        break
    end
end

%%  OUTPUT
disp('Routh table: ')
disp(R);

%checking if the system is stable by checking the sign of the first column
if any(diff(sign(R(:,1).*(R(:,1)~=0))))
    disp('Unstable system')
else
    disp('Stable system')
end
