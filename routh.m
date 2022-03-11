% ROUTH TABLE AND STABILITY CALCULATOR
% Simone Cimolato, 10/03/2022

close all
clear
clc

%%  INPUT DATA
coeffs = [8 10 21 4 9]; % in order, include 0s

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
        %routh table formula
        R(i,j) = (-1/R(i-1,2)) * det([R(i-2,2) R(i-2,j+1); R(i-1,2) R(i-1,j+1)]);
        %correcting first row zeros
        if R(i,2) == 0
            k = 2;
            %if the row is all zeros skip correction
            while R(i,k) == 0 && k < length(R(i,:))
                k = k+1;
            end
            if k ~= length(R(i,:)) 
                R(i,2) = R(i,k)*((-1)^(k-1));
            end
        end
    end
end

%removing first and last column 
R = R(1:end,2:end-1);

%checking if a row is NaN or all zeros and only saving previous rows in
%Routh table
for i = 1:size(R,1)
    if isnan(R(i,1)) || not(any(R(i,:)))
        R = R(1:(i-1),:);
        break
    end
end

%%  OUTPUT
disp('Routh table: ')
disp(R);

%checking if the system is stable or not by checking the sign of the first
%column
if any(diff(sign(R(:,1).*(R(:,1)~=0))))
    disp('Unstable system')
else
    disp('Stable system')
end
