close all
clear all

cum_aerr=zeros(3,70);
for run=1:500
  aerr=ones(3,70);
  for index0=1:70
    SNR=10^((index0-20)/20);
    aerr(:,index0)=func(SNR);
  end
  cum_aerr=cum_aerr+aerr;
end
dB=([1:70]-20)./20;
figure
plot(dB,cum_aerr/500)
xlabel('SNR (dB)')
ylabel('Variance')
title('Variance vs SNR')

function aerr=func(SNR)
%Angle of Arrivals
th1=10;
th2=15;
th3=20;
aoas=[th1 th2 th3]*pi/180;
q=length(aoas); %Number of signals

%Power
P=SNR*ones(1,q);

%Other Params
N=2000; %Number of samples
d=0.5;
p=5; %Number of antennas
K=10;
noise_var=1;
fc=28*10^5; %Frequency in Hz
t=0:1e-9:N*1e-9-1e-9; %Time
 
%Steering Vector - A is p*q, p antennas, q signals 
A=exp(-1i*2*pi*d*(0:p-1)'*sin(aoas));

%Derivative of Steering Vector -
dA=-1i*2*pi*d*(0:p-1)'*cos(aoas).*A;

%Covariance Matrix
C=zeros(p);
for index=1:q
    C=C+P(index)*A(:,index)*A(:,index)';
end
C=C+eye(p);
C_inv=inv(C);


%Fisher
FiM=zeros(q);
for index=1:q
    for index1=1:q
        FiM(index,index1)=trace(C_inv*deriv_C(index,A,dA,P)*C_inv*deriv_C(index1,A,dA,P));
    end
end

%CRLB
FiM_inv=inv(FiM);
aerr=diag(FiM_inv);
aerr=abs(aerr);
end 

function dC=deriv_C(index,A,dA,P)
    dC=P(index)*dA(:,index)*A(:,index)'+P(index)*A(:,index)*dA(:,index)';
end

