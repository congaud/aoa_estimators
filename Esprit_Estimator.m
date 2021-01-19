close all
clear all

%Plots variance
x_hat_var=zeros(3,10);
for index0=4:14
    x_hat=zeros(3,500);
    for run=1:500
        SNR=10^((index0-20)/20);
        x_hat(:,run)=func(index0);
    end
    x_hat_mean = sum(x_hat,2)./500;
    temp = (x_hat-x_hat_mean).^2;
    x_hat_var(:,index0-3) = sum(temp,2)./500;
end

%Plots MSE
% cum_aerr=zeros(3,70);
% for run=1:500
%   aerr=ones(3,70);
%   for index0=1:70
%     SNR=10^((index0-20)/20);
%     aerr(:,index0)=func(SNR);
%   end
%   cum_aerr=cum_aerr+aerr;
% end
dB=([1:70]-20)./20;
figure
plot([4:14],x_hat_var)
xlabel('Number of Antennas')
ylabel('Variance')
title('Variance vs Number of Antennas')

function x_est=func(SNR)

th1=10.001;
th2=15.001;
th3=20.001;
p1=1;
p2=1;
p3=1;
K=2000;
d=0.5;
N=SNR;
fc=28e6;
t=0:1e-9:K*1e-9-1e-9; %Time
noise_var=1;
aoas=[th1 th2 th3]*pi/180; 
P=1*[p1 p2 p3];
r=length(aoas);

% Steering vector matrix. Columns will contain the steering vectors
% of the r signals
A=exp(-i*2*pi*d*(0:N-1)'*sin([aoas(:).']));

% Signal and noise generation
%sig=round(rand(r,K))*2-1; % Generate random BPSK symbols for each of the
% r signals
sig=[sin(2*pi*fc*t);sin(2*pi*1.1*fc*t+2*pi/3);sin(2*pi*1.2*fc*t+3*pi/5)];

noise=sqrt(noise_var/2)*(randn(N,K)+i*randn(N,K)); %Uncorrelated noise

X=A*diag(sqrt(P))*sig+noise; %Generate data matrix
R=X*X'/K; %Spatial covariance matrix

%Partition size
n=3;
%ULA partition offset
espacing = 1;

k=2*pi*d;
m = size(R, 1);


% ESPRIT algorithm adapted from wikipedia
[E, ~] = eig(R);
Es = E(:,end - n + 1:end);
Es1 = Es(1:end - espacing,:);
Es2 = Es(espacing + 1:end,:);
P=Es1\Es2;
z = eig(P);
doa=-asin(angle(z)/(k*espacing));
x_est = sort(doa);
%aerr=((x_est'-aoas).^2)';
end
