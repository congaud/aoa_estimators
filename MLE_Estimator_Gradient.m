close all
clear all

%Plots variance
x_hat_var=zeros(3,70);
for index0=1:70
    index0
    x_hat=zeros(3,20);
    for run=1:20
        SNR=10^((index0-20)/20);
        x_hat(:,run)=func(SNR);
    end
    x_hat_mean = sum(x_hat,2)./20;
    temp = (x_hat-x_hat_mean).^2;
    x_hat_var(:,index0) = sum(temp,2)./20;
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
dB=([1:69]-20)./20;
figure
plot(dB,x_hat_var(:,1:69))
xlabel('SNR (dB)')
ylabel('Variance')
title('Variance vs SNR')

function x_est=func(SNR)
%Angle of Arrivals
th1=10;
th2=15;
th3=20;
aoas=[th1 th2 th3];
q=length(aoas); %Number of signals

%Power
p1=1;
p2=1;
p3=1;
P=SNR*[p1 p2 p3];

%Other Params
N=2000; %Number of samples
d=0.5;
p=5; %Number of antennas

noise_var=1;
fc=28*10^5; %Frequency in Hz
t=0:1e-9:N*1e-9-1e-9; %Time
 
%Steering Vector - A is p*q, p antennas, q signals 
A=exp(-i*2*pi*d*(0:p-1)'*sind([aoas(:).']));

%Signals - S is q*N, q signals, N samples 
S=round(rand(q,N))*2-1; % Generate random BPSK symbols for each of the r signals
%S=[sin(2*pi*fc*t); sin(2*pi*fc*t+2*pi/3); sin(2*pi*fc*t+2*pi/5)];

%Noise - N is p*N, p antennas, N samples
Noise=sqrt(noise_var/2)*(randn(p,N)+i*randn(p,N));

%Recieved signals - X is p*N, p antennas, N samples
X=A*diag(sqrt(P))*S+Noise;

%Covariance matrix
R=0;
for index=1:N
    R=R+X(:,index)*X(:,index)';
end
R=R/N;
%figure
%hold on
%plot(t,S)
%plot(t,X)
%hold off

w0=[10;15;20]+10*randn(3,1);
%Weights
w0=rand(3,1)*90;

%Model
f = @(w)g(w,R,d);
[W,fW] = gradient_descent(f,w0,0.01,3000);
x_est=sort(W(:,end)).*pi/180;
%fW(end)
end

function [W,fW] = gradient_descent(f,w0,alpha,n_iter)

k = 1;
W = w0;
fW = f(w0);

while k < n_iter
    grad = approx_grad(f,W(:,k),.00001);
    W(:,k+1) = W(:,k) - alpha*(grad')/norm(grad);
    fW(k+1) = f(W(:,k+1));
    k=k+1;
end
end

function grad = approx_grad(f,w0,delta)
%Gradient w.r.t. theta
N = length(w0); %dimension
dw1 = delta*[1 0 0]';
dw2 = delta*[0 1 0]';
dw3 = delta*[0 0 1]';
grad1 = (f(w0+dw1)-f(w0-dw1))/2/delta;
grad2 = (f(w0+dw2)-f(w0-dw2))/2/delta;
grad3 = (f(w0+dw3)-f(w0-dw3))/2/delta;
grad=[grad1 grad2 grad3];
end

function y =g(w,R,d)
%A(theta)
A = exp(-1i*2*pi*d*(0:4)'*sind(w'));
temp=eye(size(A'*A));
%Better than matlab inv
invA=A'*A\temp;
%Projection matrix
P=A*invA*A';
%Cost function
L = P*R;
y=-1.*trace(L);
end

