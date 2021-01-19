close all
clear all

%Angle of Arrivals
th1=40;
th2=10;
th3=60;
aoas=[th1 th2 th3]*pi/180;
q=length(aoas); %Number of signals

%Power
p1=500;
p2=500;
p3=500;
P=[p1 p2 p3];

%Other Params
N=2000; %Number of samples
d=0.5;
p=5; %Number of antennas

noise_var=0;
fc=28*10^5; %Frequency in Hz
t=0:1e-9:N*1e-9-1e-9; %Time
 
%Steering Vector - A is p*q, p antennas, q signals 
A=exp(-i*2*pi*d*(0:p-1)'*sin([aoas(:).']));

%Signals - S is q*N, q signals, N samples 
S=round(rand(q,N))*2-1; % Generate random BPSK symbols for each of the r signals
%S=[sin(2*pi*fc*t); sin(2*pi*fc*t+2*pi/3); sin(2*pi*fc*t+2*pi/5)];

%Noise - N is p*N, p antennas, N samples
Noise=sqrt(noise_var/2)*(randn(p,N)+i*randn(p,N));

%Recieved signals - X is p*N, p antennas, N samples
X=A*S+Noise;

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

theta_bar=[];
J=ones(1,180);
for index=0:179
    J(1,1+index)=trace(proj_oper(steer_array(index*pi/180,p))*R);
end
[~,theta_list(1)]=max(J);
for index=2:q
    J=ones(1,180);
    for index1=0:179
        b=rescol_unitvector(index1*pi/180,theta_list.*pi./180,p);
        J(1,1+index1)=b'*R*b;
    end
    [~,theta_list(index)] = max(J);
end
theta_list_0=theta_list;

while 1
    theta_list_old=theta_list
    for index=1:q
        theta_list_redu=theta_list;
        theta_list_redu(index)=[];
        J=ones(1,180);
        for index1=0:179
           if isempty(theta_list(theta_list==index1))==true
               b=rescol_unitvector(index1*pi/180,theta_list_redu.*pi./180,p); 
           end
           J(1,1+index1)=b'*R*b;
        end
        [~,theta_list(index)] = max(J);
        
    end
    if norm(theta_list-theta_list_old)<1
        break
    end
end


function B = rescol_unitvector(k,A,p)
    numer = (eye(p)-proj_oper(steer_array(A,p)))*steer_array(k,p);
    denom = norm(numer);
    B=numer/denom;
    
end

function B = steer_array(A,p)
    d=0.1;
    B=exp(-1i*2*pi*d*(0:p-1)'*sin(A));
end
    
function B = proj_oper(A)
    B=A*inv(A'*A)*A';
end

