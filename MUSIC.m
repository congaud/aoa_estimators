clear all
close all


%Angle of Arrivals
th1=40;
th2=10;
th3=60;
aoa=[th1 th2 th3];
aoas=[th1 th2 th3]*pi/180;
q=length(aoas); %Number of signals

%Power
p1=1;
p2=1;
p3=1;
P=1*[p1 p2 p3];

%Other Params
N=2000; %Number of samples
d=0.5;
p=5; %Number of antennas

noise_var=30;
fc=28e6; %Frequency in Hz
t=0:1e-9:N*1e-9-1e-9; %Time


%Steering vector - A is p*q, p antennas, q signals
A=exp(-1i*2*pi*d*(0:p-1)'*sin(aoas));


%Signal - S is q*N, q signals, N samples
sig=[sin(2*pi*fc*t);sin(2*pi*fc*t+2*pi/3);sin(2*pi*fc*t+2*pi/5)];
%sig=round(rand(q,N))*2-1;
%Noise - N is p*N, p antennas, N samples
noise=sqrt(noise_var/2)*(randn(p,N)+1i*randn(p,N));

%Recieved signals - X is p*N, p antennas, N samples
X=A*diag(sqrt(P))*sig+noise; %Generate data matrix

%Covariance Matrix
R=X*X'/N;

%Eigendecomposition
[Q ,D]=eig(R);

%Largest eigenvalues are signals, smallest noise
[D,I]=sort(diag(D),1,'descend'); %Find r largest eigenvalues

%Sort the eigenvectors classified by eigenvalues
Q=Q (:,I);
Qs=Q (:,1:q); 
Qn=Q(:,q+1:p);
angles=(-90:0.1:90);
a1=exp(-i*2*pi*d*(0:p-1)'*sin([angles(:).']*pi/180));
for k=1:length(angles)
    music_spectrum(k)=1/(a1(:,k)'*Qn*Qn'*a1(:,k));
end
figure
hold on
plot(angles,abs(music_spectrum).^2)
grid on
title('MUSIC Spectrum')
xlabel('Angle in degrees')

C=Qn*Qn';
coeff = zeros(p - 1, 1);
for index = 1:p-1
    coeff(index) = sum(diag(C, index));
end
D = [flipud(coeff); sum(diag(C)); conj(coeff)];
z = roots(D);
recieved_aoa = angle(z);

nz = length(z);
mask = true(nz, 1);

for ii = 1:nz
    absz = abs(z(ii));
    if absz > 1
        mask(ii) = false;
    elseif absz == 1
        % find the closest point and remove it
        idx = -1;
        dist = inf;
        for jj = 1:nz
            if jj ~= ii && mask(jj)
                cur_dist = abs(z(ii) - z(jj));
                if cur_dist < dist
                    dist = cur_dist;
                    idx = jj;
                end
            end
        end
        mask(idx) = false;
    end
end
z = z(mask);
[~, idx] = sort(1 - abs(z));
x_est = sort(-cm2doa(z(idx(1:q)), 2*pi*d, 'radian'));
plot(x_est,max(abs(music_spectrum).^2)/2.*[1,1,1],'x')
aerr=((x_est'-aoas).^2)';

function doa = cm2doa(z, k, unit)
%CM2DOA Converts complex exponentials to doas.
switch lower(unit)
    case 'radian'
        doa = asin(angle(z) / k);
    case 'degree'
        doa = rad2deg(asin(angle(z) / k));
    case 'sin'
        doa = angle(z) / k;
    otherwise
        error('Unkown unit ''%s''.', unit);
end
end
