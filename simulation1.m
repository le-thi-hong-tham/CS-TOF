clear all;
close all;
clc;

f=1000;
T=1/f;
fs=50;
Ts=1/fs;
k= 1; %sparse level per cycles
Nc = 100 % number sample per cycles
M = Nc*f/500; % number tranfers
N= Nc * f/fs; % length of signal

%generate signal reference
refsig = zeros(N,1);
tempref = zeros(Nc,1) % signal tranfers
tempref(1,1)= 1;

% time of the light wave flighting from the object to the imaging sensor
shiftime = 10;

%generate signal object
objsig = zeros(N,1);
tempobj = zeros(Nc,1)
tempobj (1+shiftime: k+shiftime)=1;

for i= 1:f/fs
    refsig((i-1)*Nc+1:i*Nc) = tempref(:,1);
    objsig((i-1)*Nc+1:i*Nc) = tempobj(:,1);
end

for i= 1:N
    t(i) = (i-1)/(f*Nc);
end

figure(1);
plot(t,refsig,'linewidth',2)
ylim([-0.2 1.2]);

Phi = randi([0 1],N,N);
y=Phi*refsig;
figure(2);
plot(y)


for i=1:M
   position(i,1) = (i-1) *shiftime+1
end

outputref = zeros(M,1);
outputobj = zeros(M,1);

%Making random measurements
A=zeros(M,N);
for i=1 : M
    outputref(i) = y(position(i));
    A(i,:) = Phi(position(i),:);
end
outputobj = A * objsig;


figure(3)
plot(outputref);
hold on
plot(outputobj);
title('Measurement vector');

% %Adding some measurement noise.
% SNR=15;
% outputobj=createNoise(outputobj,SNR);
% 
% 
% %Measurement vector with noise.
% outputobj = outputobj + e ;
% figure(4);
% plot(outputref);
% hold on
% plot(outputobj);
% title(sprintf('Measurement vector with noise at SNR=%d dB', SNR));



cvx_begin
    variable ref(N);
    minimize (norm(ref,1));
    subject to
    A*ref==outputref;
cvx_end

% obj = OrthogonalMatchingPursuit(A,20,outputobj);

cvx_begin
    variable obj(N);
    minimize (norm(obj,1));
    subject to
    A*obj==outputobj; 

%     norm(A*obj-outputobj,2) <= eps
%     minimize (norm(A*obj-outputobj,2)+0.01*norm(obj,1));
cvx_end

diff = objsig - obj;
recovery_error = norm(diff) / norm(objsig);
fprintf('recovery error: %0.4f\n', recovery_error);


figure(5)
plot(ref)
hold on 
plot(obj)
title('Reconstructed signal');
legend('ref','obj')

