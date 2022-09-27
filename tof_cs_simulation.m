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

% generate signal reference
refsig = zeros(N,1);
ref= zeros(Nc,1) % signal tranfers
ref(1,1)= 1;

% time of the light wave flighting from the object to the imaging sensor
shiftime = 10;

% generate signal object
objsig = zeros(N,1);
obj = zeros(Nc,1)
obj (1+shiftime: k+shiftime)=1;

for i= 1:f/fs
    refsig((i-1)*Nc+1:i*Nc) = ref(:,1);
    
end

for i= 1:N
    t(i) = (i-1)/(f*Nc);
end

figure(1);
plot(t,refsig,'linewidth',2)
plot(refsig);
ylim([-0.2 1.2]);
xlabel('sample');
ylabel('Amplitude');

% generate encode signal
Phi = randi([0 1],N,N);
y=Phi*refsig;



figure(2);
plot(y)



Tk_ref=1;
Tk_obj=10;
for i=1:M
   position_ref(i,1) = (i-1) *shiftime+ Tk_ref
   position_obj(i,1) = (i-1) *shiftime+ Tk_obj
   if (position_ref(i) > N )
       position_ref(i) = position_ref(i)-N
   end 
   if (position_obj(i) > N )
       position_obj(i) = position_obj(i)-N
   end
   
end

outputref = zeros(M,1);
outputobj = zeros(M,1);

%Making random measurements
A1=zeros(M,N);
A2=zeros(M,N);

    

for i=1 : M
    outputref(i) = y(position_ref(i));
    A1(i,:) = Phi(position_ref(i),:);
    outputobj(i) = y(position_obj(i));
    A2(i,:) = Phi(position_obj(i),:);
end

tof=10; 

outputobj(M+1:M+tof)= 0;
outputobj = circshift(outputobj,tof)

A3=zeros(size(A2,1)+tof,size(A2,2)+tof);
A3(1+tof:end,1+tof:end)= A2
% 
figure(3)
plot(outputref);
hold on
plot(outputobj);
xlabel('sample');
ylabel('Intensity');
title('Measurement vector');
legend('ref','obj')

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
    variable xp_ref(N);
    minimize (norm(xp_ref,1));
    subject to
    A1*xp_ref==outputref;
cvx_end

% obj = OrthogonalMatchingPursuit(A,20,outputobj);

cvx_begin
    variable xp_obj(N+tof);
    minimize (norm(xp_obj,1));
    subject to
    A3*xp_obj==outputobj; 

%     norm(A*obj-outputobj,2) <= eps
%     minimize (norm(A*obj-outputobj,2)+0.01*norm(obj,1));
cvx_end

% diff = refsig - xp_obj;
% recovery_error = norm(diff) / norm(refsig);
% fprintf('recovery error: %0.4f\n', recovery_error);


figure(5)
plot(xp_ref)
hold on
plot(xp_obj)
xlabel('sample');
ylabel('Amplitude');
title('Reconstructed signal');
legend('ref','obj')

