clear all;
close all;
clc;

f=1000; %frequency of light wave
T=1/f; %duty cycle of light wave
fs=50; % rate of camera sensor
Ts=1/fs; %duty cycle of camera sensor
k= 1; %sparse level per cycles
Nc = 100 % number sample per cycles
M = Nc*f/500; % number tranfers - measuments
N= Nc * f/fs; % length of signal

%generate signal reference
refsig = zeros(N,1);
ref= zeros(Nc,1) % signal tranfers per cycle
ref(1,1)= 1;

% time of the light wave flighting from the object to the imaging sensor
shiftime = 10;

for i= 1:f/fs
    refsig((i-1)*Nc+1:i*Nc) = ref(:,1);
end


figure(1);
% plot(t,refsig,'linewidth',2)
plot(refsig);
hold on
plot(objsig)
ylim([-0.2 1.2]);
xlabel('sample');
ylabel('Amplitude');

%generate encode signal
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
    A*xp_ref==outputref;
cvx_end

% obj = OrthogonalMatchingPursuit(A,20,outputobj);

cvx_begin
    variable xp_obj(N);
    minimize (norm(xp_obj,1));
    subject to
    A*xp_obj==outputobj; 

%     norm(A*obj-outputobj,2) <= eps
%     minimize (norm(A*obj-outputobj,2)+0.01*norm(obj,1));
cvx_end

%Compute error recovered
diff_ref = refsig - xp_ref;
recovery_error_ref = norm(diff_ref) / norm(refsig);
fprintf('recovery error: %0.4f\n', recovery_error_ref);

diff = objsig - xp_obj;
recovery_error = norm(diff) / norm(objsig);
fprintf('recovery error: %0.4f\n', recovery_error);


figure(5)
plot(xp_ref)
hold on 
plot(xp_obj)
xlabel('sample');
ylabel('Amplitude');
title('Reconstructed signal');
legend('ref','obj')
