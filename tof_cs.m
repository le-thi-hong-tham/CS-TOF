clear all;
close all;
clc;

f=1000; %frequency of light wave
f1=1000;
f2=1250;
T=1/f; %duty cycle of light wave
fs=50; % rate of camera sensor
Ts=1/fs; %duty cycle of camera sensor
k= 1; %sparse level per cycles
Nc1 = 100 % number sample per cycles
Nc2 = 80;
M = Nc1 * (f1+f2)/500; % number tranfers - measuments
N = Nc1 * f1/fs; % length of signal

%generate signal reference
refsig1 = zeros(N,1);
ref1= zeros(Nc1,1); % signal tranfers per cycle
ref1(1,1)= 1;

refsig2 = zeros(N,1);
ref2= zeros(Nc2,1); % signal tranfers per cycle
ref2(50,1)= 2;
% time of the light wave flighting from the object to the imaging sensor
shiftime = round(N/M);% plot(refsig);

objsig1 = zeros(N,1);
obj1= zeros(Nc1,1); % signal tranfers per cycle
obj1(1+shiftime: 1+shiftime)= 1;

objsig2 = zeros(N,1);
obj2= zeros(Nc2,1); % signal tranfers per cycle
obj2(50+shiftime: 50+shiftime)= 2;

% hold on
% plot(objsig);

% %generate signal object
% objsig = zeros(N,1);
% obj = zeros(Nc,1);
% obj (1+shiftime: 1+shiftime)=1; 
% obj (3+shiftime,1)=2; 
for i= 1:f1/fs
    refsig1((i-1)*Nc1+1:i*Nc1) = ref1(:,1);
    objsig1((i-1)*Nc1+1:i*Nc1) = obj1(:,1);
end
for i= 1:f2/fs
    refsig2((i-1)*Nc2+1:i*Nc2) = ref2(:,1);
    objsig2((i-1)*Nc2+1:i*Nc2) = obj2(:,1);
end

refsig = refsig1+refsig2;
objsig = objsig1+objsig2;
% for i= 1:N
%     t(i) = (i-1)/((f1+f2)*Nc);
% end

figure(1);
plot(refsig)
% plot(refsig);
hold on
plot(objsig);
ylim([-0.2 4]);
xlabel('ms');
ylabel('Amplitude');
title('Reconstructed signal');
legend('ref','obj')


% generate encode signal
Phi = randi([0 1],N,N);
y=Phi*refsig;
y1=Phi*objsig;
figure(2);
plot(y)
hold on
plot(y1);

for i=1:M
   position(i,1) = (i-1) *shiftime+1;
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
% 
% % %Adding some measurement noise.
% % SNR=15;
% % outputobj=createNoise(outputobj,SNR);
% % 
% % 
% % %Measurement vector with noise.
% % outputobj = outputobj + e ;
% % figure(4);
% % plot(outputref);
% % hold on
% % plot(outputobj);
% % title(sprintf('Measurement vector with noise at SNR=%d dB', SNR));
% 
% 
cvx_begin
    variable xp_ref(N);
    minimize (norm(xp_ref,1));
    subject to
    A*xp_ref==outputref;
cvx_end
% 
% % obj = OrthogonalMatchingPursuit(A,20,outputobj);
% PCR = zeros(2,1); %Reference phase element array
% PCO = zeros(2,1); %Object phase element array
cvx_begin
    variable xp_obj(N);
    minimize (norm(xp_obj,1));
    subject to
    A*xp_obj==outputobj; 

%     norm(A*obj-outputobj,2) <= eps
%     minimize (norm(A*obj-outputobj,2)+0.01*norm(obj,1));
cvx_end
% 
%Compute error recovered
diff_ref = refsig - xp_ref;
recovery_error_ref = norm(diff_ref) / norm(refsig);
fprintf('recovery error: %0.4f\n', recovery_error_ref);

diff = objsig - xp_obj;
recovery_error = norm(diff) / norm(objsig);
fprintf('recovery error: %0.4f\n', recovery_error);
% 
% 
figure(5)
plot(xp_ref)
hold on 
plot(xp_obj)
xlabel('sample');
ylabel('Amplitude');
title('Reconstructed signal');
legend('ref','obj');