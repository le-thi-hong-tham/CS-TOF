clear all;
close all;

% clc;
rng('default');
N=1000;
M=100;
K=9;
% t=randperm(N,K);
x=zeros(N,1);

for i=496:504
    x(i)=1;
end

% plot(x,'Linewidth',2);
% ylim([-0.2 1.2]);
figure(1)
stem(x, '.');
ylim([-0.2 1.2]);
xlabel('Index');
ylabel('Value');
title('Sparse vector');
find(x~=0)


% Constructing a Gaussian sensing matrix

Phi=randn(M,N);
% Make sure that variance is 1/sqrt(M)
Phi = Phi ./ sqrt(M);
%Computing norm of each column
column_norms = sqrt(sum(Phi .* conj(Phi)));
h = sqrt(sum(Phi .* conj(Phi), 1));
figure(2);
hist(h, 20);
xlabel('Norm');
ylabel('Count');

% Constructing a Gaussian dictionary with normalized columns
for i=1:N
    v = column_norms(i);
    % Scale it down
    Phi(:, i) = Phi(:, i) / v;
end
% 
% % Normalized dictionary
Phi =normalize_l2(Phi);
figure(3);
imagesc(Phi) ;
colormap(gray);
colorbar;
axis image;
title('\Phi');

% Making random measurements
% y0=Phi*x;
% stem(y0, '.');
% xlabel('Index');
% ylabel('Value');
% title('Measurement vector');
% 
% 
% %Adding some measurement noise.
% SNR=15;
% e=createNoise(y0,SNR);
% 
% % Measurement vector with noise.
% y=y0+e;
% figure(1);
% stem(y, '.');
% xlabel('Index');
% ylabel('Value');
% title(sprintf('Measurement vector with noise at SNR=%d dB', SNR));
% 
% 
% y1 = zeros ( numel (y),1);
% 
% t=4;
% y1(1+t:M+t)=y(1:M);
% figure(1), hold on ,
% stem(y1,'r.');
% 
% 
% 
% 
% col = zeros(M,t);
% 
% Phi = [col Phi]
% row = zeros(t,N+t);
% Phi=[row;Phi]
% 
% 
%           
% X_hat = OrthogonalMatchingPursuit(Phi,K,y1);
% 
% 
% 
% % cvx_begin
% % variable X_hat(N+t) ;
% % minimize norm(X_hat,1)
% % subject to
% %     Phi*X_hat - y1 <= eps;
% %     cvx_end;
% % % X_hat = OMPnorm(Phi,y1,floor(M/4),0);
% % cvx_begin
% %     variable X_hat(N+t);
% % %     minimize(norm(X_hat,1));
% %     
% %     subject to
% %         norm(Phi*X_hat-y1,2) <= eps
% %         minimize (norm(Phi*X_hat-y1,2)+0.01*norm(X_hat,1));
% % cvx_end
% % X_hat=algo_omp(K,Phi,y1);
%     figure(3)
% %     plot(x),hold on,
%     stem(X_hat,'.')
%     ylim([-0.1 1.2]);
%     
%     
% X=X_hat(5:end,:);
% diff=x-X;
% recovery_error = norm(diff) / norm(x);
% fprintf('recovery error: %0.4f\n', recovery_error);
% 
% 
% 
% 
% 
% 
% 
% 
% 
