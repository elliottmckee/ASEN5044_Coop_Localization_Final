% RBS Stat. Est. Final Project Part 2
% LKF


% F_tilde(tvec(1), dt)
%  x_perturb(:, 1)
 dx_tilde = zeros(length(x_perturb(:,1)),length(x_perturb(1,:)));
 dy_tilde = zeros(length(y_perturb(:,1)),length(y_perturb(1,:)));
 du_tilde = zeros(length(uVec(:,1)),length(uVec(1,:)));
 dx_tilde(:,1) = x_perturb(:,1);
 du_tilde(:,1) = uVec(:,1);
 dy(:,1) = y_perturb(:,1);
 % I think uVec is constant
 
 % Want to do a batch estimate of P and x0
 num = 100;
 HH = zeros(5,6);
 YY = zeros(5*num,1);
 RR = zeros(5*num);
 
 for idx = 1:num
     HH((idx-1)*5+1:5*idx,:) = H_tilde(x_nom(:,idx))*(F_tilde(tvec(idx), dt))^(idx-1);
     YY((idx-1)*5+1:idx*5) = ydata(:,idx);
     RR((idx-1)*5+1:idx*5,(idx-1)*5+1:idx*5) = Rtrue;
 end
     
 P0 = inv(HH'*inv(RR)*HH);
 x0 = inv(HH'*inv(RR)*HH)*HH'*inv(RR)*YY;
 
 
 P = zeros(6*length(x_perturb(1,:)),6);
 %P(1:6,:) = inv(H_tilde(x_nom(:,1))'*inv(Rtrue)*H_tilde(x_nom(:,1)));
 
 % K its a 6x5 matrix
 
 x_out = zeros(length(x_nom(:,1)),length(x_nom(1,:)));
 %x_out(:,1) = inv(H_tilde(x_nom(:,1))'*inv(Rtrue)*H_tilde(x_nom(:,1)))*H_tilde(x_nom(:,1))'*inv(Rtrue)*ydata(:,1);
 
 P(1:6,:) = P0;
 x_out(:,1) = x0;
 
 sig2 = zeros(6,length(x_nom(1,:)));
 sig2(:,1) = 2.*sqrt(diag(P(1:6,1:6)));
 % implimenting the kalman filter
 
 for idx = 2:length(x_nom(1,:))-1
     % dx tilde minus
     dx_tilde(:,idx) = F_tilde(tvec(idx-1), dt)*dx_tilde(:,idx-1);% + G_tilde(tvec(idx-1), dt)*uVec;
     %P minus at K+1
     P((idx-1)*6+1:6*idx,:) = F_tilde(tvec(idx-1), dt)*P((idx-2)*6+1:6*(idx-1),:)*F_tilde(tvec(idx-1), dt)';
     dy(:,idx) = ydata(:,idx) - h(x_nom(:,idx));
     K = P((idx-1)*6+1:6*idx,:)*H_tilde(x_nom(:,idx))'*inv(H_tilde(x_nom(:,idx))*P((idx-1)*6+1:6*idx,:)*H_tilde(x_nom(:,idx))' + Rtrue);
     P((idx-1)*6+1:6*idx,:) = (eye(6) - K*H_tilde(x_nom(:,idx)))*P((idx-1)*6+1:6*idx,:);
     % dx tilde plus
     dx_tilde(:,idx) = dx_tilde(:,idx) + K*(dy(:,idx) - H_tilde(x_nom(:,idx))*dx_tilde(:,idx));
     x_out(:,idx) = x_nom(:,idx) + dx_tilde(:,idx);
     sig2(:,idx) = 2.*sqrt(diag(P((idx-1)*6+1:6*idx,:)));
 end
 
 
figure()
hold on;
for k = 1:6
    ax(k) = subplot(6,1,k);
end


subplot(ax(1));
hold on;
title('LKF - State Estimation Errors')
%plot(tvec,x_out(1,:));
plot(tvec,x_out(1,:)-x_gt(1,:));
plot(tvec,x_out(1,:)+sig2(1,:),'-- r');
plot(tvec,x_out(1,:)-sig2(1,:),'-- r');
xlabel('time (s)');
ylabel('$e_{\xi_g}$','Interpreter','latex');

subplot(ax(2));
hold on;
%plot(tvec,x_out(2,:));
plot(tvec,x_out(2,:)-x_gt(2,:));
plot(tvec,x_out(2,:)+sig2(2,:),'-- r');
plot(tvec,x_out(2,:)-sig2(2,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_g}$','Interpreter','latex');

subplot(ax(3));hold on;
% plot(tvec,x_out(3,:));
plot(tvec,x_out(3,:)-x_gt(3,:));
plot(tvec,x_out(3,:)+sig2(3,:),'-- r');
plot(tvec,x_out(3,:)-sig2(3,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_g}$','Interpreter','latex');

subplot(ax(4));hold on;
% plot(tvec,x_out(4,:));
plot(tvec,x_out(4,:)-x_gt(4,:));
plot(tvec,x_out(4,:)+sig2(4,:),'-- r');
plot(tvec,x_out(4,:)-sig2(4,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\xi_a}$','Interpreter','latex');

subplot(ax(5));hold on;
% plot(tvec,x_out(5,:));
plot(tvec,x_out(5,:)-x_gt(5,:));
plot(tvec,x_out(5,:)+sig2(5,:),'-- r');
plot(tvec,x_out(5,:)-sig2(5,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\eta_a}$','Interpreter','latex');

subplot(ax(6));hold on;
% plot(tvec,x_out(6,:));
plot(tvec,x_out(6,:)-x_gt(6,:));
plot(tvec,x_out(6,:)+sig2(6,:),'-- r');
plot(tvec,x_out(6,:)-sig2(6,:),'-- r');
xlabel('time (s)')
ylabel('$e_{\theta_a}$','Interpreter','latex');
 

figure()
hold on;
plot(tvec,sig2(1,:),'-- r');
plot(tvec,-sig2(1,:),'-- r');
plot(tvec,sig2(2,:),'-- b');
plot(tvec,-sig2(2,:),'-- b');
plot(tvec,sig2(3,:),'-- g');
plot(tvec,-sig2(3,:),'-- g');
plot(tvec,sig2(4,:),'-- c');
plot(tvec,-sig2(4,:),'-- c');
plot(tvec,sig2(5,:),'-- k');
plot(tvec,-sig2(5,:),'-- k');
plot(tvec,sig2(6,:),'-- m');
plot(tvec,-sig2(6,:),'-- m');