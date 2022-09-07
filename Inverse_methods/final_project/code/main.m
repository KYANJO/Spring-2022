%%====================================================================
% Boise State University
% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
% main script
%%=====================================================================

clc
clear
close all

global g; global N;
global sig; global h;
global t; global x;
global M;
warning('off','all')

% problem
hl = 2;                          % left depth
hr = 1;                          % right depth
ul = -0.5;                       % left velocity
ur = -0.5;                       % right velocity

% Spatial domain
ax = -5;
bx = 5;
ay = -2;
by = 4;
meqn = 2;                       % Number of equations in the system


g = 1;                          % Gravity 

to = 0;                         % initial time
Tfinal = 1;                     % final time

ql = [hl; ul; hl*ul];           % left conservation variable
qr = [hr; ur; hr*ur];           % right consrvation variable
qm = [(hl+hr)/2 (ul+ur)/2 ...   % intermediate initial state
    (hl*ul+hr*ur)/2];

N = 64;                         % Number of spartial steps

dx = (bx - ax)/N;               % spartial step size

cfl = 0.9;                      % cfl number

a = 1.5;                        % maximum velocity

dt_est = cfl*dx/a;

M = (floor(Tfinal/dt_est) + 1); % number of time steps

dt = Tfinal/(M);                % temporal step size

t = linspace(to,Tfinal,M);      % temporal domain

xe = linspace(ax,bx,N+1);       % edge locations
x = xe(1:end-1) + dx/2;         % Cell-center locations

mq = 4;                         % 1-height, 
                                % 2-velocity, 
                                % 3-momentum, 
                                % 4-all fields
h = 0.9;

sig = 1e-3;                     % standard deviation

L0 = eye(3);          % generate roughening matrices
L1 = get_l_rough(3,1);
L2 = get_l_rough(3,2); 


qll = [hl ul hl*ul];            % left conservation variable
qrr = [hr ur hr*ur];            % right conservation variable
xi = x(1)/t(2);                 % initial speed
[Q(:,1),hs,us] = forestclaw(ql,qr,xi); % initial conservation variable
qmm = [hs us hs*us];                   % intermediate state
%mo = [qll;qm;qrr];
mtrue = [qll; qmm; qrr];               % true parameters

mo = 1.1*mtrue;
main_function(mtrue,mo,h,L0,x,ql,qr,M,t,mq)


function main_function(mtrue,mo,h,L,x,ql,qr,M,t,mq)
    global N; global sig;
    hr = qr(1); hl = ql(1)
    % Initial conditions
    Q = zeros(3,N);
    for i=1:N
        if (x(i)<=0) 
            Q(:,i) = ql;
        else
            Q(:,i) = qr;
        end
    end
    
    % writing a video
    v = VideoWriter('dam.avi');
    open(v);
    
    for j= 2:M % time loop
        J = [];
        for i=1:N % spartial loop
            xi = x(i)/t(j);
            [Q(:,i),hs,us] = forestclaw(ql,qr,xi);
            
            % Formulation of Jacobian Matrix
            ej = ones(3,1); 
            [Qmin(:,i),hs,us] = forestclaw(ql-h*ej,qr-h*ej,xi);
            [Qmax(:,i),hs,us] = forestclaw(ql+h*ej,qr+h*ej,xi);
            J = [J (Qmax(:,i) - Qmin(:,i))./(2*h)]; % Jacobian
            
        end
    
        noise = sig*randn(size(Q(:,:)));      % genrate noise
        d = Q(:,:) + noise;   d = d';         % add noise to the data
         
        dh = d(:,1); du = d(:,2); dhu = d(:,3);
        dnoise = [0.005*dh  0.005*du 0.005*dhu ];
        d = Q(:,:) + dnoise';
    
        G = Q'; J = J';  d = d';              % transposing variables G,d,and J
        delta = sig*sqrt(N);                  % delta

        m = occam12(G, J, L, d, mtrue, delta);    % calling the Occam model to 
                                                   % recover mtrue
         m1 = occam12(G, J, L, d, 1.2*mtrue, delta);
         m2 = occam12(G, J, L, d, 1.8*mtrue, delta);
%         m = occam( M,N,h,x,t, L, d, mtrue, delta)
        
        % plotting estimates
%         figure(1)
%         if mq == 1
%             plot(m(:,mq),'k*'); hold on
%             plot(mtrue(:,mq),'ro'); hold off
%             legend('m_{estimate}','m_{true}',Location='best')
%             title('Height field');
%             ylim([hr-0.25,hl+0.25]);
%             ylabel('height');xlabel('x')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         elseif mq == 2
%             plot(m(:,mq),'k*'); hold on
%             plot(mtrue(:,mq),'ro'); hold off
%             legend('m_{estimate}','m_{true}',Location='best')
%             title('velocity field')
%             ylabel('velocity');xlabel('x')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         elseif mq == 3
%             plot(m(:,mq),'k*'); hold on
%             plot(mtrue(:,mq),'ro'); hold off
%             legend('m_{estimate}','m_{true}',Location='best')
%             title('momentum field')
%             ylabel('Momentum');xlabel('x')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         else
%             subplot(2,2,1)
%             plot(m(:,1),'k*'); hold on
%             plot(mtrue(:,1),'ro');  hold off
%             legend('m_{estimate}','m_{true}',Location='best')
%             title('Height estimates')
%             ylim([hr-0.25,hl+0.25])
%             ylabel('height');xlabel('x')
%             subplot(2,2,2);
%             plot(m(:,2),'k*'); hold on
%             plot(mtrue(:,2),'ro'); hold off
%             legend('m_{estimate}','m_{true}',Location='best')
%             title('Velocity estimates')
%             ylabel('velocity');xlabel('x')
%             subplot(2,2,[3,4]);
%             plot(m(:,3),'k*'); hold on
%             plot(mtrue(:,3),'ro'); hold off
%             title('Momentum estimates')
%             ylabel('Momentum');xlabel('x')
%             legend('m_{estimate}','m_{true}',Location='best')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         end
%         
        for i=1:N 
            xi = x(i)/t(j);
            [Qest(:,i),hs,us] = forestclaw(m(1,:),m(3,:),xi);
            [Qest1(:,i),hs,us] = forestclaw(m1(1,:),m1(3,:),xi);
            [Qest2(:,i),hs,us] = forestclaw(m2(1,:),m2(3,:),xi);
        end
    
        % chi-square
        chi_s = zeros(3,1); ad = d';
        chi_s1 = zeros(3,1);chi_s2 = zeros(3,1);
        for k = 1:N
            chi_s = chi_s + (Qest(:,k) - ad(:,k)).^2;   
            chi_s1 = chi_s1 + (Qest1(:,k) - ad(:,k)).^2;   
            chi_s2 = chi_s2 + (Qest2(:,k) - ad(:,k)).^2;   
        end
        
        % pvalue
        dof = N - 9; %degrees of freedom
        p = 1 - chi2cdf(chi_s,3);

        p1 = 1 - chi2cdf(chi_s1,3);
        p2 = 1 - chi2cdf(chi_s2,3);
%     
%     %     % plotting data
%         figure(2)
%         if mq == 1
%             plot(Qest(mq,:),'k'); hold on
%             plot(Q(mq,:),'r'); hold off
%             legend('h_{estimate}','h_{true}',Location='best')
%             title('Height field');
%             ylim([hr-0.5,hl+0.5]);
%             ylabel('height');xlabel('x')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         elseif mq == 2
%             plot(Qest(mq,:),'k'); hold on
%             plot(Q(mq,:),'r'); hold off
%             legend('u_{estimate}','u_{true}',Location='best')
%             title('velocity field')
%             ylabel('velocity');xlabel('x')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         elseif mq == 3
%             plot(Qest(mq,:),'k'); hold on
%             plot(Q(mq,:),'r'); hold off
%             legend('hu_{estimate}','hu_{true}',Location='best')
%             title('momentum field')
%             ylabel('Momentum');xlabel('x')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         else
%             subplot(2,2,1)
%             plot(Qest(1,:),'k'); hold on
%             plot(Q(1,:),'r'); hold off
%             legend('h_{estimate}','h_{true}',Location='best')
%             title('Height field')
%             ylim([hr-0.25,hl+0.25])
%             ylabel('height');xlabel('x')
%             subplot(2,2,2);
%             plot(Qest(2,:),'k'); hold on
%             plot(Q(2,:),'r'); hold off
%             legend('u_{estimate}','u_{true}',Location='best')
%             title('Velocity field')
%             ylabel('velocity');xlabel('x')
%             subplot(2,2,[3,4]);
%             plot(Qest(3,:),'k'); hold on
%             plot(Q(3,:),'r'); hold off
%             ylim([-1.1,0])
%             title('Momentum field')
%             ylabel('Momentum');xlabel('x')
%             legend('hu_{estimate}','hu_{true}',Location='best')
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%             
%         end
        
        
        %Covariance Matrix
        C = inv(J'*J);
    
        %confidence interval
        za = 1.96; % 95% confidence interval
        
        % first parameter
        s1 = sqrt(C(1,1)); % standard deviation
        s2 = sqrt(C(2,2));
        s3 = sqrt(C(3,3));
        
        % confidence intervals 
        c1 = m(:,1)- za*s1;
        c2 = m(:,1) + za*s1;
        c11 = m(:,2)- za*s2;
        c22 = m(:,2) + za*s2;
        c13 = m(:,3)- za*s3;
        c33 = m(:,3) + za*s3;
    
        %Correlation matrix
        rho11 = C(1,1)/sqrt(C(1,1)*C(1,1));
        rho22 = C(2,2)/sqrt(C(2,2)*C(2,2));
        rho33 = C(3,3)/sqrt(C(3,3)*C(3,3));
        rho12 = C(1,2)/sqrt(C(1,1)*C(2,2));
        rho13 = C(1,3)/sqrt(C(1,1)*C(3,3));
        rho23 = C(2,3)/sqrt(C(2,2)*C(3,3));
        Correlation_matrix = [rho11 rho12 rho13;rho12 rho22 rho23;...
            rho13 rho23 rho33];
        
        
    end
    close(v);
    
    disp(['chi-square obs = [',num2str(chi_s'),']'])
    disp(['pvalue = [',num2str(p'),']'])
    disp(['pvalue1 = [',num2str(p1'),']'])
    disp(['pvalue3 = [',num2str(p2'),']'])
    L2norm = norm(mtrue - m,2)
    mtrue
    covariance_matrix = C
    % cofidence_interval_height = [c1 c2]
    % cofidence_interval_velocity = [c11 c22]
    % cofidence_interval_Momentum = [c13 c33]
    
    Correlation_matrix = Correlation_matrix
    
    
    %Linearised ellipsoid
    Delta = chi2inv(0.95,3); %Delta2
    
    figure(4)
    subplot(2,2,1)
    plot_ellipse(Delta,C(1:2,1:2),m,1,2); 
    grid on
    title('Error Ellipsoid for h and u')
    legend('','left','','middle','','right',Location='best')
    xlabel('h'); ylabel('u')
    subplot(2,2,2)
    C1 = [C(2,2) C(2,3); C(3,2) C(3,3)];
    plot_ellipse(Delta,C1,m,2,3); 
    grid on
    title('Error Ellipsoid for u and hu')
    legend('','left','','middle','','right',Location='best')
    xlabel('u'); ylabel('hu')
    subplot(2,2,[3,4]);
    C2 = [C(1,1) C(1,3); C(3,1) C(3,3)];
    plot_ellipse(Delta,C2,m,1,3); 
    grid on
    title('Error Ellipsoid for h and hu')
    legend('','left','','middle','','right',Location='best')
    xlabel('h'); ylabel('hu')
    
    Corr = corrcoef(C);

    figure(7)
    clf
    colormap('gray')
    
    imagesc(Corr)
    set(colorbar,'Fontsize',18);
    xlabel('j')
    ylabel('i')
    title('Correlation Matrix')

    figure(9)
    subplot(2,2,1)
    plot(m(:,1),'b*'); hold on
    plot(m1(:,1),'k*');
    plot(m2(:,1),'g*');
    plot(mtrue(:,1),'ro');  hold off
    legend('mo = m_{true}','mo = 1.2m_{true}','mo = 1.8m_{true}','mtrue',Location='best')
    title('Height estimates')
%     ylim([hr-0.25,hl+0.25])
    ylabel('height');xlabel('x')
    subplot(2,2,2);
    plot(m(:,2),'b*'); hold on
    plot(m1(:,2),'k*');
    plot(m2(:,2),'g*');
    plot(mtrue(:,2),'ro'); hold off
    legend('mo = m_{true}','mo = 1.2m_{true}','mo = 1.8m_{true}','mtrue',Location='best')
    title('Velocity estimates')
    ylabel('velocity');xlabel('x')
    subplot(2,2,[3,4]);
    plot(m(:,3),'b*'); hold on
    plot(m1(:,3),'k*');
    plot(m2(:,3),'g*');
    plot(mtrue(:,3),'ro'); hold off
    title('Momentum estimates')
    ylabel('Momentum');xlabel('x')
    legend('mo = m_{true}','mo = 1.2*m_{true}','mo = 1.8m_{true}','mtrue',Location='best')

    figure(10)
    subplot(2,2,1)
    plot(Q(1,:),'-b'); hold on
    plot(Qest1(1,:),'-k');
    plot(Qest2(1,:),'-g');
      hold off
    legend('mo = m_{true}','mo = 1.2m_{true}','mo = 1.8m_{true}','mtrue',Location='best')
    title('Height estimates')
%     ylim([hr-0.25,hl+0.25])
    ylabel('height');xlabel('x')
    subplot(2,2,2);
    plot(Q(2,:),'-b'); hold on
    plot(Qest1(2,:),'-k');
    plot(Qest2(2,:),'-g');
    hold off
    legend('mo = m_{true}','mo = 1.2m_{true}','mo = 1.8m_{true}','mtrue',Location='best')
    title('Velocity estimates')
    ylabel('velocity');xlabel('x')
    subplot(2,2,[3,4]);
    plot(Q(3,:),'-b'); hold on
    plot(Qest1(3,:),'-k');
    plot(Qest2(3,:),'-g');
     hold off
    title('Momentum estimates')
    ylabel('Momentum');xlabel('x')
    legend('mo = m_{true}','mo = 1.2*m_{true}','mo = 1.8m_{true}','mtrue',Location='best')
end

