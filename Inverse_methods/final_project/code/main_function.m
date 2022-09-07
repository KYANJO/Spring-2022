%%====================================================================
% Boise State University
% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
% main_fuction script
%%=====================================================================
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
        end
    
        noise = sig*randn(size(Q(:,:)));      % genrate noise
        d = Q(:,:) + noise;   d = d';         % add noise to the data
         
        dh = d(:,1); du = d(:,2); dhu = d(:,3);
        dnoise = [0.005*dh  0.005*du 0.005*dhu ];
        d = Q(:,:) + dnoise';
    
        d = d';                               % transposing data d
        delta = sig*sqrt(N);                  % delta
        fun = @(m) model(m);                  % model function handle
        jac = @(m) Jacobian(m);               % Jacobian function handle

        m = occam(fun, jac, L, d, mo, delta);    % calling the Occam model to 
                                                    % recover mtrue
        
        % plotting estimates
        figure(1)
        if mq == 1
            plot(m(:,mq),'k*'); hold on
            plot(mtrue(:,mq),'ro'); hold off
            legend('m_{estimate}','m_{true}',Location='best')
            title('Height field');
            ylim([hr-0.25,hl+0.25]);
            ylabel('height');xlabel('x')
            frame = getframe(gcf);
            writeVideo(v,frame);
        elseif mq == 2
            plot(m(:,mq),'k*'); hold on
            plot(mtrue(:,mq),'ro'); hold off
            legend('m_{estimate}','m_{true}',Location='best')
            title('velocity field')
            ylabel('velocity');xlabel('x')
            frame = getframe(gcf);
            writeVideo(v,frame);
        elseif mq == 3
            plot(m(:,mq),'k*'); hold on
            plot(mtrue(:,mq),'ro'); hold off
            legend('m_{estimate}','m_{true}',Location='best')
            title('momentum field')
            ylabel('Momentum');xlabel('x')
            frame = getframe(gcf);
            writeVideo(v,frame);
        else
            subplot(2,2,1)
            plot(m(:,1),'k*'); hold on
            plot(mtrue(:,1),'ro');  hold off
            legend('m_{estimate}','m_{true}',Location='best')
            title('Height estimates')
            ylim([hr-0.25,hl+0.25])
            ylabel('height');xlabel('x')
            subplot(2,2,2);
            plot(m(:,2),'k*'); hold on
            plot(mtrue(:,2),'ro'); hold off
            legend('m_{estimate}','m_{true}',Location='best')
            title('Velocity estimates')
            ylabel('velocity');xlabel('x')
            subplot(2,2,[3,4]);
            plot(m(:,3),'k*'); hold on
            plot(mtrue(:,3),'ro'); hold off
            title('Momentum estimates')
            ylabel('Momentum');xlabel('x')
            legend('m_{estimate}','m_{true}',Location='best')
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        
        for i=1:N 
            xi = x(i)/t(j);
            [Qest(:,i),hs,us] = forestclaw(m(1,:),m(3,:),xi);
        end
    
        % chi-square
        chi_s = zeros(3,1); ad = d';
        for k = 1:N
            chi_s = chi_s + (Qest(:,k) - ad(:,k)).^2;   
        end
        
        % pvalue
        dof = N - 9; %degrees of freedom
        p = 1 - chi2cdf(chi_s,3);
    
        % plotting data
        figure(2)
        if mq == 1
            plot(Qest(mq,:),'k'); hold on
            plot(Q(mq,:),'r'); hold off
            legend('h_{estimate}','h_{true}',Location='best')
            title('Height field');
            ylim([hr-0.5,hl+0.5]);
            ylabel('height');xlabel('x')
            frame = getframe(gcf);
            writeVideo(v,frame);
        elseif mq == 2
            plot(Qest(mq,:),'k'); hold on
            plot(Q(mq,:),'r'); hold off
            legend('u_{estimate}','u_{true}',Location='best')
            title('velocity field')
            ylabel('velocity');xlabel('x')
            frame = getframe(gcf);
            writeVideo(v,frame);
        elseif mq == 3
            plot(Qest(mq,:),'k'); hold on
            plot(Q(mq,:),'r'); hold off
            legend('hu_{estimate}','hu_{true}',Location='best')
            title('momentum field')
            ylabel('Momentum');xlabel('x')
            frame = getframe(gcf);
            writeVideo(v,frame);
        else
            subplot(2,2,1)
            plot(Qest(1,:),'k'); hold on
            plot(Q(1,:),'r'); hold off
            legend('h_{estimate}','h_{true}',Location='best')
            title('Height field')
            ylim([hr-0.25,hl+0.25])
            ylabel('height');xlabel('x')
            subplot(2,2,2);
            plot(Qest(2,:),'k'); hold on
            plot(Q(2,:),'r'); hold off
            legend('u_{estimate}','u_{true}',Location='best')
            title('Velocity field')
            ylabel('velocity');xlabel('x')
            subplot(2,2,[3,4]);
            plot(Qest(3,:),'k'); hold on
            plot(Q(3,:),'r'); hold off
            ylim([-1.1,0])
            title('Momentum field')
            ylabel('Momentum');xlabel('x')
            legend('hu_{estimate}','hu_{true}',Location='best')
            frame = getframe(gcf);
            writeVideo(v,frame);
            
        end
        
        
        %Covariance Matrix
        J = Jacobian(m); 
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
    L2norm = norm(mtrue - m,2)
    mtrue
    covariance_matrix = C
    cofidence_interval_height = [c1 c2]
    cofidence_interval_velocity = [c11 c22]
    cofidence_interval_Momentum = [c13 c33]
    
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
end
