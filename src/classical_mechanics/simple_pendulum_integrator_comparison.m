%Simple pendulum integrator comparison

%Solves the simple pendulum equation: d^2{theta}/dt^2 = - sin(theta)
%using forward Euler, backward Euler, or Euler-Cromer method

%Timofey Golubev

%Choose the method to use
%0 = all three, 1 = forward Euler, 2 = backward Euler, 3 = Euler-Cromer
% clear vars;

method = 0;

%time step and range
t_max = 25;
dt = 0.01;
num_steps = t_max/dt + 1;  %+1 to include 0th step

%initial conditions, in degrees
theta_0 = 86;
theta_dot_0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize arrays
t = zeros(num_steps,1);
t(1) = 0;  %note: b/c matlab starts indices from 1, need to use this

if(method ==1 || method ==0)
    %forward Euler
    theta_1 = zeros(num_steps,1);
    theta_dot_1 = zeros(num_steps,1);
    %Initial conditions: convert to radians and fill in array
    theta_1(1) = degtorad(theta_0);
    theta_dot_1(1) = degtorad(theta_dot_0);
    for i = 1:num_steps-1
        t(i+1) = t(i) + dt;
        theta_1(i+1) = theta_1(i) + theta_dot_1(i)*dt;
        theta_dot_1(i+1) = theta_dot_1(i) -sin(theta_1(i))*dt;
    end  
end

if(method==2 || method ==0)
    %backward Euler
    theta_2 = zeros(num_steps,1);
    theta_dot_2 = zeros(num_steps,1);
    %Initial conditions: convert to radians and fill in array
    theta_2(1) = degtorad(theta_0);
    theta_dot_2(1) = degtorad(theta_dot_0);
    for i = 1:num_steps-1
        t(i+1) = t(i) + dt;
        theta_dot_2(i+1) = theta_dot_2(i) -(sin(theta_2(i)) + cos(theta_2(i))*dt)*dt;  %use sin(theta(i+1)) ~ sin(theta(i)) + cos(theta(i))*dt
        theta_2(i+1) = theta_2(i) + theta_dot_2(i+1)*dt;    
    end
end

if(method ==3 || method ==0)
    %Euler-Cromer
    theta_3 = zeros(num_steps,1);
    theta_dot_3 = zeros(num_steps,1);
    %Initial conditions: convert to radians and fill in array
    theta_3(1) = degtorad(theta_0);
    theta_dot_3(1) = degtorad(theta_dot_0);
     for i = 1:num_steps-1
        t(i+1) = t(i) + dt;
        theta_dot_3(i+1) = theta_dot_3(i) -sin(theta_3(i))*dt;
        theta_3(i+1) = theta_3(i) + theta_dot_3(i+1)*dt;  
     end
end

%Analytic solution
%from IC's: recall that index 1 corresponds to t=0
B = degtorad(theta_0);
A = degtorad(theta_dot_0);
%initialize arrays
theta_analytic = zeros(num_steps,1);
theta_dot_analytic = zeros(num_steps,1);
%Initial conditions: convert to radians and fill in array
theta_analytic(1) = degtorad(theta_0);
theta_dot_analytic(1) = degtorad(theta_dot_0);
for i = 1:num_steps
    theta_analytic(i) = A*sin(t(i)) + B*cos(t(i));
    theta_dot_analytic(i) = A*cos(t(i)) - B*sin(t(i));
end

%Make plots
  matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual') 
     
%theta vs. t
figure 
hold on
 if(method == 0)
     plot(t, theta_1,'b','LineWidth',1)
     plot(t, theta_2, 'r','LineWidth',1)
     plot(t, theta_3, 'Color',[0.93 0.69 0.125]','LineWidth',1) %gives orange color
 elseif(method ==1)
    plot(t, theta_1, 'b','LineWidth',1)
 elseif(method ==2)
     plot(t, theta_2, 'r','LineWidth',1)
 elseif(method ==3)
     plot(t, theta_3, 'Color',[0.93 0.69 0.125]','LineWidth',1)
 end
 plot(t, theta_analytic, 'k','LineWidth',1) %black line
 set(gca, 'FontSize', 20)
 xlabel('Time(s)','interpreter','latex','FontSize',22);
 ylabel({'$\theta$ (rad/s)'},'interpreter','latex','FontSize',22);

%theta_dot vs. t
figure 
 hold on
 if(method == 0)
     plot(t, theta_dot_1,'b','LineWidth',1)
     plot(t, theta_dot_2,'r','LineWidth',1)
     plot(t, theta_dot_3,'Color',[0.93 0.69 0.125]','LineWidth',1)
 elseif(method ==1)
    plot(t, theta_dot_1,'b','LineWidth',1)
 elseif(method ==2)
     plot(t, theta_dot_2,'r','LineWidth',1)
 elseif(method ==3)
     plot(t, theta_dot_3,'Color',[0.93 0.69 0.125]','LineWidth',1)
 end
 plot(t, theta_dot_analytic,'k','LineWidth',1)
 set(gca, 'FontSize', 20)
 xlabel('Time(s)','interpreter','latex','FontSize',22);
 ylabel({'$\dot{\theta}$ (rad/s)'},'interpreter','latex','FontSize',22);


%phase space
figure 
 hold on
 if(method==0)
     plot(theta_1, theta_dot_1,'b','LineWidth',1)
     plot(theta_2, theta_dot_2,'r','LineWidth',1)
     plot(theta_3, theta_dot_3,'Color',[0.93 0.69 0.125]','LineWidth',1)
 elseif(method ==1)
    plot(theta_1, theta_dot_1,'b','LineWidth',1)
 elseif(method ==2)
     plot(theta_2, theta_dot_2,'r','LineWidth',1)
 elseif(method ==3)
     plot(theta_3, theta_dot_3,'Color',[0.93 0.69 0.125]','LineWidth',1)
 end
 plot(theta_analytic, theta_dot_analytic,'k','LineWidth',1)
 set(gca, 'FontSize', 20)
 xlabel({'$\theta$(rad)'},'interpreter','latex','FontSize',22);
 ylabel({'$\dot{\theta}$(rad/s)'},'interpreter','latex','FontSize',22);

