%Simple pendulum integrator comparison

%Solves the simple pendulum equation: d^2{theta}/dt^2 = - sin(theta)
%using forward Euler, backward Euler, or Euler-Cromer method

%Timofey Golubev

%Choose the method to use
%0 = all three, 1 = forward Euler, 2 = backward Euler, 3 = Euler-Cromer
method = 3;

%time step and range
t_max = 10;
dt = 0.001;
num_steps = t_max/dt + 1;  %+1 to include 0th step

%initialize arrays
t = zeros(num_steps,1);

%initial conditions, in degrees
theta(1) = 15;
theta_dot(1) = 0;
t(1) = 0;  %note: b/c matlab starts indices from 1, need to use this

if(method ==1 || method ==0)
    %forward Euler
    theta_1 = zeros(num_steps,1);
    theta_dot_1 = zeros(num_steps,1);
    for i = 1:num_steps
        t(i+1) = t(i) + dt;
        theta_1(i+1) = theta_1(i) + theta_dot_1(i)*dt;
        theta_dot_1(i+1) = theta_dot_1(i) -sin(theta_1(i))*dt;
    end  
end

if(method==2 || method ==0)
    %backward Euler
    theta_2 = zeros(num_steps,1);
    theta_dot_2 = zeros(num_steps,1);
    for i = 1:num_steps
        t(i+1) = t(i) + dt;
        theta_dot_2(i+1) = theta_dot_2(i) -sin(theta_2(i+1))*dt;  %BUT HOW CAN I KNOW theta(i+1) here???
        theta_2(i+1) = theta_2(i) + theta_dot_2(i+1)*dt;    
    end
end

if(method ==3 || method ==0)
    %Euler-Cromer
    theta_3 = zeros(num_steps,1);
    theta_dot_3 = zeros(num_steps,1);
     for i = 1:num_steps
        t(i+1) = t(i) + dt;
        theta_dot_3(i+1) = theta_dot_3(i) -sin(theta_3(i))*dt;
        theta_3(i+1) = theta_3(i) + theta_dot_3(i+1)*dt;  
     end
end

%Analytic solution
%from IC's: recall that index 1 corresponds to t=0
B = theta(1);
A = theta_dot(1);

for i = 1:num_steps
    theta_analytic = A*sin(t(i)) + B*cos(t(i));
    theta_dot_analytic = A*cos(t(i)) - B*sin(t(i));
end

%Make plots

%theta and theta_dot vs. t
figure 
 hold on
 if(method ==1 || method ==0)
    plot(t, theta_1)
    plot(t, theta_dot_1)
 elseif(method ==2 || method ==0)
     plot(t, theta_2)
     plot(t, theta_dot_2)
 elseif(method ==3 || method ==0)
     plot(t, theta_3)
     plot(t, theta_dot_3)
 end
 plot(t, theta_analytic)
 plot(t, theta_dot_analytic)
%  title(['Va =', str, 'V'],'interpreter','latex','FontSize',16);
 xlabel('Time','interpreter','latex','FontSize',14);
 ylabel({'Angle'},'interpreter','latex','FontSize',14);





%phase space
