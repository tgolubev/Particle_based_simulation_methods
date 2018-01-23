%Simple harmonic oscillator integrator comparison

%Solves for simple harmonic oscillator with hamiltonian
% H = p^2/2m + k/2(x-x_eq)^2

%Timofey Golubev

%time step and range
t_max = 100;
m =900;
k=100;
w = sqrt(k/m);
dt = 0.1*w;       %time step is changing the period of the solution!!--> THIS SHOULD NOT HAPPEN!
num_steps = t_max/dt + 1;  %+1 to include 0th step


%initial conditions
x_eq = 1; %equil. lenght in meters
x_0 = 1.01;  %initial position
p_0 = 0; %initial momentum


%initialize arrays
t = zeros(num_steps,1);
t(1) = 0;  %note: b/c matlab starts indices from 1, need to use this


%forward Euler
x_1 = zeros(num_steps,1);
p_1 = zeros(num_steps,1);
x_1(1) = x_0;
p_1(1) = p_0;
for i = 1:num_steps-1
    t(i+1) = t(i) + dt;
    x_1(i+1) = x_1(i) + (p_1(i)/m)*dt;
    p_1(i+1) = p_1(i) -k*(x_1(i) - x_eq)*dt;
end

%backward Euler
x_2 = zeros(num_steps,1);
p_2 = zeros(num_steps,1);
x_2(1) = x_0;
p_2(1) = p_0;
for i = 1:num_steps-1
    p_2(i+1) = p_2(i) -k*((x_2(i) + (p_2(i)/m)*dt) - x_eq)*dt;  %estimate x_2(i+1), using x_dot*dt + x(i), and x_dot = v = p/m
    x_2(i+1) = x_2(i) + (p_2(i+1)/m)*dt; 
end

%Euler-Cromer
x_3 = zeros(num_steps,1);
p_3 = zeros(num_steps,1);
x_3(1) = x_0;
p_3(1) = p_0;
for i = 1:num_steps-1
    p_3(i+1) = p_3(i) -k*(x_3(i) - x_eq)*dt;
    x_3(i+1) = x_3(i) + (p_3(i+1)/m)*dt;
end

%analytic solution: x(t) = Asin(wt) + Bcos(wt) + x_eq, where w = sqrt(k/m)
%from IC's we find the below expressions for the coefficients
B = x_0 - x_eq;
A = p_0/(m*w);
x_analytic = zeros(num_steps,1);
p_analytic = zeros(num_steps,1);
for i = 1:num_steps
    x_analytic(i) = A*sin(w*t(i)) + B*cos(w*t(i)) + x_eq;
    p_analytic(i) = m*(A*w*cos(w*t(i)) - B*w*sin(w*t(i)));  %p = x_dot*m
end
    
%period should be 2pi/w--> right now analytic solution seems right, and all
%others seem wrong!--> periods are too short!

%calculate energies
%KE
KE_1 = p_1.^2/(2*m);
KE_2 = p_2.^2/(2*m);
KE_3 = p_3.^2/(2*m);
KE_analytic = p_analytic.^2/(2*m);

%PE
U_1 = 0.5*k*(x_1-x_eq).^2;
U_2 = 0.5*k*(x_2-x_eq).^2;
U_3 = 0.5*k*(x_3-x_eq).^2;
U_analytic = 0.5*k*(x_analytic-x_eq).^2;

%plot results

%plot KE and total energy
figure
hold on
 p1 = plot(t,KE_1,'b','LineWidth',1);
 p2 = plot(t,KE_2,'r','LineWidth',1);
 p3 = plot(t,KE_3,'Color',[0.93 0.69 0.125]','LineWidth',1);
 p4 = plot(t,KE_analytic,'k','LineWidth',1);

 plot(t,KE_1+U_1,'b','LineWidth',1);
plot(t,KE_2+U_2,'r','LineWidth',1);
plot(t,KE_3+U_3,'Color',[0.93 0.69 0.125]','LineWidth',1);
plot(t,KE_analytic+U_analytic,'k','LineWidth',1);

set(gca, 'FontSize', 20)
 xlabel('Time(s)','interpreter','latex','FontSize',22);
 ylabel({'Energy(J)'},'interpreter','latex','FontSize',22);
 
 hold off
legend([p1 p2 p3 p4],'Forward Euler','Backward Euler', 'Euler-Cromer', 'Analytic','Location','northeast');
 
 %plot U and total energy
 figure
 hold on
p1 = plot(t,U_1,'b','LineWidth',1);
p2 = plot(t,U_2,'r','LineWidth',1);
p3 = plot(t,U_3,'Color',[0.93 0.69 0.125]','LineWidth',1);
p4 = plot(t,U_analytic,'k','LineWidth',1);

plot(t,KE_1+U_1,'b','LineWidth',1);
plot(t,KE_2+U_2,'r','LineWidth',1);
plot(t,KE_3+U_3,'Color',[0.93 0.69 0.125]','LineWidth',1);
plot(t,KE_analytic+U_analytic,'k','LineWidth',1);
set(gca, 'FontSize', 20)
 xlabel('Time(s)','interpreter','latex','FontSize',22);
 ylabel({'Energy(J)'},'interpreter','latex','FontSize',22);
 
 hold off
legend([p1 p2 p3 p4],'Forward Euler','Backward Euler', 'Euler-Cromer', 'Analytic','Location','northeast');

%plot positions
figure 
hold on
p1 = plot(t,x_1,'b','LineWidth',1);
p2 = plot(t,x_2,'r','LineWidth',1);
p3 = plot(t,x_3,'Color',[0.93 0.69 0.125]','LineWidth',1);
p4 = plot(t,x_analytic,'k','LineWidth',1);
set(gca, 'FontSize', 20)
 xlabel('Time(s)','interpreter','latex','FontSize',22);
 ylabel({'Position(m)'},'interpreter','latex','FontSize',22);
 hold off
legend([p1 p2 p3 p4],'Forward Euler','Backward Euler', 'Euler-Cromer', 'Analytic','Location','northeast')



 
 %plot KE U and total energy for analytic case
%  figure 
%  hold on
% plot(t,KE_analytic,'b','LineWidth',1);
% plot(t,U_analytic,'r','LineWidth',1);
%  plot(t,KE_analytic+U_analytic,'k','LineWidth',1);
%  set(gca, 'FontSize', 20)
%  xlabel('Time(s)','interpreter','latex','FontSize',22);
%  ylabel({'Energy(J)'},'interpreter','latex','FontSize',22);
%  

