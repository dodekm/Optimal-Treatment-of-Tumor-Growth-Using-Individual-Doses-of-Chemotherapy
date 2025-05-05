%% Steady states

model = nominal_parameters_Song();
[Nr,Lr,Tr] = get_ss_Song(model,0)
figure;
plot3(Nr,Lr,Tr,'or','LineWidth',2)
hold on
grid3(Nr, Lr, Tr);


for i=1:length(Tr)
    text(Nr(i),Lr(i),Tr(i)+1e8,"  N_0="+num2str(Nr(i),'%.2e'));
    text(Nr(i),Lr(i),Tr(i)+2e8,"  L_0="+num2str(Lr(i),'%.2e'));
    text(Nr(i),Lr(i),Tr(i)+3e8,"  T_0="+num2str(Tr(i),'%.2e') );
end

grid on

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'ZScale', 'log');

xlabel('N');
ylabel('L');
zlabel('T');

[Nr,Lr,Tr] = find_coexisting_stable_state_Song(model,Nr,Lr,Tr,0)

x0=[Nr,Lr,Tr];

%% Finding matrix P for Lyapunov stability

A=get_A_Song([x0';0],model);
A=A(1:3,1:3);
P = lyap(A',diag([1e2;1;1e-2]));
P=P/min(eig(P));
P=P*0.3e-9;


%%  Verifying the Lyapunov condition 
% (maximizing the time derivative subject to inequality constraint)


J_fun=@(x)(-get_V_dot_Song_((x+x0)',model,x0',P));

initialX=1.5*chol(inv(P))*randn(3,10000);
initialX=initialX(:,0.5*diag(initialX'*P*initialX)<=1)

options = optimoptions('ga','FunctionTolerance',1e-50,'MaxGenerations',20,'Display','iter','PopulationSize',size(initialX,2),'InitialPopulationMatrix',initialX');
x =x0 + ga(J_fun,3,[],[],[],[],[],[],@(x)Lyap_con((x+x0)',x0',P),options);

-J_fun(x-x0)


%% Plotting the ellipsoid stability region
figure

ellipsoid_plot(P/2,x0);
hold on
plot3(Nr,Lr,Tr,'or','LineWidth',2);
legend('V(N,L,T) \leq 1','N_0,L_0,T_0');

xlabel('N [cells]');
ylabel('L [cells]');
zlabel('T [cells]');



%% Numerical mapping of stability region
close all

x_0_stability_region=[];
x_0_unstability_region=[];

N0_=linspace(0,1e6,100);
L0_=linspace(0,1e7,100);
T0_=linspace(0,1e9,100);

tf=500;
Ts=1/10;

for i=1:length(N0_)
   
    for j=1:length(L0_)
        for k=1:length(T0_)

            N0=N0_(i);
            L0=L0_(j);
            T0=T0_(k);

            [N,L,T] = sim_Song(tf,Ts,N0,L0,T0,0,zeros(1,tf/Ts+1),model);

            if Lyap_con([N(end);L(end);T(end)],[Nr;Lr;Tr],P)<0
                x_0_stability_region=[x_0_stability_region; [N0,L0,T0]];
            else
                x_0_unstability_region=[x_0_unstability_region; [N0,L0,T0]];
            end
        end
    end
end

save('stability_region','x_0_stability_region','x_0_unstability_region')
%% Rendering the mapped stability region
figure
%shp = alphaShape(x_0_stability_region(:,1),x_0_stability_region(:,2),x_0_stability_region(:,3),3e20);
shp = alphaShape(x_0_stability_region(:,1),x_0_stability_region(:,2),x_0_stability_region(:,3),3e9);
plot(shp,'FaceAlpha',0.5,'EdgeAlpha',0.0,'FaceColor',[0,0,0.9]);
hold on

axis square

% shading interp;               % Smooth color interpolation
light('Position', [0 -1e7 5e8], 'Style', 'infinite');
lighting gouraud;             % Smooth lighting effect
material shiny;               % Shiny surface appearance
light('Position', [1e7 1e7 5e8], 'Style', 'local');
lighting phong;

grid on

hold on
plot3(Nr,Lr,Tr,'or','LineWidth',2);
legend('stability region','N_0,L_0,T_0')
xlabel('N');
ylabel('L');
zlabel('T');

%set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
%set(gca, 'ZScale', 'log');

%% Chemotherapy optimization

TD=14;
nD=5;

w_N=1e4;
w_L=1e2;
w_T=1e0;
w_u=1e16;
w_s=1e20;

tf=100;

Ts=1/100;


N0=0.5*Nr;
L0=0.5*Lr;
T0=1*1e8;

% N0=0.1*Nr;
% L0=0.1*Lr;
% T0=2*1e8;


% N0=Nr*1e-2;
% L0=Lr*1e-2;
% T0=5*1e8;

J_fun=@(D)J_Song_terminal_stability(D,TD,tf,Ts,w_N,w_L,w_T,w_u,w_s,P,Nr,Lr,Tr,N0,L0,T0,model);
options = optimoptions('ga','Display','iter','MaxGenerations',20);
D = ga(J_fun,nD,[],[],[],[],zeros(nD,1),5*ones(nD,1),[],options);
J_fun(D)

%% Simulating and plotting the treatment response

[N,L,T,u,v,t] = sim_treatment_Song(D,TD,tf,Ts,N0,L0,T0,model);

figure
subplot(5,1,1);
plot(t,N,'-k','LineWidth',1.5);
hold on
plot(t,Nr*ones(size(t)),'--g')
set(gca, 'YScale', 'log');
xlabel('t [d]')
ylabel('N [cells]')
legend('N(t)','N_0');

subplot(5,1,2);
plot(t,L,'-k','LineWidth',1.5)
hold on
plot(t,Lr*ones(size(t)),'--g')
set(gca, 'YScale', 'log');
xlabel('t [d]')
ylabel('L [cells]')
legend('L(t)','L_0');

subplot(5,1,3);
plot(t,T,'-k','LineWidth',1.5)
set(gca, 'YScale', 'log');
hold on
plot(t,Tr*ones(size(t)),'--g')
xlabel('t [d]')
ylabel('T [cells]')
legend('T(t)','T_0');

subplot(5,1,4);
plot(t,u,'-k','LineWidth',1.5)
xlabel('t [d]')
ylabel('u [IU]')
legend('u(t)')

subplot(5,1,5);
stairs(t,v,'-k','LineWidth',1.5)
xlabel('t [d]')
ylabel('v [IU/d]')
legend('v(t)')

set(gcf,'position',[0,0,450,800]);
%% Simulating and plotting the treatment state trajectory

[N,L,T] = sim_treatment_Song(D,TD,tf+100,Ts,N0,L0,T0,model);

figure

plot3(N(end-1.405e4:end),L(end-1.405e4:end),T(end-1.405e4:end),'-r','LineWidth',1);
hold on
grid on
%plot3(N0,L0,T0,'og','LineWidth',1);
plot3(N(end),L(end),T(end),'og','LineWidth',1);
plot3(Nr,Lr,Tr,'ob','LineWidth',1);

hold on
ellipsoid_plot(P/2,x0);

%legend('N(t),L(t),T(t)','N(0),L(0),T(0)','N_0,L_0,T_0')
legend('N(t),L(t),T(t)','N(t_f),L(t_f),T(t_f)','N_0,L_0,T_0','V(N,L,T) \leq 1')

xlabel('N [cells]');
ylabel('L [cells]');
zlabel('T [cells]');

%set(gca, 'ZScale', 'log');





