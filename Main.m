m=50;n=5000;%m=25;n=2500
Targeting=0.2;%0.2;%0.8
[s,d,p,w,link,theta]=GenerateGDA(m,n,Targeting);
% Set up lambda, and run solver 
[u_v,v_v,Objective_v] = Gurobi_GDA(s,d,p,w,link,theta);
% if Objective_v<=1e-8
%    Objective_v=0;
% end
%index=find(u_v<=1e-8);
%u0=u_v.*1.2;
%u0(index)=0;
u0=zeros(m,n);
iteration = 10000;%4000;
%delta=0.5;
t00=cputime;
[x_VAPP,y_VAPP,h_VAPP] = VAPP_GDA(s,d,p,w,link,theta,iteration,u0,Objective_v);
t11=cputime-t00;

figure(1); 
semilogy(1:iteration,h_VAPP.obj,'b');
legend('VAPP-GDA')
xlabel('iteration'); ylabel('suboptimality');
figure(2);
semilogy(1:iteration,h_VAPP.constraint,'b');
legend('VAPP-GDA');
xlabel('iteration'); ylabel('feasibility');
figure(3);
semilogy(1:iteration,h_VAPP.relation,'b');
legend('VAPP-GDA');
xlabel('iteration'); ylabel('relation');
figure(4); 
semilogy(1:iteration,h_VAPP.plus,'b');
%axis([0 iteration 1e-8 1e2]);
legend('VAPP-GDA')
xlabel('iteration'); ylabel('|F(U^k,v^k)-F(U^*,v^*)|+||max\{0,\Theta(U^k,v^k)\}||');