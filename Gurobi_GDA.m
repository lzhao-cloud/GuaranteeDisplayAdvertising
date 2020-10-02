function [u_v,v_v,Objective_v] = Gurobi_GDA(s,d,p,w,link,theta)
[m,n]=size(link);
u=sdpvar(m,n,'full');
v=sdpvar(m,1,'full');
Objective=norm(((w./(2*theta))*s').^(1/2).*link.*(u-theta*ones(1,n)),'fro')^2+p*sum(v);
% for j=1:m
%     Objective=Objective+p*v(j,1);
%     for i=1:n
%         if link(j,i)==1
%            Objective=Objective+(w(j,1)/(2*theta(j,1)))*s(i,1)*(u(j,i)-theta(j,1))^2;
%         end
%     end
% end
CX=[];
CX=[CX;d-v-(link.*u)*s<=0];
CX=[CX;(link.*u)'*ones(m,1)-ones(n,1)<=0];
CX=[CX;u>=0];
CX=[CX;v>=0];
options=sdpsettings('verbose',1,'solver','gurobi');
sol=optimize(CX,Objective,options);
Objective_v=value(Objective);
u_v=value(u);
v_v=value(v);
sol.info