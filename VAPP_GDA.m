function [x,y,h]=VAPP_GDA(s,d,penal,w,link,theta,iteration,u_v,fstar)
[m,n]=size(link);
x=u_v;
y=zeros(m,1);
p1=zeros(m,1);
p2=zeros(n,1);
gamma=1;
epsilon=0.86;%0.00016;%50-0.2;%0.000044;%25-0.8;
for i=1:iteration
    xk=x;
    yk=y;
    p1k=p1;
    p2k=p2;
    q1k=max(0,p1k+gamma*(d-yk-(link.*xk)*s));
    q2k=max(0,p2k+gamma*((link.*xk)'*ones(m,1)-ones(n,1)));
    y=max(yk-epsilon*(penal*ones(m,1)-q1k),0);
    x=max(xk-epsilon*(((w./theta)*s').*link.*(xk-theta*ones(1,n))-q1k*ones(1,n).*link.*(ones(m,1)*s')+ones(m,1)*q2k'.*link),0);
    %x=max((-delta/alpha),min((delta/alpha),x));
    p1=max(0,p1k+gamma*(d-y-(link.*x)*s));
    p2=max(0,p2k+gamma*((link.*x)'*ones(m,1)-ones(n,1)));
    %save informations;
    h.obj(i)=abs(norm(((w./(2*theta))*s').^(1/2).*link.*(xk-theta*ones(1,n)),'fro')^2+penal*sum(yk)-fstar);
    h.constraint(i)=norm(max(0,d-yk-(link.*xk)*s),2)+norm(max(0,xk'*ones(m,1)-ones(n,1)),2);
    h.relation(i)=norm(x-xk,2)/max(norm(x,2),1);
    h.plus(i)=h.obj(i)+h.constraint(i);
end
%h.obj(1)=abs(norm(((w./(2*theta))*s').^(1/2).*link.*(zeros(m,n)-theta*ones(1,n)),'fro')^2+0-fstar);
%h.plus(1)=h.obj(1)+h.constraint(1);