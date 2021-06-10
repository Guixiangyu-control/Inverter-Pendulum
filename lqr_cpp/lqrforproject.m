A=[0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   44.7761 -36.2687 2.68657 0 0 0;
   -73.8806 91.3433 -13.4328 0 0 0;
   -26.8657 -32.2388 22.3881 0 0 0];

B=[0 0 0;
   0 0 0;
   0 0 0;
   44.7761 -36.2687 2.68657;
   -73.8806 91.3433 -13.4328;
   -26.8657 -32.2388 22.3881];

C=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0];

D=zeros(3);

Q=diag([1,1,1,1,1,1]);
R=diag([1,1,1]);

[K,s,e]=lqr(A,B,Q,R);

S=inv(C*inv(B*K-A)*B);

t=0:0.01:6;

u=[pi*ones(size(t))/2;pi*ones(size(t))/2;pi*ones(size(t))/2];
x0=[0 0 0 0 0 0];
sys=ss(A-B*K,B*S,C,D);

[Y,T,X]=lsim(sys,u,t,x0);

subplot(3,2,1)
plot(t,X(1:end,1))
ylabel('\theta_1'); 
xlabel('t'); 

subplot(3,2,2)
plot(t,X(1:end,2))
ylabel('\theta_2'); 
xlabel('t'); 

subplot(3,2,3)
plot(t,X(1:end,3))
ylabel('\theta_3'); 
xlabel('t'); 

subplot(3,2,4)
plot(t,X(1:end,4))
ylabel('a'); 
xlabel('t'); 

subplot(3,2,5)
plot(t,X(1:end,5))
ylabel('b'); 
xlabel('t'); 

subplot(3,2,6)
plot(t,X(1:end,6))
ylabel('c'); 
xlabel('t'); 