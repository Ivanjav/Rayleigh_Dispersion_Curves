function F_R = dispersion(vs,vp,rho,h,fq,k)
n=length(vs);
K=zeros(2*n);
w=2*pi*fq;
c=w/k;
r=sqrt(1-(c^2)./(vp.^2));
s=sqrt(1-(c^2)./(vs.^2));
D=2*(1-cosh(k*r(1:end-1).*h).*cosh(k*s(1:end-1).*h))...
    +(1./(r(1:end-1).*s(1:end-1))+(r(1:end-1).*s(1:end-1))).*...
    sinh(k*r(1:end-1).*h).*sinh(k*s(1:end-1).*h);
for j=1:n-1
    k11=((k*rho(j)*c^2)/D(j))*((s(j)^(-1))*...
         cosh(k*r(j)*h(j))*sinh(k*s(j)*h(j))-r(j)*...
         sinh(k*r(j)*h(j))*cosh(k*s(j)*h(j)));
    k12=((k*rho(j)*c^2)/D(j))*(cosh(k*r(j)*h(j))*...
        cosh(k*s(j)*h(j))-r(j)*s(j)*sinh(k*r(j)*h(j))...
        *sinh(k*s(j)*h(j))-1)-k*rho(j)*(vs(j)^2)*(1+s(j)^2);
    k13=((k*rho(j)*c^2)/D(j))*(r(j)*sinh(k*r(j)*h(j))-...
        (s(j)^(-1))*sinh(k*s(j)*h(j)));
    k14=((k*rho(j)*c^2)/D(j))*(-cosh(k*r(j)*h(j))+...
        cosh(k*s(j)*h(j)));
    k21=k12;
    k22=((k*rho(j)*c^2)/D(j))*((r(j)^(-1))*...
         sinh(k*r(j)*h(j))*cosh(k*s(j)*h(j))-s(j)*...
         cosh(k*r(j)*h(j))*sinh(k*s(j)*h(j)));
     k23=-k14;
     k24=((k*rho(j)*c^2)/D(j))*(-(r(j)^(-1))*...
         sinh(k*r(j)*h(j))+s(j)*sinh(k*s(j)*h(j)));
     k31=k13;
     k32=k23;
     k33=k11;
     k34=-k12;
     k41=k14;
     k42=k24;
     k43=-k21;
     k44=k22;
     
 K11=[k11 k12;k21 k22];
 K12=[k13 k14;k23 k24];
 K21=[k31 k32;k41 k42];
 K22=[k33 k34;k43 k44];

      if j==1
          K(1:2,1:4)=[K11 K12];
      else
          rows=(3:4)+2*(j-2);
          columns=(1:6)+2*(j-2);
          K(rows,columns)=[K21a K22a+K11 K12];
      end
K21a=K21;
K22a=K22;     
end    
k11e=(k*rho(n)*(vs(n)^2))*((r(n)*(1-s(n)^2))/(1-r(n)*s(n)));
k12e=(k*rho(n)*(vs(n)^2))*((1-s(n)^2)/(1-r(n)*s(n)))-...
    2*k*rho(n)*(vs(n)^2);
k21e=k12e;
k22e=(k*rho(n)*(vs(n)^2))*((s(n)*(1-s(n)^2))/(1-r(n)*s(n)));

Ke=[k11e k12e;k21e k22e];
K(end-1:end,end-3:end)=[K21 K22+Ke];

F_R=real(det(K));
end