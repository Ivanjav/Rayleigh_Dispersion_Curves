function cR = dispersion_modes(vs,vp,rho,h,f,c0)

for l=1:length(f)
fq=f(l);
wq=2*pi*fq;
k0=wq./c0;
for i=1:length(k0)
Phi(i)=dispersion(vs,vp,rho,h,fq,k0(i));
end

dPhi=diff(sign(Phi));
dPhi2=[0 diff(diff(abs(Phi)))];
dPhi2=interp1(1:length(k0)-1,dPhi2,(1:length(k0)-1)+1/3);
kk=k0((abs(dPhi)==2)&(dPhi2>0));
cR(1:length(kk),l)=(wq./kk)*1e3;
end

cR(cR==0)=NaN;

