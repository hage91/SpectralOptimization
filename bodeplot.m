function [] = bodeplot(E,A,b,c,Col_val,type)

%[p,k]=get_bandwidth(E,A,b,c)

f=linspace(10^5,10^9,10^4);
w=linspace(0,0,10^4);


f_im=1i.*f;


for i=1:length(f)
w(i)=20*log10(abs(c*inv(f_im(i)*E-A)*b)); %%sometimes prefactor 20
end
w=w-w(1);


semilogx(f,w,type,'Color',Col_val,'LineWidth',3);
ylim([-10 20])




end

