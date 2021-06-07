function [band,peak] = get_bandwidth(E,A,b,c)

f=linspace(1*10^7,2*10^8,2*10^3);
w=linspace(0,0,2*10^3);


f_im=1i.*f;

%c=zeros(1,12); c(6)=1;
%b=zeros(12,1); b(12)=1;
%c=zeros(1,12); c(19)=1;  %%%entspricht Knoten 26 im Bild
%b=zeros(12,1); b(59)=1;  %%%entspricht Knoten 5 im Bild; das ist irgendwie komisch und passt nicht zu meiner nummerierung


%w=20*log(abs((4*(f_im+d))./(f_im.^2+d*f_im+c)));


for i=1:length(f)
w(i)=20*log10(abs(c*inv(f_im(i)*E-A)*b)); %% mit oder ohne Vorfaktor 20????
end
w=w-w(1);


        for i=2:length(w)
                 if w(i-1)>-3%1/sqrt(2)
                  if w(i)<-3%1/sqrt(2)
                 band=f(i);
                  end
                end
        end

            peak=max(w);
end

