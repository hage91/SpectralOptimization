  %  isens=isens+1;
   % sens=zeros(sz,1);        
%for b=1:sz
 %   if abs(Ute(b,b))>domsens
  %      sens(b)=b;
   % end
    
%end


%%%Sensitiven Faktoren mit dazu nehmen
%for i1=1:num
 %   for l=1:sz
    
  %      if l~=iplus(1:num)
            
   %         if l~=iminus(1:num)
          
    %        fac(i1)=fac(i1)*V1(l,l)/(V1(l,l)+V1(l,l)*step*uinf(Iinv(l))*vinf(Iinv(l)))*(V1(iplus(i1),iplus(i1))-V1(l,l)-step*V1(l,l)*uinf(Iinv(l))*vinf(Iinv(l)))/(V1(iplus(i1),iplus(i1))-V1(l,l));
            
     %       end
            
      %  end
        
    %end
%end


%fac=1;
%for i1=1:num   %%%%V1 überall ersetzen durch EW, weil wir die Indexmenge vereinheitlichen müssen
%for l=1:sz
%if abs(V1(l,l))<10^(16) %%%EW bei oo ausschließen
 %   if sens(Iinv(l))>0  %%%für alle sensitiven Eigenwerte nehme ich die Faktoren in das Produkt auf
 %   if l==iplus(i1)
 %   fac(i1)=fac(i1)/m(i1)*V1(iminus(i1),iminus(i1))*(V1(l,l)-m(i1))*(V1(l,l)-conj(m(i1)))/((V1(l,l)-V1(iminus(i1),iminus(i1)))*conj(m(i1)));
    %%%Hier wird der zweite sensitive Faktor auch berücksichtigt (oder allgemeiner alle sensitiven EW)
 %   end
  %      if l~=iplus(1:num)
            
   %         if l~=iminus(1:num)
          
    %        fac(i1)=fac(i1)*V1(l,l)/(V1(l,l)+V1(l,l)*step*uinf(Iinv(l))*vinf(Iinv(l)))*(V1(iplus(i1),iplus(i1))-V1(l,l)-step*V1(l,l)*uinf(Iinv(l))*vinf(Iinv(l)))/(V1(iplus(i1),iplus(i1))-V1(l,l));
            
     %       end
            
      %  end
      
   % end
    
  %  for i2=1:num    %%%setze die nicht sensitiven, die ich aber vorgebe
  %  fac(i1)=fac(i1)/m(i2)*V1(iminus(i2),iminus(i2))*V1(iplus(i2),iplus(i2))*(V1(iplus(i1),iplus(i1))-m(i2))*(V1(iplus(i1),iplus(i1))-conj(m(i2)))/((V1(iplus(i1),iplus(i1))-V1(iplus(i2),iplus(i2)))*(V1(iplus(i1),iplus(i1))-V1(iminus(i2),iminus(i2)))*conj(m(i2)));    
    
  %  end     
    
%end
%end
%end

