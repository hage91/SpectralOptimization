function y = stepsize(gloquo,gmax,phi_min)

%syms y(x);

h=gmax/3;
x1=1; x2=1+(gmax-1)/3; x3=1+2*(gmax-1)/3;
a=(h-phi_min)/(x3-x2); %%geradenanstieg
b=phi_min-a*x2;

  y = phi_min;
  if gloquo > x2 && gloquo < x3
    y = a*gloquo+b;
  elseif gloquo > x3 
    y = h;
  end
end

    
      

