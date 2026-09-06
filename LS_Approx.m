function [s,F_s,G_s,I,f_eval,g_eval] = LS_Approx(gamma,eps,x,d,w_opt,P_Name,ff,f_eval,g_eval,options,i,k_out)
% Some initializations
t_bar=0.01; 

Beta_R=.1; Beta_A=sqrt(eps)*Beta_R/4; Beta_T=Beta_A; Beta_L=1e-6*sqrt(eps)*Beta_R/5;

t=1; t_U=1; t_L=0; eta=0.25; s=x+t*d;  

I=-1;


h_0=1e-8;
h=h_0;



% k is a counter
k=0; 

while k<50 && t<1e+2 
    
    F_s=P_Name(s);
    G_s=Grad_Eval(P_Name,s,h);

     if any( isinf(G_s)) || any(isnan(G_s)) || isinf(F_s) || isnan(F_s) 

        break;

    end

    f_eval=f_eval+1; g_eval=g_eval+1;
    
    alpha=max([abs(ff-(F_s+dot(G_s,x-s))),gamma*sum((x-s).^2)]);
    
% First conditional block
    if F_s-ff<=-Beta_T*t*w_opt 
        t_L=t;
    else 
        t_U=t;
    end

    
    
% Second conditional block    
    if F_s-ff<=-Beta_L*t*w_opt && (t >= t_bar || alpha>Beta_A*w_opt)
      
        I=1;
        break;
    end
    
% Third conditional block    
    if -alpha+dot(G_s,d)> -Beta_R*eps*w_opt
        
        I=0;
        break;
    end
    
    if t_L==0
        
        F_U=P_Name(x+t_U*d);
        f_eval=f_eval+1;
        
        t=max(eta*t_U, (0.5*w_opt*t_U^2)/(t_U*w_opt+ff-F_U));
        
    else
        t=(t_U+t_L)/2;
        
    end
    s=x+t*d;

    k=k+1;
    if options.grad_eval==1 && h>1e-12
    h=h_0/sqrt(k+i+k_out);
    end

end



end

function G=Grad_Eval(P_Name,x,h)
n=size(x,1);
G=zeros(n,1);
F=P_Name(x);
for i=1:n
    e=zeros(n,1); 
    e(i)=1;
    G(i)=(P_Name(x+h*e)-F)/h;
end
G=G(:);
end

