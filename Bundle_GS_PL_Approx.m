function output = Bundle_GS_PL_Approx(P_Name, x_0, options)

tic

%some initializations

x=x_0;

n=size(x,1);

f_rec=inf*ones(1,options.max_iter);

w_rec=inf*ones(1,options.max_iter);

i_w=1;

N_F_LS=1;

Max_of_Failures=25;

m=options.m;

eps_max=4;

eps_min=1e-6;

tau=0.5;

t_bar=0.01;

Delta= 0.5;

kesi=0.9;

eps=1;

w_opt=inf;

F_s=inf;

theta=1;

T=0;

gamma=0.0;

Lambda_Indicator=0;

zeta_bar=1;

safe_guard=-1;

SM_p=[];

GOSM_p=[];

FOSM_p=[];

SM=[];


if options.grad_eval==1
h_0=1e-7;
h=h_0;
else
h_0=1e-8;
h=h_0;
end

% Counters

k=0;

k_i=0;

f_eval=0;

g_eval=0;

%Outer loop

f_rec(k+1)=P_Name(x);

while  (w_opt > options.s_tol || safe_guard~=4) && N_F_LS < Max_of_Failures && F_s >= options.f_target   && k < options.max_iter && toc < options.max_time && T~=1
    
    
    if size(SM,2)>=100

        SM_p=[];

        FOSM_p=[];

        GOSM_p=[];

    end
    
% Sample m point from B(x,eps)

    RGP=2*rand(n,m)-1;

    SM=[x,kron(x,ones(1,m))+eps*RGP];

 %Appending active points from previous iteration

 SM=[SM,SM_p]; 
   
% Function and Gradient evaluations at sampled points

 for j=1:m+1

     F=P_Name(SM(:,j));

     FOSM(j)=F;

     G=Grad_Eval(P_Name,SM(:,j),h);

     GOSM(:,j)=G;
 end

 g_eval=g_eval+(m+1);

 f_eval=f_eval+(m+1);
 
 FOSM=[FOSM,FOSM_p];

 GOSM=[GOSM,GOSM_p];
   
 ff=FOSM(1);


% Gradient locality measures(2-row matrix)

alpha=max([abs(ff-(FOSM+sum(GOSM.*(x-SM),1)));gamma*(sum((x-SM).^2,1))]); 
 
% Solving the quadratic subproblem of outer iteration 
try
lambda=QP_SOLVER(GOSM,alpha,tau,eps,x);
catch
    N_F_LS=N_F_LS+1;

    eps=0.5*eps;

    k=k+1;

    f_rec(k+1)=ff;

    continue

end
    
% Updating variables of outer iteration

sigma_t=dot(lambda,sqrt(sum((((x-SM).^2)),1)));

f_t=dot(lambda,FOSM+sum(GOSM.*(x-SM),1));

alpha_t=max(abs(ff-f_t),gamma*sigma_t^2);

g_t=GOSM*lambda;

d=-eps^tau*g_t;

alpha_hat=dot(lambda,alpha);

w_opt=.5*norm(g_t)^2+1/(eps^tau)*alpha_hat;

w_rec(i_w)=w_opt; i_w=i_w+1;

if w_opt< options.s_tol
    
    safe_guard=safe_guard+1;

    k=k+1;

    if options.grad_eval==1
        h_0=.1*h_0;
       else


        h=.1*h;
     end

end

%Applying gradient selection strategy in outer iteration (reducing the index set)

[GOSM,alpha,SM,FOSM]=Adaptive_Grad_Select(lambda,theta,GOSM,alpha,SM,FOSM);
    
%Beginning of inner iterations

        i=0;

        while w_opt > options.s_tol && toc < options.max_time
            
                  k_i=k_i+1;

                  %Line search procedure

                  [s,F_s,G_s,I,f_eval,g_eval]=LS_Approx(gamma,eps,x,d,w_opt,P_Name,ff,f_eval,g_eval,options,i,k);

                  if I==-1

                       eps=0.5*eps;

                       N_F_LS=N_F_LS+1;

                       k=k+1;

                       f_rec(k+1)=ff;
                      
                      break
                      
                  end
                  
                  if I==1
                      
                      if (ff-F_s)/(abs(F_s)+1) < options.f_tol && norm(x-s,inf)/((norm(s,inf))+1)<options.x_tol

                          T=1;
                      end

                      r_2=max(-alpha+(GOSM'*(d))');

                      r_3=max(-alpha);

                      r=(F_s-ff)/(r_2-r_3);

                      if abs(P_Name(x+t_bar*d)-P_Name(x))<= zeta_bar*w_opt && Lambda_Indicator==0

                          Lambda_Indicator=1;

                          f_eval=f_eval+1;

                      end
                      
                      
                      if r>= Delta && 2*eps<=eps_max 

                          x=s;
                          
                          eps=2*eps;

                          k=k+1;

                          f_rec(k+1)=F_s;
                          
                          
                          if Lambda_Indicator==1

                              eps=max(eps,eps_min);

                          end

                          break;

                      else

                          x=s;
                          
                          k=k+1;

                          f_rec(k+1)=F_s;


                          if Lambda_Indicator==1

                              eps=max(eps,eps_min);

                          end

                          break;
                          
                      end
                      
                   end
                  
            if I==0
                
            % Computing the gradient locality measure of the new auxiliary point
               
              alpha_new=max(ff-(F_s+dot(G_s,x-s) ), gamma*norm(x-s)^2);
              
               % Checking a criterion for enriching the model function

               f_eval=f_eval+1;

              if alpha_new <= kesi*alpha_t || abs(P_Name(x+t_bar*d)-P_Name(x))<= zeta_bar*w_opt 
                  
                 % The model function is enriched 
                 
                 GOSM=[GOSM,G_s]; alpha=[alpha,alpha_new]; SM=[SM,s]; FOSM=[FOSM,F_s];
                 
                 % Solving the quadratic subproblem of the inner iteration
                    try
                    lambda=QP_SOLVER([GOSM,g_t],[alpha,alpha_t],tau,eps,x);
                    catch
                                            
                      eps=0.5*eps;
                      
                      N_F_LS=N_F_LS+1;

                      k=k+1;

                      f_rec(k+1)=ff;

                    break

                    end
                     
                 
                 % Updating the variables of the inner iteration

                 sigma_t=dot(lambda,[sqrt(sum((x-SM).^2,1)),sigma_t]);

                 f_t=dot(lambda,[FOSM+sum(GOSM.*(x-SM),1),f_t]  );

                 alpha_t=max(abs(ff-f_t),gamma*sigma_t^2);

                 alpha_hat=dot(lambda,[alpha,alpha_t]);

                 g_t=[GOSM,g_t]*lambda;

                 d=-eps^tau*g_t;

                 w_opt=.5*norm(g_t)^2+1/(eps^tau)*alpha_hat;

                 w_rec(i_w)=w_opt; i_w=i_w+1;

                 if w_opt< options.s_tol
    
                        safe_guard=safe_guard+1;

                        k=k+1;

                        if options.grad_eval==1
                            h_0=.1*h_0;
                        else


                        h=.1*h;
                        end

                 end
                 
                                  
                 %Applying the adaptive gradient selection strategy in inner iteration (reducing the index set)

                 [GOSM,alpha,SM,FOSM]=Adaptive_Grad_Select2(lambda,theta,GOSM,alpha,SM,FOSM);

              else % The sampling radius is reduced

                    eps=.5*eps;

                    k=k+1;

                    f_rec(k+1)=ff;

                    
                  break;
              end
           end
                          
         i=i+1;
         
         

         
         if i==100

             eps=0.5*eps;

             k=k+1;

             f_rec(k+1)=ff;
             
             break;
           
          end
 
         
         
        end
    
      


if options.grad_eval==1 && h>1e-13

h=h_0/sqrt(k);

end

SM_p=SM;

GOSM_p=GOSM;

FOSM_p=FOSM;

GOSM=[];

FOSM=[];

if mod(k,50)==0
    fprintf("The best obtained value of the objective function up to iteration %d is %d \n",k, min(f_rec));
end

end


fprintf("The best obtained value of the objective function is %d \n", min(f_rec));
        
if w_opt <= options.s_tol
    flag=1;
elseif F_s <= options.f_target
    flag=3;
elseif k >= options.max_iter
    flag=4;
elseif toc >= options.max_time 
    flag=5;
elseif T==1
    flag=2;
elseif N_F_LS==Max_of_Failures
    flag=0;
end
    
output.cpu_time=toc;

output.iters=k;

output.f_end=min(f_rec);

output.x_end=x; 

output.f_eval=f_eval;

output.g_eval=g_eval;

output.exit_flag=flag;

output.inner_iters=k_i;

output.f_rec=f_rec(1:k+1);

output.optimality_certificate=min(w_rec);

disp("End of the optimization process")

end


function G=Grad_Eval(P_Name,x,h)

n=size(x,1);

G=zeros(n,1);

ff=P_Name(x);

for i=1:n

    e=zeros(n,1); 

    e(i)=1;

    G(i)=(P_Name(x+h*e)-ff)/h;

end

G=G(:);

end
