function lambda = QP_SOLVER(GOSM,alpha,tau,eps,x)

m=size(GOSM,2);

H=GOSM'*GOSM; f=eps^(-tau)*alpha; Aeq=ones(1,m); beq=1; lb=zeros(1,m);

H=(H+H')/2;

if isempty(H) || ~isreal(H)
    error('An unexpected error occured while using quadprog solver. Make sure that function values and subgradient vectors remain real. Trying another initial point is also recommended');
end

qp_options = optimset('Display','off','TolX', 1e-12, 'TolFun', 1e-12);


x_0=ones(size(H,2),1)/size(H,2);
lambda=quadprog(H, f, [], [], Aeq, beq, lb, [],x_0,qp_options);




if isempty(lambda) || ~isreal(lambda) || sum(lambda)>1.1 || sum(lambda)<.9 || min(lambda)<-0.1
    error('An unexpected error occured while using quadprog solver. Make sure that function values and subgradient vectors remain real. Trying another initial point is also recommended');
end

end

