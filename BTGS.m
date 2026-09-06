function output= BTGS(P_Name,x_0,type,options)

x_0=x_0(:);

n=size(x_0,1);

x_0=x_0+1e-6*rand(n,1);

disp('Checking input parameters...')

if nargin<3
    error('Error in input parameters: "P_Name", "type" and "starting point" are required input parameters')
end

if nargin==3
    options=[];
end

if type~=0 && type~=1 && type~=2
    error('Error in input parameters: type must be chosen from the set {0,1,2}')
end

if isfield(options, 'm')
   if options.m-floor(options.m)>0 || options.m < 1
      error('Error in input parameters: options.m (size of the sample) must be a positive integer') 
   end
end

if ~isfield(options, 'm')

    if n<=25
        options.m=2*n;
    elseif (25<n) && (n<=50)
        options.m=n;
    elseif (50<n) && (n<=100)
       options.m=floor(n/4);
    elseif (100<n) && (n<=1000)
        options.m=10;
    elseif (1000<n) && (n<=2000)
        options.m=5;
    else
        options.m=2;
    end
end

if ~isfield(options, 'max_iter')
    options.max_iter=50000;
elseif options.max_iter - floor(options.max_iter)>0 || options.max_iter < 1
    error('Error in input parameters: options.max_iter (maximum number of outer iterations) must be a positive integer')
end

if ~isfield(options, 'grad_eval')
    options.grad_eval=2;
elseif options.grad_eval~=0 && options.grad_eval~=1 && options.grad_eval~=2 
    
        error("Error in input parameters: options.grad_eval must be chosen from the set {0, 1, 2}.")
end


if ~isfield(options, 'f_tol')
    options.f_tol=1e-6;
elseif  options.f_tol < 0 || options.f_tol > 0.1
    error('Error in input parameters: options.f_tol must be non-negative and less than 0.1')
end

if ~isfield(options, 'x_tol')
    options.x_tol=1e-6;
elseif  options.x_tol < 0 || options.x_tol > 0.1
    error('Error in input parameters: options.x_tol must be non-negative and less than 0.1')
end


if ~isfield(options, 's_tol')
    options.s_tol=1e-8;
elseif  options.s_tol < 0
    error('Error in input parameters: options.s_tol must be non-negative')
end

if ~isfield(options, 'f_target')
    options.f_target=-inf;
end 

if ~isfield(options, 'max_time')
    options.max_time=inf;
end



if options.grad_eval==2

    if nargout(P_Name) < 2
        error('The user-supplied function must provide both f and g since options.grad_eval=2')
    end
    [f_0,g_0]=P_Name(x_0); g_0=g_0(:);

else
    f_0=P_Name(x_0);

end

if ~isreal(f_0) || size(f_0,1)*size(f_0,2)~=1
    error('Error in input parameters: The user-supplied function must return a real value at the inital point x_0')
end
if isnan(f_0) || abs(f_0)==inf 
    error('Error in input parameters: The objective function is not defined at the inital point x_0')
end

if options.grad_eval==2
if ~isreal(g_0) || size(g_0,1)~=n
    error('Error in input parameters: The user-supplied function must return a real vector of size n as a subgradient of f at the inital point x_0 ')
end

if sum(isnan(g_0))~=0 || ~isempty(find(abs(g_0)==inf))
    error('Error in input parameters: an arbitrary subgradient of the objective function is not defined at inital point x_0')
end
end

s=functions(P_Name);
name=s.function;

disp('Input parameters are valid')
fprintf('minimizing %s with dimension n=%d... \n',name,n);

if type==2
    if options.grad_eval==2
         output=Bundle_GS_PL(P_Name,x_0,options);
    else
         output=Bundle_GS_PL_Approx(P_Name,x_0,options); 
    end
else
    if options.grad_eval==2
         output=Bundle_GS(P_Name,x_0,type,options);
    else
         output=Bundle_GS_Approx(P_Name,x_0,type,options);
    end
end

disp("Preparing a report...")




Table_of_Results=table(output.f_end,output.iters,output.g_eval,output.f_eval,output.cpu_time,output.exit_flag,output.optimality_certificate,'VariableNames',...
{'f_end','Iters','Subgradient_evaluations', 'Function_evaluations','Elapsed_time','Exit_flag','Optimality_certificate'})

end



    






