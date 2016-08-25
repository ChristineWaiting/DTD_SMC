%IMPVALUE

%Purpose : Compute the implied asset value given equity value. A bisection
%search is used until the upper and lower bounds are within bi_eps of each
%other (relative).
%Then Newton-Raphson iterations are performed until the step taken
%is less than bi_nr * av. If a Newton-Raphson iterations take a step outside of
%the bisection bounds, the bisection iterations resume until convergence is
%achieved.
%
%Note that the bisection iterations are stopped after 53 iterations since
%the precision of a double is 2^-52


%Input   : parval - Nx3 matrix (Row elements: 1 (equity value),
%          2 (debt value), 3 (interest rate) )
%           sig - volatility
%           t - maturity,expressed in daily.
%           av0 - initial input of asset values

%Output  : av     - Nx1 implied asset values

%Input and output: no missing values

function [av, Nd1] = impvalue_annual(eqval, debt, r, sig, t, av_ini)
%fprintf('-')

% convergence parameters
BI_EPS = 0.1;
tol = 1e-6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxiter_newton = 10; % max iterations for Newton's method
maxiter_bisection = ceil(log(tol)/log(0.5)); % max iterations for bisection method

[n Nf] = size(eqval);

if strcmp('gpuArray',class(eqval))
    
    zeroMTX=gpuArray.zeros(n,Nf);
else
    zeroMTX=zeros(n,Nf);
end

if size(sig,1)==1
    sig=repmat(sig,n,1);
end

%precompute fixed values
rsigsq = bsxfun(@plus, r , sig .* sig / 2) * t;

sigsqrtt = sig * sqrt(t);
discfact = exp(-r*t);

%initialize upper and lower bounds
if ~isnan(av_ini) 
    av = av_ini;
else
    d_upper = 0.75 * (eqval + debt);
    d_lower = 0.8 * d_upper;
    
    %increase the upper bound until the model equity value is greater than eqval
    success_u = 0;
    av = d_upper; %initialize asset value
    while (success_u~=1) %set the upper bound
        d = bsxfun( @rdivide,(log(av./debt) + rsigsq) , sigsqrtt);
        Nd1 = normcdf(d);
        fv = av .* Nd1 - discfact .* debt .* normcdf(bsxfun( @minus,d , sigsqrtt));
        
        zz = (fv<eqval);
                  
        if all(all((zz == zeroMTX)))
            d_upper = av;
            success_u = 1;
        else
            d_upper = av + zz .* av;
            av = d_upper;
        end;
    end;
      
    %decrease the lower bound until the model equity value is less than eqval
    success_l = 0;
    av = d_lower;
    while (success_l ~= 1)                       %set the lower bound
        d = bsxfun( @rdivide, (log(av ./ debt) + rsigsq), sigsqrtt );
        Nd1 = normcdf(d);
        fv = av .* Nd1 - discfact .* debt .* normcdf(bsxfun( @minus,d , sigsqrtt));
        zz = (fv > eqval);
        
        if all(all(zz == zeroMTX));
            d_lower = av;
            success_l = 1;
        else
            d_lower = av - 0.5 * zz .* av;
            av = d_lower;
        end;
    end;
     
    av = (d_upper + d_lower) / 2;
end

for iter = 1:maxiter_newton 
    
    d1 = (log(av ./ debt) + rsigsq) ./ sigsqrtt;
    d2 = d1 - sigsqrtt; 
    
    % compute f'(av)
    
    Nd1 = normcdf(d1);
    Nd2 = normcdf(d2);
    
    deriv = Nd1; 
    
    fv = av .* Nd1 - discfact .* debt .* Nd2; 
      
    % compute the step to take from av 
    
    av_step = -(fv - eqval) ./ deriv;  
    
    
    if any(any(isinf(av_step)))==1
        break;
    end
    
    if  max(max(abs(av_step)./av)) < tol 
                      
        return
    else
        av = av + av_step; 
    end
end


% If Newton's method fails, switch to bi-section method
% disp('Using Bisection')
if ~isnan(av_ini)
    
    d_upper = 0.75 * (eqval + debt);
    d_lower = 0.8 * d_upper;
    
    %increase the upper bound until the model equity value is greater than eqval
    success_u = 0;
    av = d_upper; %initialize asset value
    
    while (success_u~=1) %set the upper bound
        d = (log(av./debt) + rsigsq) ./ sigsqrtt;
        Nd1 = normcdf(d);
        fv = av .* Nd1 - discfact .* debt .* normcdf(d - sigsqrtt);
        
        zz = (fv<eqval);
        
        if all(all(zz ==zeroMTX))
            d_upper = av;
            success_u = 1;
        else
            d_upper = av + zz .* av;
            av = d_upper;
        end
    end
    
    %decrease the lower bound until the model equity value is less than eqval
    success_l = 0;  %set the lower bound
    av = d_lower;
    while (success_l ~= 1)
        d = (log(av ./ debt) + rsigsq) ./ sigsqrtt;
        Nd1 = normcdf(d);
        fv = av .* Nd1 - discfact .* debt .* normcdf(d - sigsqrtt);
        zz = (fv > eqval);
        
        if all(all(zz == zeroMTX));
            d_lower = av;
            success_l = 1;
        else
            d_lower = av - 0.5 * zz .* av;
            av = d_lower;
        end
    end
end


for iter = 1:maxiter_bisection
    
    av = (d_upper + d_lower)/2;
    
    d = (log(av ./ debt) + rsigsq) ./ sigsqrtt;
    Nd1 = normcdf(d);
    fv = av .* Nd1 - discfact .* debt .* normcdf(d - sigsqrtt);
    
    zz = (fv > eqval);
    
    d_upper = zz .* av + (1 - zz) .* d_upper;
    d_lower = zz .* d_lower + (1 - zz) .* av;
    
    if (max(max((d_upper - d_lower) ./ d_lower)) <= tol);
        av = (d_upper + d_lower) /2;
        return
    end
end
 
 

