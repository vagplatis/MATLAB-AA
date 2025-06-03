function erg2_3
    % Ergasia 2 - 05_2025

    clc;

    f = @(x) x^2 - 4*sin(x);
    df = @(x) 2*x - 4*cos(x);

    tol = 10^-7;
    max_it = 50;

    a = 1; b = 6; i = 1;
    fprintf('\n\n\n');
    fprintf('---- Bisect ----\n'); bisect(f,a,b,tol,max_it);
    fprintf('---- Bisect Rec ----\n'); bisect_r(f,a,b,i,tol,max_it);

    x0 = 6; i = 1;
    fprintf('\n\n\n');
    fprintf('--- NewtRaph ---\n'); newton_raphson(f,df,x0,tol,max_it);
    fprintf('--- NewtRaph Rec ---\n'); newton_raphson_r(f,df,x0,i,tol,max_it);

    x0 = 6; x1 = 7; i = 1;
    fprintf('\n\n\n');
    fprintf('---- Secant ----\n'); secant(f,x0,x1,tol,max_it);
    fprintf('---- Secant Rec ----\n'); secant_r(f,x0,x1,i,tol,max_it);
end

function x = bisect(f, a, b, tol, max_it)
% Description: Finds a root of f in the interval [a, b] with precision tol
% using the bisect method
% Input: 
%   - f: the function
%   - a: the lower value of the interval
%   - b: the upper value of the interval
%   = tol: the tolerance
%   - max_it: the maximum number of iterations
% Output: 
%   - x: the root of f
% Author: Markos Tsipouras

    fa = feval(f, a);
    fb = feval(f, b);
    
    i = 1;
    x_previous = b;
	while i <= max_it
        x = (a + b) / 2;
        fprintf('Iteration %3d: %10.8f\n', i, x);
        i = i + 1;
        fx = feval(f, x);
        if fx == 0
            return;
        end
        if fb * fx < 0
            a = x;
            fa = fx;
        else
            b = x;
            fb = fx;
        end
        if abs(x - x_previous) < tol
            return;
        end
        x_previous = x;
    end
end

function x = bisect_r(f, a, b, i, tol, max_it)
% Description: Finds a root of f in the interval [a, b] with precision tol
% using the recurrent bisect method
% Input: 
%   - f: the function
%   - a: the lower value of the interval
%   - b: the upper value of the interval
%   - i: num of iteration 
%   = tol: the tolerance
%   - max_it: the maximum number of iterations
% Output: 
%   - x: the root of f
% Author: omadara 3

    x =(a+b)/2;
    fprintf('Iteration %3d: %10.8f\n', i, x);%meta to x gia na pairnei thn timh
    fa=feval(f,a);%den xreiazetai alla ama thes mporeis na kaneis kai fa*fx<0 kai to allazeis analoga sthn anadromh
    fb=feval(f,b);
    fx=feval(f,x);

    if fx == 0 || abs((a-b)/2) < tol || i>=max_it %(a-b)/2 giati h riza einai mesa sto diasthma [a,b]
        return;
    end

    if fb * fx < 0
        x = bisect_r(f, x, b, i+1, tol, max_it);% a=x
    else
        x = bisect_r(f, a, x, i+1, tol, max_it);% b=x
    end
end

function x = newton_raphson(f, df, x0, tol, max_it)
% Description: Finds a root of f starting from x0 with precision tol
% using the Newton-Raphson method
% Input: 
%   - f: the function
%   - df: the derivative of the function
%   - x0: the initial guess
%   = tol: the tolerance
%   - max_it: the maximum number of iterations
% Output: 
%   - x: the root of f
% Author: Markos Tsipouras

    i = 1;
    x = x0;
    x_previous = x0;
    while i <= max_it
        dx = feval(f, x_previous) / feval(df, x_previous);
        x = x_previous - dx;
        fprintf('Iteration %3d: %10.8f\n', i, x);
        i = i + 1;
        if abs(x - x_previous) < tol
            return;
        end
        x_previous = x;
    end
end

function x = newton_raphson_r(f, df, x0, i, tol, max_it)
% Description: Finds a root of f starting from x0 with precision tol
% using the recurrent Newton-Raphson method
% Input:
%   - f: the function
%   - df: the derivative of the function
%   - x0: the initial guess
%   - i: num of iteration 
%   = tol: the tolerance
%   - max_it: the maximum number of iterations
% Output:
%   - x: the root of f
% Author: pio omadara 3
    dx = feval(f, x0) / feval(df, x0);    
    x = x0 - dx;
    fprintf('Iteration %3d: %10.8f\n', i, x);
    if i>=max_it || abs(x - x0) < tol %giati x-x0 einai to h ara oso pio mikro toso pio konta sthn riza
        return
    end
     x = newton_raphson_r(f, df, x, i + 1, tol, max_it);
end

function x = secant(f, x0, x1, tol, max_it)
% Description: Finds a root of f starting from x0, x1 with precision tol
% using the secant method
% Input: 
%   - f: the function
%   - x0: the initial guess
%   - x1: the initial guess
%   = tol: the tolerance
%   - max_it: the maximum number of iterations
% Output: 
%   - x: the root of f
% Author: Markos Tsipouras

    i = 1;
	while i <= max_it
        f0 = feval(f, x0);
        f1 = feval(f, x1);
        x = x1 - f1 * (x1 - x0) / (f1 - f0);
        fprintf('Iteration %3d: %10.8f\n', i, x);
        i = i + 1;
        if abs(x - x1) < tol
            return;
        end
        x0 = x1;
        x1 = x;
    end
end

function x = secant_r(f, x0, x1, i, tol, max_it)
% Description: Finds a root of f starting from x0, x1 with precision tol
% using the recurrent secant method
% Input: 
%   - f: the function
%   - x0: the initial guess
%   - x1: the initial guess
%   - i: num of iteration 
%   = tol: the tolerance
%   - max_it: the maximum number of iterations
% Output: 
%   - x: the root of f
% Author: min ta ksanaleme
    f0 = feval(f, x0);
    f1 = feval(f, x1);
    x = x1 - f1 * (x1 - x0) / (f1 - f0);
    fprintf('Iteration %3d: %10.8f\n', i, x);

    if i>=max_it || abs(x - x1) < tol%pali opws sthn NR x-x1 deixnei metavolh sthn proseggish
        return
    end
    x = secant_r(f, x1, x, i+1, tol, max_it); %x0=x1 , x1=x
end