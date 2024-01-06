clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = a_(n-1) * H_n / H_(n-1).
% H_n=psi(n+1)-psi(1) is the n-th term of the Harmonic sequence.
% The solution is obivously H_n itself.
% The guesses are x!, H_x, 1, x, ..., x^5 and pairwise products, such as 
% x^2*H_x.
syms x;
functions = {factorial(sym(x)), psi(sym(x)+1)- psi(sym(1))};
precision = 200;
start = 1;
rec_degree = 2;
poly_degree = 5;
[func_name, solutions] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true, false); 
% solutions has for each solution a vector with the coefficients
% func_name has the corresponding functions
disp(solutions)
disp(func_name)

function ret = substitute(n, func_val)
    % ret = a_n - a_(n-1) * H_n / H_(n-1)
    h = (psi(n+1)- psi(sym(1))) / (psi(n)- psi(sym(1)));
    ret = func_val(:,1) - h * func_val(:,2);
end
