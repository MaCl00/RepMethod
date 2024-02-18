clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = n * (n-1) * a_(n-2).
% The solution is n!.
% The guesses are x!, H_x, 1, x, ..., x^3 and pairwise products, such as 
% x^2*H_x.
syms x;
functions = {gamma(sym(x)+1), psi(sym(x)+1)- psi(sym(1))};
precision = 200;
start = 0;
rec_degree = 3;
poly_degree = 3;
[func_name, solutions] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, "verbose", "quitSearch", "addMul"); 
% solutions has for each solution a vector with the coefficients
% func_name has the corresponding functions
disp(solutions)
disp(func_name)

function ret = substitute(n, func_val)
    % ret = a_n - n*(n-1)*a_(n-2)
    ret = func_val(:,1) - n * (n-1) * func_val(:,3);
end