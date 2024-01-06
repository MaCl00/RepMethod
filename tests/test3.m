clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = 4*(n+1)/n*a_(n-1)+4*(n+1)/(n-1).
% The solution are (1+x)*2^x and (x+x^2)*2^x.
% The guesses are x!, 2^x, 1, x, ..., x^6 and pairwise products, such as 
% x^2*H_x.
syms x;
functions = {gamma(sym(x)+1), 2^sym(x), 3^sym(x), 6^sym(x)};
precision = 300;
start = 0;
rec_degree = 3;
poly_degree = 5;
[func_name, solutions] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true, false); 
% solutions has for each solution a vector with the coefficients
% func_name has the corresponding functions
disp(solutions)
disp(func_name)

function ret = substitute(n, func_val)
    ret = func_val(:,1) - 4/sym(n)*sym(n+1) * func_val(:,2) + 4/sym(n-1)*sym(n+1) * func_val(:,3);
end
