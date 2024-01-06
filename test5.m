clc;
clearvars;

syms x;
functions = {factorial(sym(x)), 3^sym(x), (-1)^sym(x), (psi(sym(x)+1)- psi(sym(1))), sym(x), sym(x)^2, sym(x)^3};
precision = 30; 
start = 0;
rec_degree = 3;
poly_degree = 1;
solution = repertoire(functions, @substitute, 1, precision, start, rec_degree, poly_degree, true);

function ret = substitute(~, func_val)
    ret = func_val(:,1) - 4 * func_val(:,2) + 3 * func_val(:, 3);
end
