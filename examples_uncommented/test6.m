clc;
clearvars;
syms x;
functions = {factorial(sym(x))};
precision = 100; 
start = 0;
rec_degree = 3;
poly_degree = 2;
[function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true, false);
function ret = substitute(n, func_val)
    ret = (4*sym(n)^2-8*sym(n)+3)*func_val(:,1) - 2*(5*sym(n)^2+1) * func_val(:,2) + 4*(sym(n)^2-sym(n)+1) * func_val(:,3);
end
