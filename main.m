clc;
clearvars;
syms x;
functions = {factorial(sym(x)-1)};
precision = 500; 
start = 1;
rec_degree = 3;
poly_degree = 6;
[function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, "verbose");
disp(function_name)
disp(solution)
function ret = substitute(n, func_val)
    ret = func_val(:,1) - (sym(n)+1)*func_val(:,2)+(n-2)*func_val(:,3);
end
