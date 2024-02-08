clc;
clearvars;
syms x;
functions = {(-1)^sym(x)};
start = 0;
rec_degree = 3;
poly_degree = 4;
precision = 2000;
[function_name, solution] = repertoire(functions, @substitute, sym(x)^2+sym(x)+3, precision, start, rec_degree, poly_degree, false, true);
for i=201:300
    precision = 10*i; 
    disp(precision);
    tic
    [function_name, solution] = repertoire(functions, @substitute, sym(x)^2+sym(x)+3, precision, start, rec_degree, poly_degree, false, true);
    toc
end
function ret = substitute(n, func_val)
    ret = func_val(:,1) - (-1)^sym(n) * func_val(:,2)+sym(n)*func_val(:, 3);
end
