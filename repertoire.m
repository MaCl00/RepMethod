function solution = repertoire(func_guess, recursive_relation, precision, start_point, rec_degree, num_poly, verbose)
    digits(precision);
    syms x;
    
    num_func = length(func_guess);
    
    % The number of required data points
    % 1 + 2p + p*f where p is the number of polynomials and f is the number
    % of functions
    num_function_values = 1 + num_poly * (2 + num_func);

    % Calculating function-values
    range = num2cell(start_point:start_point+num_function_values+rec_degree-1);
    sequences = cellfun(@(j) cellfun(@(i) subs(func_guess{j}, x, sym(i)), range), num2cell(1:num_func), 'UniformOutput', false);
    N = length(sequences);
    result = cell(1, (N*N-N)/2);
    counter = 1;
    prod_func = cell(1, (N*N-N)/2);
    for i = 1:N
        for j = i+1:N
             result{counter} = sequences{i} .* sequences{j};
             prod_func{counter} = func_guess{i} * func_guess{j};
             counter = counter + 1;
        end
    end
    sequences = [sequences, result];
    func_guess = [func_guess, prod_func];
    polynomes = cell(1, num_poly);
    poly_name = cell(1, num_poly);
    for j=0:num_poly-1
        polynomes{j+1} = cellfun(@(i) (sym(i)/fraction_value)^j, range);
        poly_name{j+1} = sym(x)^j;
    end
    N = length(sequences);
    result = cell(1, N*num_poly);
    counter = 1;
    function_name = cell(1, N*num_poly);
    for i = 1:N
        for j= 1:num_poly
            result{counter} = sequences{i} .* polynomes{j};
            function_name{counter} = func_guess{i} *poly_name{j};
            counter = counter + 1;
        end
    end
    function_name = [function_name, poly_name];
    sequences = [result, polynomes];
    fprintf("Finished calulating function-values\n");

    % Substituting function values into the recursive equation
    matrix2 = vertcat(sequences{:})';
    N = length(sequences);
    M = length(sequences{1})-rec_degree*fraction_value;
    substituted = sym(zeros(N, M));
    fprintf("Calculating substitution(")
    for i = 1:M
        fprintf('%5d of %5d)\n', i, M);
        for j = 1:N
            n = (sym(i)-1)/fraction_value+rec_degree;
            s = recursive_relation(n, arrayfun(@(k) sequences{j}(i+(rec_degree-k)*fraction_value), 0:rec_degree));
            substituted(j, i) = s;
        end
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    end
    fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
    disp("Finished calulating substitution in recursive equation");

    % Finding a vector that minimizes the matrix of the substituted values
    matrix = substituted';
    [~, ~, V] = svd(matrix);
    x = V(:, end);
    y = matrix * x;
    y3 = matrix2 * x;
    if verbose
        [~, ~, V2] = svd(matrix2);
        x2 = V2(:, end);
        y2 = matrix2 * x2;
        disp("How good is the best null? (larger values are better) " + string(norm(y2, 'fro')));
        disp("How good is the solution? (smaller values are better) " + string(norm(y, 'fro')));
        disp("How close is the solution to null? (larger values are better) " + string(norm(y3, 'fro')));
    end
    solution = 0;
    for i=1:length(function_name)
        if abs(norm(x(i)*sequences(i), 'fro')) > 0.1
            solution = solution + function_name{i} * x(i);
        end
    end
    
    disp(solution);
end