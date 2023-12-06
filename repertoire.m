function solution = repertoire(func_guess, recursive_relation, nonHomogeneous, precision, start_point, rec_degree, num_poly, varargin)
    digits(precision);
    num_func = length(func_guess);
    % The number of required data points
    % 1 + 2p + p*f where p is the number of polynomials and f is the number
    % of functions
    num_substitution = 1 + num_poly * ((num_func + num_func^2) / 2 + 1);
    disp(num_substitution);
    num_function_values = num_substitution + rec_degree;
    % ---------------------Calculating function-values---------------------
    range = num2cell(start_point:start_point+num_function_values-1);
    verbose = false;
    checkNull = false;
    nVarargs = length(varargin);
    if nVarargs > 0
        checkNull = varargin{1};
    end
    if nVarargs > 1
        verbose = varargin{2};
    end
    % A matrix where enties are the input functions evaluated at the points
    [function_values, function_name] = generateGuesses(func_guess, range, num_poly, verbose);
    %--------------------Checking best trivial solution--------------------
    % Will probably change, is not relavant right now
    if checkNull
        if verbose
            fprintf("Checking for null in input");
        end
        A = function_values';
        column_norms = sqrt(sum(A.^2, 1));
        A = A ./ column_norms;
        [~, ~, V2] = svd(A);
        x2 = V2(:, end);
        y2 = A * x2;
        n_best = norm(y2, 'fro');
        if log10(n_best) < -precision
            disp("Found a null!");
            disp("Consider raising the precision to at least " + num2str(ceil(-double(log10(n_best)))) + " if this is not a null.")
        end

        if verbose
            fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
            fprintf("All null eleminated\n")
        end
    end
    %-------Substituting function values into the recursive equation-------
    [N, M] = size(function_values);
    M = M-rec_degree;
    substituted = sym(zeros(N, M));
    if verbose
        fprintf("Calculating substitution")
    end
    idx = 1:rec_degree;
    for i = 1:M
        n = (sym(i)-2)+rec_degree + start_point;
        s = recursive_relation(n, function_values(:,i+(rec_degree-idx)));
        substituted(:, i) = s;
    end
    if verbose
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
        disp("Finished calulating substitution in recursive equation");
        fprintf("Calculating results");
    end
    %--------------------------------Result--------------------------------
    % Finding a vector that minimizes the matrix of the substituted values
    matrix = substituted';
    % Adding a column to the matrix in case of a nonhomogeneous equation
    if nonHomogeneous ~= 0
        nonHColumn = zeros(M, 1);
        for i = 1:M
            n = sym(i)-1+rec_degree;
            nonHColumn(i) = subs(nonHomogeneous, x, n);
        end
        matrix = [matrix, nonHColumn];
    end
    [~, ~, V] = svd(matrix);
    x = V(:, end);
    column_norms = sqrt(sum(matrix.^2, 1)); 
    [value, index] = min(column_norms);
    if log10(value) < -precision/2
        solution = zeros(size(matrix, 2), 1);
        solution(index) = 1;
        if verbose
            fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            disp("---------------------Results---------------------")
            if checkNull
                disp("How good is the best null? (larger values are better) " + sprintf('%.2e', n_best));
            end
            disp("How good is the solution? (smaller values are better) " + sprintf('%.2e', norm(matrix*solution, 'fro')));
            disp("How close is the solution to null? (larger values are better) " +  sprintf('%.2e', norm(function_values' * solution, 'fro')));
            disp(function_name(index));
        end
        return;
    end
    x_normalized = x .* column_norms.';  
    selectedColumns = [];
    remainingColumns = 1:size(matrix, 2);
    
    for i=1:size(matrix, 2)
        [~, max_index] = max(abs(x_normalized(remainingColumns)));
        most_influential_index = remainingColumns(max_index);
        remainingColumns(max_index) = [];
        selectedColumns = [selectedColumns, most_influential_index];
        reduced_matrix = matrix(:, selectedColumns);
        [~, ~, candidate_V] = svd(reduced_matrix);
        candidate_x = candidate_V(:, end);
        result = log10(norm(reduced_matrix * candidate_x, 'fro'));
        if result < -precision/2
            break;
        end
    end
    x = zeros(size(matrix, 2), 1);
    x(selectedColumns) = candidate_x;
    solution = x;
    y = matrix * x;
    if nonHomogeneous ~= 0
        x = x(1:end-1);
    end
    y3 = function_values' * x;
    example = 0;
    x = x ./ candidate_x(1);
    for i=selectedColumns
        formatted_x = sprintf('%.2e', x(i));
        example = example + function_name{i} * str2double(formatted_x);
    end
    if verbose
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        disp("---------------------Results---------------------")
        if checkNull
            disp("How good is the best null? (larger values are better) " + sprintf('%.2e', n_best));
        end
        disp("How good is the solution? (smaller values are better) " + sprintf('%.2e', norm(y, 'fro')));
        disp("How close is the solution to null? (larger values are better) " +  sprintf('%.2e', norm(y3, 'fro')));
        disp(example);
    end
end