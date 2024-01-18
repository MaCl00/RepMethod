function [solution, solution_found] = SolutionExtraction(matrix, x, precision, checkSolution, verbose)
    selectedColumns = [];
    m = size(matrix, 2);
    remainingColumns = 1:m;
    solution_found = false;
    for i=1:m
        if verbose
            fprintf('Completed %d/%d', i, m);
        end
        [~, max_index] = max(abs(x(remainingColumns)));
        most_influential_index = remainingColumns(max_index);
        remainingColumns(max_index) = [];
        selectedColumns = [selectedColumns, most_influential_index];
        reduced_matrix = matrix(:, selectedColumns);
        [~, ~, candidate_V] = svd(reduced_matrix);
        candidate_x = candidate_V(:, end);
        result = log10(norm(reduced_matrix * candidate_x, 'fro'));
        if verbose
            fprintf(repmat('\b', 1, strlength("Completed /"+num2str(i)+num2str(m))));
        end
        if result < -precision/2
            solution_found = true;
            break;
        end
    end
    solution = zeros(m, 1);
    solution(selectedColumns) = candidate_x;
    if (~isequal(checkSolution, 0) && size(candidate_x, 1) == m) || size(candidate_x, 1) == 1 || ~solution_found
        if verbose
            fprintf("Finished search\n");
        end
        return;
    end
    [solution(selectedColumns), ~] = SolutionExtraction(reduced_matrix, candidate_x, precision, solution, verbose);
end