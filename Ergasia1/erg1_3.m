% ergasia.m
% ΟΜΑΔΑ: 01
% ΟΝΟΜΑΤΕΠΩΝΥΜΑΤΑ: ΠΛΑΤΗΣ ΕΥΑΓΓΕΛΟΣ, ΞΥΛΟΥΡΗΣ ΕΜΜΑΝΟΥΗΛ
% ΑΜ: 2229, 2216
% Μάθημα: Αριθμητική Ανάλυση - Εργασία 1
% Θέμα: Αλγόριθμος PTRANS-II για πενταδιαγώνια συστήματα
function erg1_3()
    erg1();
end
%%Συνάρτηση  κατασκευής πενταδιαγώνιου πίνακα
function p = pendatiagonal(e,c,d,a,b)
    N = length(d);
    p = diag(d);
    if N > 1
        p = p + diag(c, -1) + diag(a, 1);
    end
    if N > 2
        p = p + diag(e, -2) + diag(b, 2);
    end
end
%%Συνάρτηση PTRANSII
function [x, psi] = PTRANSII(N, e, c, d, a, b, y)
    if N < 4
        error('Η PTRANSII απαιτεί N >= 4.');
    end
    eta = zeros(N, 1);
    sigma = zeros(N, 1);
    phi = zeros(N, 1);
    psi = zeros(N, 1);
    rho = zeros(N, 1);
    w = zeros(N, 1);
    x = zeros(N, 1);
    psi(N) = d(N);
    if N > 1
        sigma(N) = c(N-1) / psi(N);
    end
    if N > 2
        phi(N) = e(N-2) / psi(N);
    end
    w(N) = y(N) / psi(N);
    rho(N-1) = a(N-1);
    psi(N-1) = d(N-1) - sigma(N) * rho(N-1);
    sigma(N-1) = (c(N-1) - phi(N) * rho(N-1)) / psi(N-1);
    phi(N-1) = e(N-1) / psi(N-1);
    w(N-1) = (y(N-1) - w(N) * rho(N-1)) / psi(N-1);
    
    for i = N-2:-1:3
        rho(i) = a(i) - sigma(i+2) * b(i);
        psi(i) = d(i) - phi(i+2) * b(i) - sigma(i+1) * rho(i);
        sigma(i) = (c(i) - phi(i+1) * rho(i)) / psi(i);
        phi(i) = e(i) / psi(i);
        w(i) = (y(i) - w(i+2) * b(i) - w(i+1) * rho(i)) / psi(i);
    end
    rho(2) = a(2) - sigma(4)*b(2);
    psi(2) = d(2) - phi(4)*b(2) - sigma(3)*rho(2);
    sigma(2) = (c(2) - phi(4)*rho(2)) / psi(2);
    rho(1) = a(1) - sigma(3)*b(1);
    psi(1) = d(1) - phi(3)*b(1) - sigma(2)*rho(1);
    w(2) = (y(2) - w(4)*b(2) - w(3)*rho(2))/psi(2);
    w(1) = (y(1) - w(3)*b(1) - w(2)*rho(1))/psi(1);
    x(1) = w(1);
    x(2) = w(2) - sigma(2)*x(1);
    for i = 3:1:N
        x(i) = w(i) - sigma(i)*x(i - 1) - phi(i)*x(i - 2);
    end    
end
%Συνάρτηση Πειράματος 
function erg1()
    N_small = 4:35;
    results_small = run_experiments(N_small);
    plot_results(N_small, results_small, 'Μικρά Συστήματα (Ν=4:35)');
    N_large = 50:50:2000;
    results_large = run_experiments(N_large);
    plot_results(N_large, results_large, 'Μεγάλα συστήματα (N=50:50:2000)');
end
function results = run_experiments(N_values)
    num_tests = length(N_values);
    results = struct( ...
        'PTRANS_times', zeros(num_tests, 1), ...
        'Cramer_times', NaN(num_tests, 1), ...  % Initialize as NaN
        'Gauss_times',  NaN(num_tests, 1), ...  % Initialize as NaN
        'MATLAB_times', zeros(num_tests, 1) ...
    );
    for i = 1:num_tests
        N = N_values(i);
        fprintf('Testing N = %d...\n', N);
        
        t_PTRANS = 0; 
        t_Cramer = 0; 
        t_Gauss = 0; 
        t_MATLAB = 0;
        
        for r = 1:10  % Average over 10 trials
            % Generate a stable pentadiagonal system (diagonally dominant)
            [e, c, d, a, b, y] = generate_stable_system(N);
            P = pendatiagonal(e, c, d, a, b);
            % --- PTRANS-II ---
            tic;
            x_ptrans = PTRANSII(N, e, c, d, a, b, y);
            t_PTRANS = t_PTRANS + toc;
            % --- MATLAB backslash (reference) ---
            tic;
            x_matlab = P \ y;
            t_MATLAB = t_MATLAB + toc;
            % Validate PTRANS-II against MATLAB
            err = norm(x_ptrans - x_matlab);
            if err > 1e-6
                warning('PTRANS-II error for N=%d: %e', N, err);
            end
            % --- Cramer's Rule (only for tiny systems) ---
            if N <= 10
                tic;
                x_cramer = cramer(P, y);
                t_Cramer = t_Cramer + toc;
                
                % Validate Cramer (optional)
                err_cramer = norm(x_cramer - x_matlab);
                if err_cramer > 1e-6
                    warning('Cramer error for N=%d: %e', N, err_cramer);
                end
            end
            % --- Gaussian Elimination (for small systems) ---
            if N <= 35
                tic;
                x_gauss = gaussianElimination(P, y);
                t_Gauss = t_Gauss + toc;
                
                % Validate Gauss (optional)
                err_gauss = norm(x_gauss - x_matlab);
                if err_gauss > 1e-6
                    warning('Gauss error for N=%d: %e', N, err_gauss);
                end
            end
        end
        % Store average times
        results.PTRANS_times(i) = t_PTRANS / 10;
        results.MATLAB_times(i) = t_MATLAB / 10;
        
        if N <= 10
            results.Cramer_times(i) = t_Cramer / 10;
        end
        if N <= 35
            results.Gauss_times(i) = t_Gauss / 10;
        end
    end
end
function [e, c, d, a, b, y] = generate_stable_system(N)
    if N < 4
        error('Το N πρέπει να είναι τουλάχιστον 4 για PTRANS-II');
    end
    e = randi([1,11], max(0, N-2), 1);
    c = randi([1,11], max(0, N-1), 1);
    d = randi([1,11], N, 1);
    a = randi([1,11], max(0, N-1), 1);
    b = randi([1,11], max(0, N-2), 1);
    y = randi([1,101], N, 1);
end
function x = cramer(A, b)
% Description: Solves Ax = b using the Cramer rule
% Input: 
%   - A: the coefficient matrix
%   - b: the right-hand side values
% Output: 
%   - x: the solution vector
% Author: Markos G. Tsipouras
	[m, n] = size(A);
	if m ~= n
		error('Matrix A must be square!');
	end
	n1 = length(b);
	if n1 ~= n
		error('Vector b should be equal to the number of rows and columns of A!');
	end
	p = det(A); 
    if p == 0
		error('The determinant is zero! The system has "infinitely many solutions" or "no solutions"!');
    end
	x = zeros(n, 1);
	for j = 1:n
		x(j) = det([A(:,1:j-1) b A(:,j+1:end)]) / p;
	end
end
function x = gaussianElimination(A, b)
% Description: Solves Ax = b using the Gaussian elimination with partial pivoting
% Input: 
%   - A: the coefficient matrix
%   - b: the right-hand side values
% Output: 
%   - x: the solution vector
% Author: Markos G. Tsipouras
	[m, n] = size(A);
	if m ~= n
		error('Matrix A must be square!');
	end
	n1 = length(b);
	if n1 ~= n
		error('Vector b should be equal to the number of rows and columns of A!');
	end
	Aug = [A b]; % build the augmented matrix
	C = zeros(1, n + 1);
	
	% elimination phase
	for k = 1:n - 1
		% ensure that the pivoting point is the largest in its column
		[pivot, j] = max(abs(Aug(k:n, k)));
		C = Aug(k, :);
		Aug(k, :) = Aug(j + k - 1, :);
		Aug(j + k - 1, :) = C;
		if Aug(k, k) == 0
			error('Matrix A is singular');
		end
		for i = k + 1:n
			r = Aug(i, k) / Aug(k, k);
			Aug(i, k:n + 1) = Aug(i, k:n + 1) - r * Aug(k, k: n + 1);
		end
	end
	
	% back substitution phase
	x = zeros(n, 1);
	x(n) = Aug(n, n + 1) / Aug(n, n);
	for k = n - 1:-1:1
		x(k) = (Aug(k, n + 1) - Aug(k, k + 1:n) * x(k + 1:n)) / Aug(k, k);
	end
end
function plot_results(N_values, results, title_str)
    figure;
    if max(N_values) <= 35
        semilogy(N_values, results.PTRANS_times, 'b-o', 'DisplayName', 'PTRANS-II');
        hold on;
        semilogy(N_values, results.Cramer_times, 'r-*', 'DisplayName', 'Cramer');
        semilogy(N_values, results.Gauss_times, 'g-s', 'DisplayName', 'Gauss');
        semilogy(N_values, results.MATLAB_times, 'k-d', 'DisplayName', 'MATLAB \\');
    else
        semilogy(N_values, results.PTRANS_times, 'b-o', 'DisplayName', 'PTRANS-II');
        hold on;
        semilogy(N_values, results.MATLAB_times, 'k-d', 'DisplayName', 'MATLAB \\');
    end
    xlabel('Διάσταση συστήματος N');
    ylabel('Μέσος χρόνος εκτέλεσης (s)');
    title(title_str);
    legend('Location','northwest');
    grid on;
    hold off;
end
