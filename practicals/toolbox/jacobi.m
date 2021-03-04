function yy = jacobi(A, b, FGS, abs_tol, maxiter, SOR, show_error, for_loops)
    % Solves a linear system using the Jacobi method
    % with Successive-Over-Relaxation
    %
    % Arguments:
    % - A: (N x N) the system square matrix
    % - b: (N x 1)the known vector
    % - FGS: (N x 1) a First-Guess-Solution vector
    %
    % Optional arguments:
    % - abs_tol: (float) absolute tolerance
    % - maxiter: (int) maximum number of iterations
    % - SOR: (float) Successive-Over-Relaxation factor (1 < SOR < 2)
    % - show_error: (bool) print error on screen
    % - for_loops: (bool) if set to 1 uses for loops instead of vector
    %   computations (it might be faster in some cases)
    
    if ~exist('abs_tol','var')
        abs_tol = 1e-6;      % Absolute Tolerance
    end
    
    if ~exist('maxiter','var')
        maxiter = 1000000;   % Maximum number of iterations
    end
    
    if ~exist('show_error','var')
        show_error = 1;     % Print error on screen
    end
    
    if ~exist('SOR','var')
        SOR = 1.5;          % Successive-Over-Relaxation factor
    end
    
    if ~exist('for_loops','var')
        for_loops = 0;      % Use for loops instead of vector algebra
    end
    
    d = diag(A);        % Vector containing the diagonal elements of A
    N = length(b);      % Length of b
    B = A - eye(N).*d;  % A matrix without the diagonal
    
    yy = FGS;
    
    % Iterative solution
    for ii = 1:maxiter
        
        if for_loops
            for jj = 1:N

                % Update the vector value
                yy(jj) = SOR*(b(jj) - B(jj, :)*yy)./d(jj) + (1- SOR)*yy(jj);

            end
        else
            
            % Update the vector value
            yy = SOR*(b - B*yy)./d + (1 - SOR)*yy;
            
        end
        
        % Evaluate the error
        error = sum(abs((A*yy - b)));
     
        % Print the error every 100 iterations
        if show_error
            if mod(ii, 100) == 0
                fprintf('iter: %d, abs_error = %f\n', ii, error)
            end
        end
        
        % Check if the error is small enough
        if error <= abs_tol
            break
        end
        
    end
    
    % Raise a warning if the maximum number of iterations was reaches
    if ii == maxiter
        disp('Warning: jacobi reached the maximum number of iterations')
    end
    
end