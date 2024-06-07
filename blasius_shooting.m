function [f, f_prime, f_double_prime] = blasius_shooting(eta, scheme)
    % Shooting method to solve Blasius equation
    
    % ODE function
    blasius_ode = @(eta, F) [F(2); F(3); -0.5*F(1)*F(3)];
    
    % Boundary conditions
    f_0 = 0;
    f_prime_0 = 0;
    f_prime_inf = 1;
    
    % Initial guess for f''(0)
    f_double_prime_guess = 0.469600;  % Literature value
    
    % Tolerance for shooting method
    tol = 1e-6;
    
    % Maximum iterations
    max_iterations = 100;
    
    % Shooting method
    for iteration = 1:max_iterations
        fprintf('Iteration %d, f_double_prime_guess = %.6f\n', iteration, f_double_prime_guess);
        
        % Solve ODE
        [~, F] = ode45(blasius_ode, eta, [f_0; f_prime_0; f_double_prime_guess]);
        
        % Extract solutions
        f = F(:, 1);
        f_prime = F(:, 2);
        f_double_prime = F(:, 3);
        
        % Check if boundary condition at infinity is satisfied
        error = abs(f_prime(end) - f_prime_inf);
        fprintf('  f_prime(inf) = %.6f, error = %.6f\n', f_prime(end), error);
        
        if error < tol
            fprintf('Converged in %d iterations!\n', iteration);
            break;
        end
        
        % Adjust guess using secant method
        f_double_prime_guess = adjust_guess(f_double_prime_guess, f_prime(end), scheme);
    end
    
    if iteration == max_iterations
        warning('Did not converge in %d iterations!', max_iterations);
    end
end

function new_guess = adjust_guess(old_guess, f_prime_inf, scheme)
    % Secant method to adjust the guess
    persistent prev_guess prev_error;
    
    error = f_prime_inf - 1;
    
    if isempty(prev_guess) || isnan(prev_guess)
        % First call or after a NaN, use a good initial guess
        new_guess = 0.332047;  % Known good value for f''(0)
        fprintf('  Using known good guess: %.6f\n', new_guess);
    elseif abs(error - prev_error) < 1e-10
        % Avoid division by near-zero
        new_guess = old_guess * 0.9;  % Simple decrease
        fprintf('  Avoiding division by zero. New guess: %.6f\n', new_guess);
    else
        % Secant method
        new_guess = old_guess - error * (old_guess - prev_guess) / (error - prev_error);
        fprintf('  Secant method. New guess: %.6f\n', new_guess);
    end
    
    % Update for next call
    prev_guess = old_guess;
    prev_error = error;
end
