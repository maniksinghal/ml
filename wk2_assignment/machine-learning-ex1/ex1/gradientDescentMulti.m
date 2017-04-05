function [theta, J_history] = gradientDescentMulti(X, y, theta, alpha, num_iters)
%GRADIENTDESCENTMULTI Performs gradient descent to learn theta
%   theta = GRADIENTDESCENTMULTI(x, y, theta, alpha, num_iters) updates theta by
%   taking num_iters gradient steps with learning rate alpha

% Initialize some useful values
m = length(y); % number of training examples
J_history = zeros(num_iters, 1);

for iter = 1:num_iters

    % ====================== YOUR CODE HERE ======================
    % Instructions: Perform a single gradient step on the parameter vector
    %               theta. 
    %
    % Hint: While debugging, it can be useful to print out the values
    %       of the cost function (computeCost) and gradient here.
    %
    theta_prev = theta;
    theta = theta - ((X' * (X * theta - y)) * (alpha/m));

    % ============================================================

    % Save the cost J in every iteration    
    J_history(iter) = computeCost(X, y, theta);

    fprintf('Cost-function now [%d]: %f\n', iter, J_history(iter))

    if iter > 1,
    	if J_history(iter) > J_history(iter-1),
		% Cost function started increasing, stop here
		fprintf('Cost function started to increase: %f -> %f, Stop here', J_history(iter-1), J_history(iter))
		theta = theta_prev
		break
	end
    end

end

end
