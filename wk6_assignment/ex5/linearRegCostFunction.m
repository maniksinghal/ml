function [J, grad] = linearRegCostFunction(X, y, theta, lambda)
%LINEARREGCOSTFUNCTION Compute cost and gradient for regularized linear 
%regression with multiple variables
%   [J, grad] = LINEARREGCOSTFUNCTION(X, y, theta, lambda) computes the 
%   cost of using theta as the parameter for linear regression to fit the 
%   data points in X and y. Returns the cost in J and the gradient in grad

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost and gradient of regularized linear 
%               regression for a particular choice of theta.
%
%               You should set J to the cost and grad to the gradient.
%

% Compute cost:
h = X * theta - y;

% Regularize
reg_sum = theta .* theta;
reg_sum = sum(reg_sum) - reg_sum(1);

J = (sum(h .^ 2) / (2 * m)) + ((reg_sum * lambda) / (2*m));

% Gradient
grad_theta = (lambda / m) * theta;
grad_theta(1) = 0;   % Don't regularize theta0
grad = ((1/m) * (X' * h)) + grad_theta;











% =========================================================================

grad = grad(:);

end
