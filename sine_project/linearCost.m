function [J, grad] = linearCost(theta, X, y, alpha, lambda)
% Function parameters:
% X      => mxn matrix with m training examples and n features
%          First column of X should be 1
% theta  => nx1 matrix
% y      => mx1 matrix
% alpha  => learning rate

% Output variables
% J      => Cost as a real number
% grad   => nx1 vector (Use theta = theta - grad for descent)

J = 0;

m = size(X,1);
n = size(X,2);

if size(theta,1) != n || size(theta,2) != 1,
  fprintf('ERROR: theta should be %dx%d, found %dx%d\n', n, 1, size(theta, 1), size(theta, 2));
  pause;
  return;
end

if size(y,1) != m || size(y,2) != 1,
  fprintf('ERROR: y should be %dx%d, found %dx%d\n', m, 1, size(y, 1), size(y, 2));
  pause;
  return;
end

% Validations complete

% Calculate hypothesis
% h(x) = X * Theta
h = X * theta; 

% Calculate cost
% J = 1/m * SUM ((h(x) - y)^2) + (lambda/m) * SUM [theta .* theta]
J = (1/m) * sum( (h - y) .^ 2);

% Regularize J => + (lambda/m) * (sum(theta^2) => not counting theta0)
J = J + ((lambda/m) * sum(theta .* theta));
J = J - ((lambda/m) * theta(1) * theta(1));  % Cancel theta0

% Calculate gradient descent
% Theta_j = Theta_j - (alpha/m) * SUM(1..m) [ (h(x) - y) * x_j ]
grad = (alpha / m ) * (X' * (h - y));

% Regularize grad
% + (lambda/m) * theta_j  (not applying for theta0)
temp = grad(1);
grad = grad + ((lambda * alpha/m) * theta);
grad(1) = temp;


