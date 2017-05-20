function [J, thetaMin, J_slope] = minimizeLinearCost(f, X, theta, y, alpha, lambda, iterations)
% Minimize the cost calculated by the input function f by iterating in a loop and applying
% Gradient descent.
%
% Function parameters:
% f      => Cost function which must return [J, grad]
% X      => mxn matrix with m training examples and n features
%          First column of X should be 1
% theta  => nx1 matrix
% y      => mx1 matrix
% alpha  => learning rate
% lambda => regularization parameter

% Output variables
% J        => Minimized Cost as a real number
% thetaMin => nx1 matrix of final theta
% J_slope  => Nx1 matrix (N: Number of iterations to find min cost) of how J descent

m = size(X,1);
n = size(X,2);

J = 0;
J_slope = zeros(iterations, 1);
thetaMin = theta;

%Verify that first column of input should be all ones
if sum(X(:,1)) != m,
  fprintf('ERROR[minimizeLinearCost]: First column of X should be a bias vector (all ones)\n');
  pause;
  return;
end

% Theta should be nx1 (where n is number of features in X)
if size(theta,1) != n || size(theta,2) != 1,
  fprintf('ERROR: theta should be %dx%d, found %dx%d\n', n, 1, size(theta, 1), size(theta, 2));
  pause;
  return;
end

% Y should be mx1 where m is the number of training examples
if size(y,1) != m || size(y,2) != 1,
  fprintf('ERROR: y should be %dx%d, found %dx%d\n', m, 1, size(y, 1), size(y, 2));
  pause;
  return;
end

% Validations complete

i = 0;
theta_prev = thetaMin;
J_prev = 0;
while i < iterations,
	[J, grad] = f(thetaMin, X, y, alpha, lambda);
	if i > 0,
		if J > J_prev,
			% We have started to overshoot, stop here
			break;
		end
	end

	J_prev = J;
	J_slope(i+1) = J_prev;
	theta_prev = thetaMin;

	thetaMin = thetaMin - grad;
	i = i + 1;
		
end

thetaMin = theta_prev;
J = J_prev;
J_slope = J_slope(1:i,1);

