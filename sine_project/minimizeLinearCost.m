function [J, thetaMin, J_slope] = minimizeLinearCost(X, theta, y, alpha, lambda, iterations)
% Function parameters:
% X      => mxn matrix with m training examples and n features
%          First column of X should be 1
% theta  => nx1 matrix
% y      => mx1 matrix
% alpha  => learning rate
% lambda => regularization parameter

% Output variables
% J        => Minimized Cost as a real number
% thetaMin => nx1 matrix of final theta
% J_slope  => Nx1 matrix (N, upto iterations) of how J descent

m = size(X,1);
n = size(X,2);

J = 0;
J_slope = zeros(iterations, 1);
thetaMin = theta;


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

i = 0;
theta_prev = thetaMin;
J_prev = 0;
while i < iterations,
	[J, grad] = linearCost(thetaMin, X, y, alpha, lambda);
	%fprintf('Running iteration %d, cost: %f\n', i, J);
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

