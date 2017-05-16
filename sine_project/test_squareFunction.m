clear; close all; clc;

% Use a test vector X
X = [0.01 : 0.01 : 5];
X = X';

% y = X^2 + 3
%y = (X .* X) + 3;
y = sin(X);
%y = (X + 3) - X;


% Classic case of overfitting due to too many features
% Use y = 3 (constant)
% and in X, use features like x, x^2, x^3, 1/x 1/x^2, log(x)
% We find the curve starts from very complex shape
% and eventually ends up not matching very well within 10K iterations
%
% At the same time when we reduce the features to
% just x, x^2, the algorithm is able to learn and match the output
% quite well.
%
% Alternatively, use regularization (lambda = 10) and observed that
% the algorithm ends up fitting very nicely

% @todo: 
% Compare fmincg vs usual iteration loop
% Plot with different alpha
% Plot with different lambda

% Add more features to X
Xf = [X, X .* X];   % X = [x, x^2]
Xf = [Xf, X .* X .* X];   % X = [x, x^2, x^3]
Xf = [Xf, log(X)];   % X = [x, x^2, x^3, log(x)]
Xf = [Xf, 1 ./ X];   % X = [x, x^2, x^3, log(x), 1/x]
Xf = [Xf, 1 ./ (X .* X)];   % X = [x, x^2, x^3, log(x), 1/x, 1/x^2]

%Perform feature scaling
% If we don't do this, we find that gradient descent overshoots
% at first iteration itself, as some features are lot more
% impacting than others
Xf = (Xf - mean(Xf)) ./ (max(Xf) - min(Xf));

% Prepare for linearCost
%    X1 = vector of 1s added to X
X1 = [ones(size(Xf,1),1), Xf];
theta = rand(size(X1,2), 1);
alpha = 1.0;
lambda = 0;

%Plot initial X
plot(X,y);
fprintf('Initial plot');
pause;

iterations = 100000;
step_size = iterations / 10; 
i = 0;
theta_prev = theta;
J_prev = 0;
while i < iterations,
	[J, grad] = linearCost(theta, X1, y, alpha, lambda);
	%fprintf('Running iteration %d, cost: %f\n', i, J);
	if i > 0,
		if J > J_prev,
			% We have started to overshoot, stop here
			break;
		end
	end

	J_prev = J;
	theta_prev = theta;

	theta = theta - grad;
	i = i + 1;

	if mod(i, step_size) == 0,
		clf;
		plot(X,y);
		hold on;
		h = X1 * theta;
		plot(X, h);
		hold off;
		pause(0.5);
	end
		
end

fprintf('final theta at iteration %d: \n', i);
theta

clf;
plot(X,y);
hold on;
h = X1 * theta;
plot(X,h);
hold off;
