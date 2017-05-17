clear; close all; clc;


% Parameters for this code
CompareFmincgWithManualIterations = 0;
ComputeNormalEquation = 1;

% Use a test vector X
X = [0.01 : 0.01 : 8];   %1xm
X = X';  %Make it mx1

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

% -------------------   Compare fmincg with manual iterations for linear regression ----------------
% Compare the algorithms for manually performing
% descent and using fmincg
% Both the algorithms internally call linearCost to compute
% the cost and gradient for a given X,theta
fprintf('CompareFmincgWithManualIterations = %d\n', CompareFmincgWithManualIterations); 
if CompareFmincgWithManualIterations == 1, 

alpha = 1;
lambda = 0;

%Plot initial X
plot(X,y);
fprintf('Initial plot\nPress any key to continue');
pause

% Minimize linear cost using manual loop 
iterations = 1000;
[J, thetaMin, J_slope] = minimizeLinearCost(X1, theta, y, alpha, lambda, iterations);
fprintf('theta with manual minimizing loop (cost:%f)\n', J); 
thetaMin
fprintf('Press any key to continue\n');
pause

%Minimize linear cost using fmincg
options = optimset('GradObj', 'on', 'MaxIter', iterations);
[thetaMin2] = fmincg(@(t)(linearCost(t, X1, y, 1, lambda)), theta, options);
fprintf('theta with fmincg\n'); 
thetaMin2
fprintf('Press any key to continue\n');
pause

% Plot comparison of the above two algorithms
h = X1 * thetaMin;
hold on;
plot(X,h);
fprintf('Plotting with manual minimizing\n');
fprintf('Press any key to continue\n');
pause
h = X1 * thetaMin2;
plot(X,h);
fprintf('Plotting with fmincg\n');
fprintf('Press any key to continue\n');
pause;

% Conclusion -
% Manual loop uses far more iterations while learning rate alpha is small (0.05), and 
% converges faster with alpha = 1
% fmincg didn't perform well with small alpha, but gave converged to matching results 
% (as with manual loop) in very less time
end


%----------------------------------------- Minimize using Normal equation method -------------------
fprintf('ComputeNormalEquation = %d\n', ComputeNormalEquation); 
if ComputeNormalEquation == 1, 
plot(X,y);
hold on;
fprintf('Initial plot\nPress any key to continue');
pause;
thetaMin = (inv(X1' * X1) * X1') * y;
fprintf('thetaMin\n');
thetaMin
h = X1 * thetaMin;
plot(X, h);
hold off;
fprintf('Plot with normal equation\nPress any key to continue');
pause

% Conclusion -
% Looks good, as number of features are very small
end
