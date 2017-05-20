clear; close all; clc;


% Parameters for this code
CompareFmincgWithManualIterations = 1;
ComputeNormalEquation = 1;
BinaryClassificationManual  = 1;

% Use a test vector X
X = [0.01 : 0.01 : 8];   %1xm
X = X';  %Make it mx1

% y = X^2 + 3
%y = (X .* X) + 3;
y = sin(X);
%y = (X + 3) - X;


% EXPERIMENT: Classic case of overfitting due to too many features
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
theta = zeros(size(X1,2), 1);

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
title('Initial plot (Linear regression)');
fprintf('Initial plot\nPress any key to continue');
pause

% Minimize linear cost using manual loop 
iterations = 1000;
f = @linearCost;
[J, thetaMin, J_slope] = minimizeLinearCost(f, X1, theta, y, alpha, lambda, iterations);
fprintf('theta with manual minimizing loop (cost:%f)\n', J); 
thetaMin

%Minimize linear cost using fmincg
options = optimset('GradObj', 'on', 'MaxIter', iterations);
[thetaMin2] = fmincg(@(t)(linearCost(t, X1, y, 1, lambda)), theta, options);
fprintf('theta with fmincg\n'); 
thetaMin2

% Plot comparison of the above two algorithms
h = X1 * thetaMin;
hold on;
plot(X,h);
title('Learnt plot with manual learning');
fprintf('Plotting with manual minimizing\n');
pause
h = X1 * thetaMin2;
fprintf('Plotting with fmincg\n');
plot(X,h);
title('Comparison with fmincg');
pause
hold off;

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
title('Initial plot (Normal equation)');
hold on;
fprintf('Initial plot\nPress any key to continue');
pause;
thetaMin = (inv(X1' * X1) * X1') * y;
fprintf('thetaMin\n');
thetaMin
h = X1 * thetaMin;
plot(X, h, 'r');
hold off;
title('Learnt plot with Normal equation (in red)');
fprintf('Plot with normal equation (In RED)\nPress any key to continue');
pause

% Conclusion -
% Looks good, as number of features are very small
end





%-------------------------- Binary classification using logistic regression --------------------------------
% We set input to be a circular region characterized by two features x1 and x2
% If the input vector lies within the circular region, y = 1, else y = 0

centre_x1 = 10;
centre_x2 = 10;
radius = 5;
X = rand(1000, 2);   % 2 features, 1000 training examples, all in range 0-1
X = X .* 20;    % Bump up range from 0-1 to 0-20


y = ((X(:,1) .- centre_x1) .^ 2)   .+  ((X(:,2) .- centre_x2) .^ 2);
y = y < (radius * radius);
positive = find(y == 1); 
negative = find(y != 1);
X_pos = X(positive,:);
X_neg = X(negative,:);

% Add more features to X to make it a polynomial equation

% EXPERIMENT: To test high bias, remove these additional features
% and the algorithm will not be able to match

X = [X, X(:,1) .* X(:,1)];   %[x1, x2, x1^2]
X = [X, X(:,2) .* X(:,2)];   %[x1, x2, x1^2, x2^2]
X = [X, log(X)];   % X = [x, x^2, x^3, log(x)]
X = [X, 1 ./ X];   % X = [x, x^2, x^3, log(x), 1/x]
X = [X, 1 ./ (X .* X)];   % X = [x, x^2, x^3, log(x), 1/x, 1/x^2]

%Perform feature scaling
% If we don't do this, we find that gradient descent overshoots
% at first iteration itself, as some features are lot more
% impacting than others

% EXPERIMENT: When removed feature scaling, the fmincunc and
% manual algorithms did not converge, as too many features above
% without scaling had big impact.
% However, on reducing the features to just x1,x2,x1^2,x2^2, the
% fmincunc was able to converge.
%Xf = X;

Xf = (X - mean(X)) ./ (max(X) - min(X));

m = size(X,1);
n = size(X,2);

%Add bias
X1 = [ones(m,1), Xf];
n = n + 1;

theta = zeros(n,1);

fprintf('BinaryClassificationManual = %d\n', BinaryClassificationManual); 
if BinaryClassificationManual == 1, 

hold off;
fprintf('Initial plot (Binary classification)\n');
plot(X_pos(:,1), X_pos(:,2), "ob", X_neg(:,1), X_neg(:,2), "or");
xlabel('x1');
ylabel('x2');
title('Initial plot (Binary classification)');
printf('Press any key to continue\n');
pause

lambda = 0;
alpha = 1;
iterations = 10000;


% Trying fminunc
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);

%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[thetaMin, J] = ...
	fminunc(@(t)(logisticRegression(t, X1, y, alpha, lambda)), theta, options);

z = X1 * thetaMin;
h = 1 ./ (1 + exp(-z));

positive = find(h >= 0.5);
negative = find(h < 0.5);
X_pos = X(positive,:);
X_neg = X(negative,:);
fprintf('Learned plot\n');
plot(X_pos(:,1), X_pos(:,2), "ob", X_neg(:,1), X_neg(:,2), "or");
xlabel('x1');
ylabel('x2');
title('Learnt plot using fminunc');
pause

% Now derive hypothesis using manual iterations
alpha = 5;            % Needed higher alpha rate to converge the cost
iterations = 1000;   % Need high number of iterations for converging the cost
f = @logisticRegression;
[J, thetaMin, J_slope] = minimizeLinearCost(f, X1, theta, y, alpha, lambda, iterations);

plot(J_slope);
fprintf('Press any keys please\n'); 
fprintf('J-slope (iterations:%d)\n', size(J_slope,1));
fprintf('Press any key to continue\n');
title('J-slope using manual iterations');
pause

z = X1 * thetaMin;
h = 1 ./ (1 + exp(-z));


positive = find(h >= 0.5);
negative = find(h < 0.5);
X_pos = X(positive,:);
X_neg = X(negative,:);
fprintf('Learned plot\n');
plot(X_pos(:,1), X_pos(:,2), "ob", X_neg(:,1), X_neg(:,2), "or");
xlabel('x1');
ylabel('x2');
title('Learnt plot using manual loop'); 
pause

% Now plot the desicion boundary
% Plot from training examples where z ~ 0
% (h = 0.5 when z = 0)
negative = find(z < -0.2 | z > 0.2);
positive = find(z > -0.2 & z < 0.2);
X_pos = X(positive,:);
X_neg = X(negative,:);
fprintf('Learned plot\n');
plot(X_pos(:,1), X_pos(:,2), "ob", X_neg(:,1), X_neg(:,2), "or");
%plot(X_neg(:,1), X_neg(:,2), "or");
xlabel('x1');
ylabel('x2');
title('Decision boundary');
pause

end


close
