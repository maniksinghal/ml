#!/usr/local/bin/octave -qf
clear; close all; clc;

arglist = argv();
addpath("/Users/maniksin/Documents/manik/code/ml/sine_project")
% Parameters for this code
BinaryClassificationManual  = 1;
MultiClassClassification = 0;

data = dlmread(arglist{1}, ' ', 0, 0);

X = data(:,2:end);  % Separate input, first column is output
y = data(:,1);


%-------------------------- Binary classification using logistic regression --------------------------------

%positive = find(y == 1);   % match all inside the circle
%negative = find(y != 1);

% Add more features to X to make it a polynomial equation

% EXPERIMENT: To test high bias, remove these additional features
% and the algorithm will not be able to match

%X = [X, X(:,1) .* X(:,1)];   %[x1, x2, x1^2]
%X = [X, X(:,2) .* X(:,2)];   %[x1, x2, x1^2, x2^2]
%X = [X, log(X)];   % X = [x, x^2, x^3, log(x)]
%X = [X, 1 ./ X];   % X = [x, x^2, x^3, log(x), 1/x]
%X = [X, 1 ./ (X .* X)];   % X = [x, x^2, x^3, log(x), 1/x, 1/x^2]

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
%Xf = X;

m = size(X,1);
n = size(X,2);

%Add bias
X1 = [ones(m,1), Xf];
n = n + 1;

theta = zeros(n,1);

fprintf('BinaryClassificationManual = %d\n', BinaryClassificationManual); 
if BinaryClassificationManual == 1, 

lambda = 0;
alpha = 1;
iterations = 10000;


% Trying fminunc
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);

% Train test subset (80% of training data)
num_rows = int32(m * 0.8)

%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[thetaMin, J] = ...
	fminunc(@(t)(logisticRegression(t, X1(1:num_rows,:), y(1:num_rows), alpha, lambda)), theta, options);

% validation set
X_val = X1(num_rows+1:end,:)
y_val = y(num_rows+1:end,:)

z = X_val * thetaMin;
h = 1 ./ (1 + exp(-z));

positive = find(h >= 0.7);
err = sum (1 .- y_val(positive));

negative = find(h < 0.7);
err = err + sum(y_val(negative) .* 1);

fprintf('Test-set entries: %d, error: %d', size(h,1), err)

%result = [y_val,h]

pause
%positive = find(h >= 0.5);
%negative = find(h < 0.5);

%X_pos = X(positive,:);
%X_neg = X(negative,:);


end


%-------------------------- Multi-class classification using logistic regression --------------------------------
% We set input to be a range of numbers between 0 and 1000...
% X-axis shall be devided into K equal regions
% Output Y shall be the region-number (class) where X lies (y = k)

fprintf('MultiClassClassification= %d\n', MultiClassClassification); 
if MultiClassClassification == 1, 
RANGE = 100;          % X axis range
K = 20;               % 20 classes
X = rand(5000, 2);   % 2 features (x1,x2) to form coordinates, 1000 training examples, all in range 0-1
X = X .* RANGE;    % Bump up range from 0-<RANGE>
K = 4;   % Keeping K small for easy plotting

%x2 does not have any effect
%Classification depends only on value of x1 (X-portion of coordinate x1,x2)
y = floor(X(:,1)) ./ (RANGE/K);  % Gives y in classes 0,1,2,...K-1

if K == 4,
%Plot if number of classes is known and fixed
temp = find(X(:,1) < ((RANGE * 1)/K));
X0 = X(temp,:);

temp = find(X(:,1) < ((RANGE * 2)/K)); 
X1 = X(temp,:);
temp = find(X1(:,1) > ((RANGE * 1)/K));
X1 = X1(temp,:);

temp = find(X(:,1) < ((RANGE * 3)/K)); 
X2 = X(temp,:);
temp = find(X2(:,1) > ((RANGE * 2)/K));
X2 = X2(temp,:);

temp = find(X(:,1) < ((RANGE * 4)/K)); 
X3 = X(temp,:);
temp = find(X3(:,1) > ((RANGE * 3)/K));
X3 = X3(temp,:);

hold off;
fprintf('Initial plot (MultiClass classification)\n');
xlabel('x1');
ylabel('x2');
plot(X0(:,1), X0(:,2), "ob", X1(:,1), X1(:,2), "or", X2(:,1), X2(:,2), "og", X3(:,1), X3(:,2), "om");
title('Initial plot (Multiclass classification)');
printf('Press any key to continue\n');
pause
end


% EXPERIMENT: Not much difference by changing number of iterations from
% 10K to 50K
lambda = 0;
alpha = 0.5;
iterations = 4000;


% Add more features to X to make it a polynomial equation
% Tried different combinations but algorithm did not converge
% very accurately

X = [X, X(:,1) .* X(:,1)];   %[x1, x2, x1^2]
X = [X, X(:,2) .* X(:,2)];   %[x1, x2, x1^2, x2^2]
%X = [X, X(:,1) .* X(:,1) .* X(:,1)];   %[x1, x2, x1^2, x2^2, x1^3]
%X = [X, log(X(:,1))];   % X = [x, x^2, x^3, log(x1)]
%X = [X, 1 ./ X(:,1)];   % X = [x, x^2, x^3, log(x1), 1/x1]
%X = [X, 1 ./ (X .* X)];   % X = [x, x^2, x^3, log(x), 1/x, 1/x^2]


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

theta = zeros(n,K);  %Store theta for each classification in separate column

% Trying fminunc
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', iterations);

%Run fminunc K times and get theta for each
for k = 1:K,

output = (y == k-1);
%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[theta(:,k), J] = ...
	fminunc(@(t)(logisticRegression(t, X1, output, alpha, lambda)), theta(:,k), options);
end

z = X1 * theta;            %mxk => Each column is result for class k
h = 1 ./ (1 + exp(-z));    % also mxk, with probability results for class 1,2,3...K in each row

[value, index] = max(h');
% Now index is an array of indexes of max probability values in each column
index = index';

if K == 4,  % Plot only when dimensions are fixed and known
X0 = find(index == 1);
X0 = X(X0,:);

X1 = find(index == 2);
X1 = X(X1,:);

X2 = find(index == 3);
X2 = X(X2,:);

X3 = find(index == 4);
X3 = X(X3,:);

fprintf('Learnt plot (MultiClass classification)\n');
xlabel('x1');
ylabel('x2');
plot(X0(:,1), X0(:,2), "ob", X1(:,1), X1(:,2), "or", X2(:,1), X2(:,2), "og", X3(:,1), X3(:,2), "om");
title('Learnt plot (Multiclass classification)');
printf('Press any key to continue\n');
pause
end

%CONCLUSION:
% Algorithm learnt a similar pattern as expected, but the pattern was a little shifted.
% Increasing number of iterations or the training data did not help much
% Probably the features need to change for better results

end


close
