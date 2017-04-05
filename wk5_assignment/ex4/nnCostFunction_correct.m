function [J grad] = nnCostFunction_correct(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% nn_params => Total number of \thetas (unrolled)
% input_layer_size => Number of feature variables in input (count of x1,x2,x3... xn)
% hidden_layer_size => Number of nodes in hidden layer
% num_labels => number of output classes (number of nodes in output layer)
% X => Input (without the bias node)
% y => Output
% lambda => Regularization parameter

% Based on above
% Theta1 = s(j+1) x s(j) + 1    => hidden_layer_size x (input_layer_size + 1)
% Theta2 = s(j+1) x s(j) + 1    => num_labels x (hidden_layer_size + 1)
% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

% ====================== YOUR CODE HERE ======================
% Instructions: You should complete the code by working through the
%               following parts.
%
% Part 1: Feedforward the neural network and return the cost in the
%         variable J. After implementing Part 1, you can verify that your
%         cost function computation is correct by verifying the cost
%         computed in ex4.m

% To compute the cost, first need to find h(x)
% For that, first compute hidden layer

% Add bias layers to X
X = [ones(m,1), X];

% a = X * (Theta1)'     (m x n+1) * (n+1 x A)  where A = num nodes in hidden layer
a = X * (Theta1)';   % m x A
a = sigmoid(a);

% Add bias layers to A
a = [ones(m,1), a];

% a2 = a * (Theta2)'   (m x A+1) * (A+1 x num_labels) 
a2 = a * (Theta2)';    % m x num_labels
a2 = sigmoid(a2);

% y is currently an mx1 matrix with values equals to labels
% Convert it into m x num_labels matrix with each column having the corresponding 1
Y = zeros(m, num_labels);
for i = 1:m,
	Y(i,y(i)) = 1;
end

% Now compute cost
J = Y .* log(a2) + (1 - Y) .* log(1 - a2);
J = sum(sum(J));
J = (-1 / m) * J;

% Regularize the cost
% Sum over Theta*Theta for all excluding Theta0s
temp = Theta1(:, 2:end);
temp = temp .* temp;
temp2 = Theta2(:, 2:end);
temp2 = temp2 .* temp2;
reg = (lambda / (2 * m)) * (sum(sum(temp)) + sum(sum(temp2)));
J = J + reg;



%
% Part 2: Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. You should return the partial derivatives of
%         the cost function with respect to Theta1 and Theta2 in Theta1_grad and
%         Theta2_grad, respectively. After implementing Part 2, you can check
%         that your implementation is correct by running checkNNGradients
%
%         Note: The vector y passed into the function is a vector of labels
%               containing values from 1..K. You need to map this vector into a 
%               binary vector of 1's and 0's to be used with the neural network
%               cost function.
%
%         Hint: We recommend implementing backpropagation using a for-loop
%               over the training examples if you are implementing it for the 
%               first time.
for i = 1:m,
	% Extract i'th training set (bias is already added to X)
	x = X(i,:);   %1 x n+1

	%Theta1 is Axn
	a1 = sigmoid(x * Theta1');   % 1 x A
	a1 = [1,a1];	    % becomes 1 x A+1

	%Theta2 is num_labels x (A+1)
	a2 = sigmoid(a1 * Theta2');  % 1 x num_labels 

	% Convert y to a vector of zeros having 1 in the label index
	Y = zeros(1,num_labels);   % 1 x num_labels
	Y(y(i)) = 1;

	% delta2 => error in last layer
	d2 = a2 - Y;     % 1 x num_labels

	% delta1 => error in hidden layer
	% => deltaL * Theta2  .* a1 .* (1 - a1)   => (1 x num_labels) x (num_labels x A+1) .* (1 x A+1) .* (1 - (1 x A+1))
	d1 = (d2 * Theta2) .* a1 .* (1 - a1);    % 1 x A+1

	%Skip bias-error
	d1 = d1(1,2:end);   % 1 x A

	% Theta2_grad += d2' * a1  => (num_labels x 1) * (1 x A+1)
	Theta2_grad = Theta2_grad + (d2' * a1);

	% Theta1_grad += d1' * x  => (A x 1) * (1 x n+1)
	Theta1_grad = Theta1_grad + (d1' * x);
	if i < 3,
	%	fprintf('activations in correct: %dx%d\n', size(a1,1), size(a1,2));
	%	fprintf('activations in correct: %dx%d\n', size(x,1), size(x,2));
	end
	
end

Theta1_grad = (1/m) * Theta1_grad;
Theta2_grad = (1/m) * Theta2_grad;




%
% Part 3: Implement regularization with the cost function and gradients.
%
%         Hint: You can implement this around the code for
%               backpropagation. That is, you can compute the gradients for
%               the regularization separately and then add them to Theta1_grad
%               and Theta2_grad from Part 2.
%
% Save theta1_row0
temp = Theta1 * (lambda/m);
temp(:,1) = zeros(size(Theta1,1),1);

Theta1_grad = Theta1_grad + temp;

temp = Theta2 * (lambda/m);
temp(:,1) = zeros(size(Theta2,1),1);
Theta2_grad = Theta2_grad + temp;


% -------------------------------------------------------------

% =========================================================================

% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];


end
