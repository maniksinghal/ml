function [J grad] = nnCostFunction(nn_params, ...
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

% We generalize this further by adding another parameter num_hidden_layers (=1 here)
% assuming hidden_layer_size same for each hidden layer
num_hidden_layers = 1;

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;

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

% Our neural network is
% x0 h10 h20 h30... h<num_hidden_layers>0 y0
% x1 h11 h21 h31....h<num_hidden_layers>1 y1
% ....................
% xk h1k h2k h3k....h<num_hidden_layers>k yk
% ......................
% xm h1n h2n h3n....h<num_hidden_layers>m 
% ...........................
%    h1n h2n h3n....h<num_hidden_layers>n
% where m, n, k may not be equal

params_out = 0;
prev_layer = X;
reg_sum = 0;    % regularized sum
for l = 1:num_hidden_layers+1,
	
	% Find theta dimensions to compute layer l
	n_rows = hidden_layer_size;
	n_cols = hidden_layer_size + 1;
	if l == 1,
		% For first layer, n_cols = input
		n_cols = size(X,2);
	elseif l == num_hidden_layers+1,
		% For last hidden layer, num_rows = num-nodes in output layer 
		n_rows = num_labels;
	end

	% Extract theta to compute layer l
	n_elems = n_rows * n_cols;
	theta = reshape(nn_params(params_out+1:params_out+n_elems), n_rows, n_cols);
	params_out = params_out + n_elems;

	next_layer = sigmoid(prev_layer * theta');
	if l != num_hidden_layers+1,
		% Not computing last layer, add bias for next iteration
		next_layer = [ones(m,1),next_layer];
	end

	prev_layer = next_layer;

	% Compute regularization for cost function
	% Ignore first column as theta0 is not regularized
	theta_reg = theta(:,2:end);
	theta_reg = theta_reg .* theta_reg;
	reg_sum = reg_sum + sum(sum(theta_reg));
end

reg_sum = (lambda / (2*m)) * reg_sum;
output_layer = prev_layer;


% Convert [num_labels] to matrix with [label_value]=1 and rest all num_labels-1 entries as 0s
Y = zeros(m, num_labels);
for i = 1:m,
	Y(i,y(i)) = 1;
end

% Now compute cost  (-1/m) * sum(ylog(hx) + (1-y)log(1-hx))
J = Y .* log(output_layer) + (1 - Y) .* log(1 - output_layer);
J = sum(sum(J));
J = (-1 / m) * J;

% regularized cost
J = J + reg_sum;

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

% Now we start forward propogation followed by back propogation for each training example
theta_grad_stack = zeros(1,0);
for i = 1:m,

	%
	% Forward propogation - 
	%

	% Extract i'th training set (bias is already added to X)
	x = X(i,:);   %1 x n+1

	% Initialize activation stack with x
	% Keep pushing activation values on the stack. They would be needed
	% (in reverse order) during backward propogation.
	activation_stack = [1, size(x,2), x];  % Start by pushing a0 (input layer) on stack
	params_out = 0;
	theta_stack = zeros(1,0);

	prev_layer = x;
	for l = 1:num_hidden_layers+1,
		% Find theta dimensions to compute layer l
		n_rows = hidden_layer_size;
		n_cols = hidden_layer_size + 1;
		if l == 1,
			% For first layer, n_cols = input
			n_cols = size(X,2);
		elseif l == num_hidden_layers+1,
			% Last hidden layer has n_rows = output-layer nodes
			n_rows = num_labels;
		end
	
		% Extract theta from parameters input
		n_elems = n_rows * n_cols;
		theta = reshape(nn_params(params_out+1:params_out+n_elems), n_rows, n_cols);

		%Initialize theta_grad stack to zeros in first iteration 
		if i == 1,
			theta_grad = zeros(size(theta));
			theta_grad_stack = [n_rows, n_cols, theta_grad(:)', theta_grad_stack];
		end

		params_out = params_out + n_elems;
	
		next_layer = sigmoid(prev_layer * theta');
		if l != num_hidden_layers+1,
			% Not computing last layer, add bias for next iteration
			next_layer = [1,next_layer];

			%act_stack = [layer_rows, layer_cols, layer_itself...., act_stack]
			% Push activation for this layer on to the stack, it will be required (in reverse order)
			% during backward propogation
			activation_stack = [size(next_layer,1), size(next_layer,2), next_layer, activation_stack];
		end

		% Keep pushing theta for each layer on the stack
		% It will be required in reverse order during backward propogation
		theta_stack = [size(theta,1), size(theta,2), theta(:)', theta_stack];
	
		prev_layer = next_layer;
	end

	output_layer = prev_layer;
	
	
	% Convert y to a vector of zeros having 1 in the label index
	Y = zeros(1,num_labels);   % 1 x num_labels
	Y(y(i)) = 1;

	% Start back propogation
	delta_prev = 0;
	for l = 1:num_hidden_layers+1,
		if l != 1,

			% Extract theta (in reverse order) from stack
			% to compute delta
			% Push it again at the end (to maintain the reverse order) as it 
			% will be required again during regularization
			t_row = theta_stack(1);
			t_col = theta_stack(2);
			t_el = t_row * t_col;
			theta = reshape(theta_stack(3:2+t_el),t_row,t_col);
			theta_stack = theta_stack(3+t_el:end);
			theta_stack = [theta_stack, t_row, t_col, theta(:)'];

			delta_prev = (delta_prev * theta) .* act_prev .* (1 - act_prev);

			%Ignore the bias error
			delta_prev = delta_prev(1,2:end);

		else,
			% For the last layer (first iteration here), delta computation
			% is different => deltaL = h(x) - y
			delta_prev = output_layer - Y;
		end

		% Pop activations for the layers (in reverse order)
		% They are used to compute delta for 1 to L-1 layers and for computing
		% gradients
		% No need to push them back as they are not needed further
		act_prev_rows = activation_stack(1);
		act_prev_cols = activation_stack(2);
		act_prev_num_el = act_prev_rows * act_prev_cols;
		act_prev = reshape(activation_stack(3:2+act_prev_num_el),act_prev_rows,act_prev_cols);
		activation_stack = activation_stack(3+act_prev_num_el:end);

		% Gradients were initially set to 0 (during first iteration of forward propogation)
		% Keep updating them for every training example and keep pushing them at the end 
		% of the stack for cyclic access through every training iteration.
		tg_row = theta_grad_stack(1);
		tg_col = theta_grad_stack(2);
		theta_num_el = tg_row * tg_col;
		theta_grad = reshape(theta_grad_stack(3:2+theta_num_el),tg_row,tg_col);
		theta_grad_stack = theta_grad_stack(3+theta_num_el:end);


		theta_grad = theta_grad + delta_prev' * act_prev;

		%Store theta_grad again
		theta_grad_stack = [theta_grad_stack, tg_row, tg_col, theta_grad(:)'];
	end

	%deltaL (last layer) did not require thetaL-1
	%deltal (l=1 to L-1) requires thetal
	%But delta0 is not computed (input layer error does not make sense)
	%So theta0 was left out on the stack
	%Pop it and push it back in the end so that we are left with correct
	%order of theta on the stack (thetaL-1 thetaL-2 thetaL-3 .... theta0)
	% This would be needed for regularization below.
	t_row = theta_stack(1);
	t_col = theta_stack(2);
	t_el = t_row * t_col;
	theta = reshape(theta_stack(3:2+t_el),t_row,t_col);
	theta_stack = theta_stack(3+t_el:end);
	theta_stack = [theta_stack, t_row, t_col, theta(:)'];

end

%Fetch theta_grads and combine them into a single grad vector
grad = zeros(0,1);

%
% Part 3: Implement regularization with the cost function and gradients.
%
%         Hint: You can implement this around the code for
%               backpropagation. That is, you can compute the gradients for
%               the regularization separately and then add them to Theta1_grad
%               and Theta2_grad from Part 2.
%
for l = 1:num_hidden_layers+1,

	% grad = (1/m) * grad
	tg_row = theta_grad_stack(1);
	tg_col = theta_grad_stack(2);
	tg_el = tg_row * tg_col;
	theta_grad = reshape(theta_grad_stack(3:2+tg_el),tg_row,tg_col);
	theta_grad = (1/m) * theta_grad;
	theta_grad_stack = theta_grad_stack(3+tg_el:end);

	%Regularize => + (lambda/m)*theta, except for theta0s
	t_row = theta_stack(1);
	t_col = theta_stack(2);
	t_el = t_row * t_col;
	theta = reshape(theta_stack(3:2+t_el),t_row,t_col);
	theta_stack = theta_stack(3+t_el:end);

	theta = theta * (lambda/m);
	theta(:,1) = zeros(size(theta,1),1);
	theta_grad = theta_grad + theta;

	grad = [theta_grad(:) ; grad];

end



% -------------------------------------------------------------

% =========================================================================

% Unroll gradients
% Done above!!


end
