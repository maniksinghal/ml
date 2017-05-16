function [C, sigma] = dataset3Params(X, y, Xval, yval)
%DATASET3PARAMS returns your choice of C and sigma for Part 3 of the exercise
%where you select the optimal (C, sigma) learning parameters to use for SVM
%with RBF kernel
%   [C, sigma] = DATASET3PARAMS(X, y, Xval, yval) returns your choice of C and 
%   sigma. You should complete this function to return the optimal C and 
%   sigma based on a cross-validation set.
%

% You need to return the following variables correctly.
C = 1;
sigma = 0.3;

% ====================== YOUR CODE HERE ======================
% Instructions: Fill in this function to return the optimal C and sigma
%               learning parameters found using the cross validation set.
%               You can use svmPredict to predict the labels on the cross
%               validation set. For example, 
%                   predictions = svmPredict(model, Xval);
%               will return the predictions on the cross validation set.
%
%  Note: You can compute the prediction error using 
%        mean(double(predictions ~= yval))
%

chosen_C = 0.01;
chosen_sigma = 0.01;
chosen_cost = 9999999999;

C = 0.01;
while C < 30
	sigma = 0.01;
	while sigma < 30
		model= svmTrain(X, y, C, @(x1, x2) gaussianKernel(x1, x2, sigma)); 
		pred = svmPredict(model, Xval);
		cost = sum((pred - yval) .* (pred - yval));
		if cost < chosen_cost,
			chosen_C = C;
			chosen_sigma = sigma;
			chosen_cost = cost;
		endif
		sigma = sigma * 3;
	endwhile
	C = C * 3;
endwhile

C = chosen_C;
sigma = chosen_sigma;





% =========================================================================

end
