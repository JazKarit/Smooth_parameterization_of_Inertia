function r = sample_rand(lower_bound, upper_bound, size_vec)
% SAMPLE_RAND Samples random values between lower_bound and upper_bound.
%
%   r = SAMPLE_RAND(lower_bound, upper_bound)
%       returns a single random value between the bounds.
%
%   r = SAMPLE_RAND(lower_bound, upper_bound, size_vec)
%       returns an array of random values of size `size_vec` between the bounds.

    if nargin < 3
        size_vec = [1, 1]; % default to single value
    end
    
    r = lower_bound + (upper_bound - lower_bound) .* rand(size_vec);
end