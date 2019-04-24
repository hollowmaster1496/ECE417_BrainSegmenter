% Copyright (C) 2019 Haseeb
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.

% -*- texinfo -*- 
% @deftypefn {} {@var{retval} =} fcm (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Haseeb Syed <s2haseeb@edu.uwaterloo.ca>
% Created: 2019-04-22

function [U, V] = fcm (data, num_clusters, m)
  % Step 2: Execute FCM completely
  U = ones(num_clusters, size(data, 1)); % Initialize Membership matrix U
  U = U./sum(U);                      %    with random but normalized values
  U = imnoise(U, 'gaussian', 0, 0.7);
  
  MAX_ITERATION = 5;
  EPSILON = 0.02;    % Control exit condition for iterative updates to U and V 
  J = zeros(MAX_ITERATION, 1);  % Initialize the cost function
  
  % Execute FCM up to 'MAX_ITERATION' times
  i = 1;
  while i <= MAX_ITERATION,
    
    % check termination condition
    if i > 1,      
        if abs(U(i) - U(i-1)) < EPSILON, break; end,
    end
    
    i = i+1; % increment iteration counter
    
    membership_function = U.^m; 
    
    % ensure matching dimensions for pointwise division of j elements
    V_numerator = membership_function*data;
    V_denominator = (ones(size(data, 1), 1)*sum(membership_function'))';
    V = V_numerator./V_denominator;
    
    % compute intensity differences between points and centroid |x - v|    
    for k = 1:size(V, 1),
      dist(k, :) = abs(data - V(k))';
    end
    
    J(i) = sum(sum((dist.^2).*membership_function)); % i'th cost function result
    
    % compute new degree of fuzziness and update membership matrix
    fuzziness = dist.^(2/(m-1));
    for k=1:num_clusters,
      U = U + sum(fuzziness);
    end
    
    U = fuzziness./U;
  end
  
end
