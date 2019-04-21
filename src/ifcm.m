## Copyright (C) 2019 Haseeb Syed
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} ifcm (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Haseeb Syed <s2haseeb@edu.uwaterloo.ca>
## Created: 2019-04-21

% data: column-stacked input image
% num_clusters: number of segments to group pixels into
% m: exponent to control degree of fuzziness 
function [U, V] = ifcm (data, num_clusters, m)
    % Step 1: Determine num clusters and degree of fuzziness m
    %         Note: These are passed as parameters of function
  
  % Step 2: Execute FCM completely
  U = randn(num_clusters, size(data)); % Initialize Membership matrix U
  U = U./sum(U);              %    with random but normalized values
  
  MAX_ITERATION = 5;
  EPSILON = 0.69;    % Control exit condition for iterative updates to U and V 
  J = zeros(MAX_ITERATION, 1);  % Initialize the cost function
  
  % Execute FCM up to 'MAX_ITERATION' times
  for i = 1:MAX_ITERATION,
    membership_function = U.^m; 
    V = membership_function*data./((ones(size(data, 2), 1)*sum(membership_function'))'); % new center
    
    %dist = distfcm(V, data);       
    for k = 1:size(V, 1),
      dist(k, :) = abs(data - V(k))';  % fill the distance matrix |x - v|
    end
    
    J(i) = sum(sum((dist.^2).*membership_function)); % i'th cost function result
    
    tmp = dist.^(-2/(m-1));      % calculate new U
    U = tmp./(ones(num_clusters, 1)*sum(tmp));
    
    % check termination condition
    if i > 1,
      if abs(J(i) - J(i-1)) < EPSILON, break; end,
    end
  end
  
  
  % Step 3: Utilize final membership of FCM as initial membership of IFCM
    
  % Step 4: At l'th iteration, calculate cluster center v^l using membership u_ij^l
  
  % Step 5: Calculate the improved similarity measurement d^2(x, v)
  
  % Step 6: Compare u_ij and u_ij^(l-1)
  % .  6a) if | u_ij^l - u_ij^(l-1) | < epsilon, then STOP
  % .  6b) otherwise repeat from step 4
  
endfunction
