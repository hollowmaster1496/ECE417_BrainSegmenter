## Copyright (C) 2019 Haseeb
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
## @deftypefn {} {@var{retval} =} kernel_distEuclidean (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Haseeb Syed <s2haseeb@edu.uwaterloo.ca>
## Created: 2019-04-22

% Computes sum of q_jk assuming a square 3x3 neighborhood
function [X] = kernel_distEuclidean (neighborhood)
    % Note: This has been hardcoded for a 3x3 neighborhood
    %       (a_j - a_k)^2 + (b_j - b_k)^2
    X = [sqrt(2), 1, sqrt(2);
               1, 0, 1;      
         sqrt(2), 1, sqrt(2)];
endfunction
