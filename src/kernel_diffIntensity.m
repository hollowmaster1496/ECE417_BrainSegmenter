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
## @deftypefn {} {@var{retval} =} kernel_diffIntensity (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Haseeb Syed <s2haseeb@edu.uwaterloo.ca>
## Created: 2019-04-22

% Computes g_jk (ie. | x_j - x_k | for a neighborhood of pixels
function [X] = kernel_diffIntensity (neighborhood)
  center = neighborhood(5);
  sum = 0;
  for i = 1:size(neighborhood)
    sum = sum + abs(neighborhood(i) - center);
  end
  
  X = sum;
endfunction
