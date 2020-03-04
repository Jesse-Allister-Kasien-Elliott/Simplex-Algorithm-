


%
% EXAMPLE 1 - bounded
c = [1; 4; 0; 0];
A = [1 2 1 0; 1 1 0 1];
b = [2;4];
BAS = [1;4];
% sol = [0;1;0;3]
%

%{
% EXAMPLE 2 - bounded
c = [0; 1; -2; 0; -3];
A = [1 -2 1 0 2; 0 1 -1 1 3];
b = [2;4];
BAS = [3;4];
% sol = [10;4;0;0;0]
%}

%{
% EXAMPLE 3 - unbounded
c = [1; 0];
A = [1 -1];
b = [0];
BAS = [1];
% sol = [1;1]
%}

%{
% EXAMPLE 4 - unbounded
c = [1; 3; 0; 0];
A = [1 -2 1 0; 1 -1 0 1];
b = [2;2];
BAS = [1;2];
% sol = [1;1;1;0]
%}

%

[stat, sol] = simplex(c,A,b,BAS);
