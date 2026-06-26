function feigenbaum
% feigenbaum.m
% Jeremie Korta
% December 1, 2003
% The following program will compute the Feigenbaum ratio using parameter
% values r(n) at superstable fixed points of f^(2^n).  As a check, the
% program can be used to estimate the same ratio using period doubling
% cascade from 3-cycle, which is superstable at r = 1.7549... taking the
% unimodular map to be r-x^2.

clear n r R c;
global n q c

max=10;    
q = 2; % if you want to find the Feigenbaum ratio for higher order maxima
% change this parameter (the 3-cycle will no longer be stable at the same r!)

% Known superstable r for 1 and 2-cycles for any q.
c = 0; % signifies we're searching from 2-cycle
r(1)=0;   % for 1-cycle to be superstable, r = 0
r(2)=1;  % for 2-cycle, r = 1.

% We can also calculate the Feigenbaum ratio at the 3-cycle
% period-doubling cascade (for q=2).  Hopefully, it's still 4.6692... !
%c = 1; % signifies we're searching from 3-cycle
%r(1) = 1.7549;  % superstable 3-cycle occurs at r = 1.7549...
          
for i=2:max
    n=i-1;
    del=10^(-n-c); % smaller for 3-cycle window
    next=r(i-1)+del;
    prev_sign=sign(func(next));
    while (sign(func(next))==prev_sign)
       next=next+del;
    end

    r(i)=fzero(@func, [r(i-1)+10^(-n-1-c), next]);

end

for j=1:(max-2)
    R(j)=(r(j)-r(j+1))/(r(j+1)-r(j+2));
end

fprintf('\n Best estimate of Feigenbaum ratio = %1.6f', R(max-2))




% File: func.m
%
% func returns the f^(2^n) map of 0 where f = r-x^(2q) q=1,2,3... if
% c = 0, or the f^(3*2^n) map of 0 if c = 1.  The zero of this function
% occurs at the superstable value r(n).

function x=func(r)
   
global n q c
x=0;
for i=1:((c > 0)*2+1)*2^n
    x = r - x^q;
end