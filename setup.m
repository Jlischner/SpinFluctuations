%S=[ 25; 25; 25]
%R=[ 0.500000    0.5000      0.5000;
%   -0.50000     0.5000      0.5000;
%   -0.50000    -0.5000      0.5000]
%
%R *= 8.1037;

Vcell=abs(det(R));
dV=abs(det(R))/prod(S);
lenS = prod(S);

ms=[0:prod(S)-1];
ms=ms.';
m1=rem(ms,S(1));
m2=rem(floor(ms/S(1)),S(2));
m3=rem(floor(ms/(S(1)*S(2))),S(3));

n1=m1-(m1>S(1)/2)*S(1);
n2=m2-(m2>S(2)/2)*S(2);
n3=m3-(m3>S(3)/2)*S(3);
N=[n1, n2, n3];
M=[m1, m2, m3] * inverse(diag(S));

r=M*R.';
tlc = norm(r(1,:)-r(2,:)) *2;

#Ms = M;
#for ii=1:3;
#  vec = M(:,ii);
#  vec(vec>0.5) -= 1;
#  Ms(:,ii) = vec;
#endfor;
#rs = Ms*R.';
#absrs = sqrt( sum( rs.^2,2) );

G=2.*pi*N*inverse(R);
G2=sum(G.^2,2);

global gbl_S = S; 
global gbl_N = N;