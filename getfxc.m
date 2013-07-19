function fxc = getfxc(n);

alpha = 3/4 * (3/2/pi)^(2/3);
a =  0.0311;
b = -0.0480;
c =  0.0020;
d = -0.0116;
A = -0.1423;
B =  1.0529;
C =  0.3334;

rs = (3/4/pi)^(1/3) * n.^(-1/3);
indxS = (rs < 1.0);
indxL = (rs >= 1.0);
rsS = rs(indxS);
rsL = rs(indxL);

%# get first derivate: excp = d exc / dn = (d exc / drs) * drs/dn
rsp = (3/4/pi)^(1/3) * (-1/3.) * n.^(-4/3);
%# calculate dexc/drs:
dedrs = alpha./rs.^2;
dedrs(indxS) += a./rsS + c*(log(rsS)+1)+ d;
dedrs(indxL) += -A./(1+B*sqrt(rsL)+C*rsL).^2 .*( B/2./sqrt(rsL)+C );

excp = dedrs .* rsp;

%# second derivative: excpp = d^2 exc / dn^2 = (d^2 exc / drs^2) * (drs/dn)^2 + dexc/drs * (d^2 rs/dn^2)
rspp = (3/4/pi)^(1/3) * (4/9.) .* n.^(-7/3);
%# get second derivative: exc
ddedrs = -2*alpha./rs.^3;
ddedrs(indxS) += -a./rsS.^2 + C./rsS;
ddedrs(indxL) += ( 2*(B/2./sqrt(rsL) +C ).^2 + B/4./rsL.^(3/2) .* (1+B*sqrt(rsL) + C*rsL) ) *A./(1+B*sqrt(rsL)+C*rsL).^3;

excpp = ddedrs .* rsp.^2 + dedrs .* rspp;

fxc = 2*excp + n.*excpp;
endfunction;