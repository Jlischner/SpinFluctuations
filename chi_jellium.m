function out = chi_jellium(rs);

printf(" \n");
printf("rs   = %f \n",rs);
kF = 1.92/rs; %# from Ashcroft book
chi0 = kF/pi^2;

printf("chi0 = %f a.u. \n",chi0);
printf("note: chi0(q=0,w=0) equals DOS at EF! \n");
printf("note: BGW output in rydberg, multiply by -2! \n \n");

n = 3/4/pi/rs^3;
getIxc;
%# divide by 4 below (factor of 2 for spin-polarized, another for rydberg units)
%# also: chi0 has minus in BGW
chi = chi0/(1+chi0/4*Ixc);
printf("chi = %f a.u. \n",chi);
printf("chi/chi0 = %f \n",chi/chi0);

endfunction;