more off;
addpath("/auto/jlischner/chiM/octave_code");
input;
setup;

nmtxdat = load("nmtx.dat");
%# factor of two - see Rousseau paper
epsinv = load("epsinv");
vcdat  = load("vcoul");
n = load("cd.dat");

Nq = length(nmtxdat);
epsinv = epsinv(:,1) + I*epsinv(:,2);

%# calculate Ixc:
getIxc;

chi0h = zeros(Nq,1);
chih  = zeros(Nq,1);

for iq = 1:Nq;

  nmtx  = nmtxdat(iq);
  ind1 = sum(nmtxdat(1:iq-1));
  ind2 = sum(nmtxdat(1:iq)) ;
  gvecs = vcdat( ind1+1 : ind2 , 4:6 );  
  vc = vcdat( ind1+1 : ind2, 7) ;

  indx = FFTbox(gvecs,S);

  ind1 = sum(nmtxdat(1:iq-1).^2);
  ind2 = sum(nmtxdat(1:iq).^2) ;
  epsinvw = epsinv( ind1+1 : ind2 ) ;
  epsinvw = reshape(epsinvw,nmtx,nmtx);
  chi0w = diag(1./vc) * (eye(nmtx) - inv(epsinvw));
  chi0w /= 2;

  epsmat = zeros(nmtx,nmtx);
  for nrow = 1:nmtx;
    
    %# transform each col of chi0 to real space
    %# multiply by Ixc and transform back
    row_chi0 = zeros(lenS,1);
    row_chi0(indx) = chi0w(nrow,:);
    row_chi0   = cJ( row_chi0);
    epsmat(nrow,:) = cI( row_chi0 .* Ixc)(indx);
    
  endfor;
  
  epsmat = eye(nmtx) - epsmat;
  chi = inv(epsmat) * chi0w;

  chi0h(iq) = chi0w(1,1);
  chih(iq)  = chi(1,1);

endfor;  
more on;