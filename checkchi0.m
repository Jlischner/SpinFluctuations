more off;
addpath("/auto/jlischner/chiM/octave_code");

%# read input.m file:
input;
setup;

%# read stuff (use epsinvomega code):
nmtxdat = load("nmtx.dat");
%# factor of two because spin-resolved chi (see Rousseau/Bergara)
chi0 = load("fulleps.dat")/2;
vcdat  = load("gvecs.dat");
n = load("cd.dat");

Nq = length(nmtxdat);
chi0 = chi0(:,1) + I*chi0(:,2);
chi02 = zeros(size(chi0));
Nfreq = length(chi0)/sum(nmtxdat.^2);

%# calculate Ixc:
getIxc;

chi0h = zeros(Nfreq,Nq);
chih  = zeros(Nfreq,Nq);
eig_chi  = zeros(Nfreq,Nq);
eig_chi0 = zeros(Nfreq,Nq);
eig_inv  = zeros(Nfreq,Nq);
Ws = zeros(size(chi0));

for iq = 1:Nq;

  printf("doing %d of %d qpoints \n",iq,Nq);
  nmtx  = nmtxdat(iq);
  ind1 = sum(nmtxdat(1:iq-1));
  ind2 = sum(nmtxdat(1:iq)) ;
  gvecs = vcdat( ind1+1 : ind2 , : );  
  indx = FFTbox(gvecs,S);

  ind1 = sum(nmtxdat(1:iq-1).^2);
  ind2 = sum(nmtxdat(1:iq).^2) ;
  chi0w = chi0( ind1*Nfreq+1 : ind2*Nfreq ) ;
  %# fortran/octave convention: row index runs fastest
  %# freqs written out fastest => should be rows
  chi0w = reshape(chi0w,Nfreq,nmtx^2);
  Ww = zeros(size(chi0w));
  chi0w2 = zeros(size(chi0w));

  for ifreq = 1:Nfreq;

    printf("doing %d of %d freqs \n",ifreq,Nfreq);
    chi0f = reshape( chi0w(ifreq,:),nmtx,nmtx);

    epsmat = zeros(nmtx,nmtx);
    for nrow = 1:nmtx;
      
      %# transform each row of chi0 to real space
      %# multiply by Ixc and transform back
      #row_chi0 = zeros(lenS,1);
      #row_chi0(indx) = chi0f(nrow,:);
      #row_chi0   = cJ( row_chi0);
      #epsmat(nrow,:) = cI( row_chi0 .* Ixc)(indx);
    endfor;
    
    epsmat = eye(nmtx) - epsmat;
    epsinv = inv(epsmat);
    %# chi = chi0 + chi0 * Ixc * chi = epsInv * chi0
    chi = epsinv * chi0f;
 
    %#------------------------------------------------------
    %# calculate screened interaction: W = Ixc * epsInv - Ixc
    %# note: epsilon = 1 - chi0 * Ixc
    Wmat = zeros(nmtx,nmtx);
    for ncol = 1:nmtx;
      
      %# transform each col of epsinv to real space
      %# multiply by Ixc and transform back
      #col_epsI = zeros(lenS,1);
      #col_epsI(indx) = epsinv(:,ncol);
      #col_epsI   = cI( col_epsI);
      #Wmat(:,ncol) = cJ( col_epsI .* Ixc)(indx);
    endfor;
    
    chi0w2(ifreq,:) = reshape(chi0f,nmtx^2,1);
    Ww(ifreq,:) = reshape(Wmat,nmtx^2,1); 
    %#------------------------------------------------------

    #eig_inv(ifreq,iq)  = eigs(epsinv,1);
    #eig_chi(ifreq,iq)  = eigs(chi-chi',1);
    #eig_chi0(ifreq,iq) = eigs(chi0f-chi0f',1); 
    chi0h(ifreq,iq) = chi0f(1,1);
    chih(ifreq,iq)  = chi(1,1);
  endfor; %# freq loop
chi02(ind1*Nfreq+1 : ind2*Nfreq) = reshape(chi0w2,Nfreq*nmtx^2,1);
Ws(ind1*Nfreq+1 : ind2*Nfreq) = reshape(Ww,Nfreq*nmtx^2,1);
endfor; %# iq loop

Ws *= 2; %# Rousseau factor

save output chi0h chih;
 
more on;