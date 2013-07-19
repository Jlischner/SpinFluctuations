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

if( coreflag == 1);
  nc = load("ncore.dat");
  nci = interp1(nc(:,1),nc(:,2),absrs)./(4*pi*absrs.^2);
  n += nci;
  n(1) = n(2);
endif;

Nq = length(nmtxdat);
chi0 = chi0(:,1) + I*chi0(:,2);
Nfreq = length(chi0)/sum(nmtxdat.^2);

%# calculate Ixc:
getIxc;

chi0h = zeros(Nfreq,Nq);
chih  = zeros(Nfreq,Nq);
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

  for ifreq = 1:Nfreq;

    printf("doing %d of %d freqs \n",ifreq,Nfreq);
    chi0f = reshape( chi0w(ifreq,:),nmtx,nmtx);

    epsmat = zeros(nmtx,nmtx);
    for nrow = 1:nmtx;
      
      %# transform each row of chi0 to real space
      %# multiply by Ixc and transform back
      row_chi0 = zeros(lenS,1);
      row_chi0(indx) = chi0f(nrow,:);
      row_chi0   = cJ( row_chi0);
      epsmat(nrow,:) = cI( row_chi0 .* Ixc)(indx);
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
      col_epsI = zeros(lenS,1);
      col_epsI(indx) = epsinv(:,ncol);
      col_epsI   = cI( col_epsI);
      Wmat(:,ncol) = cJ( col_epsI .* Ixc)(indx);
    endfor;
    
    Ww(ifreq,:) = reshape(Wmat,nmtx^2,1); 
    %#------------------------------------------------------

    chi0h(ifreq,iq) = chi0f(1,1);
    chih(ifreq,iq)  = chi(1,1);
  endfor; %# freq loop
Ws(ind1*Nfreq+1 : ind2*Nfreq) = reshape(Ww,Nfreq*nmtx^2,1);
endfor; %# iq loop

Ws *= 3/2; %# factor 3: Overhauser, factor 1/2: double chi, half each Ixc
save output chi0h chih;
#save Ws Ws
fid = fopen("Wdat","w");
fprintf(fid,"%f %f \n",[real(Ws) imag(Ws)]');
fclose(fid); 

more on;