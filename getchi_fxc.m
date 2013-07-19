more off;
addpath("/auto/jlischner/chiM/octave_code");

%# read input.m file:
input;
setup;

%# read stuff (use epsinvomega code):
nmtxdat = load("nmtx.dat");
chi0 = load("fulleps.dat");
vcdat  = load("gvecs.dat");
n = load("cd.dat");
vcoul = load("vcoul")(:,7);

Nq = length(nmtxdat);
chi0 = chi0(:,1) + I*chi0(:,2);
Nfreq = length(chi0)/sum(nmtxdat.^2);

%# calculate fxc - multiply by two to get rydberg units:
fxc = 2*getfxc(n);

chi0h = zeros(Nfreq,Nq);
chih  = zeros(Nfreq,Nq);
Ws = zeros(size(chi0));

for iq = 1:Nq;

  printf("doing %d of %d qpoints \n",iq,Nq);
  nmtx  = nmtxdat(iq);
  ind1 = sum(nmtxdat(1:iq-1));
  ind2 = sum(nmtxdat(1:iq)) ;
  gvecs = vcdat( ind1+1 : ind2 , : ); 
  vc   = vcoul( ind1+1 : ind2);
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
      %# multiply by fxc and transform back
      row_chi0 = zeros(lenS,1);
      row_chi0(indx) = chi0f(nrow,:);
      row_chi0   = cJ( row_chi0);
      epsmat(nrow,:) = cI( row_chi0 .* fxc)(indx);
    endfor;
    
    %# chi = chi0 + chi0 * (v+fxc) * chi = epsInv * chi0
    %# epsInv = [ 1 - chi0 * (v+fxc) ]^-1
    epsmat = eye(nmtx) - epsmat - chi0f*diag(vc);
    epsinv = inv(epsmat);
    chi = epsinv * chi0f;
 
    %#------------------------------------------------------
    %# calculate epsilonv inverse: epsinv = 1 + v * chi and write to Ww
    Wmat = eye(nmtx) + diag(vc)*chi;    
    Ww(ifreq,:) = reshape(Wmat,nmtx^2,1); 
    %#------------------------------------------------------

    chi0h(ifreq,iq) = chi0f(1,1);
    chih(ifreq,iq)  = chi(1,1);
  endfor; %# freq loop
Ws(ind1*Nfreq+1 : ind2*Nfreq) = reshape(Ww,Nfreq*nmtx^2,1);
endfor; %# iq loop

save output chi0h chih;
fid = fopen("Wdat","w");
fprintf(fid,"%f %f \n",[real(Ws) imag(Ws)]');
fclose(fid); 

more on;