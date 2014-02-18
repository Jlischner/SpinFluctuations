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
  for ii = 1:size(atoms)(1);
    absrs = sqrt( sum( ( r - ones(lenS,1)*atoms(ii,:) ).^2,2 ) );
    absrs(absrs<tlc) = tlc;
    nci = interp1(nc(:,1),nc(:,2),absrs) ./(4*pi*absrs.^2);
    n += nci;
  endfor;
endif;

Nq = length(nmtxdat);
chi0 = chi0(:,1) + I*chi0(:,2);
Nfreq = length(chi0)/sum(nmtxdat.^2);

%# calculate Ixc:
getIxc;

chi0h = zeros(Nfreq,Nq);
chih  = zeros(Nfreq,Nq);
eigchi0 = zeros(Nfreq,Nq);
eigchi  = zeros(Nfreq,Nq);
eigeps  = zeros(Nfreq,Nq);
Ws = zeros(size(chi0));

for iq = 1:Nq;

  printf("doing %d of %d qpoints \n",iq,Nq);
  nmtx  = nmtxdat(iq);
  ind1 = sum(nmtxdat(1:iq-1));
  ind2 = sum(nmtxdat(1:iq)) ;
  gvecs = vcdat( ind1+1 : ind2 , : );  %# get gvecs for current qpoint
  indx = FFTbox(gvecs,S);

  ind1 = sum(nmtxdat(1:iq-1).^2);
  ind2 = sum(nmtxdat(1:iq).^2) ;
  chi0w = chi0( ind1*Nfreq+1 : ind2*Nfreq ) ; %# get chi0 for current freq+qpoint

  %# fortran/octave convention for converting vector-> matrix: row index runs fastest
  %# e.g. [a11 a12; a21 a22]->[a11 a21 a12 a22] or [1 2 3 4]->[1 3;2 4]
  %# freqs written out fastest => should be rows
  chi0w = reshape(chi0w,Nfreq,nmtx^2);
  Ww = zeros(size(chi0w));

  for ifreq = 1:Nfreq;

    printf("doing %d of %d freqs \n",ifreq,Nfreq);
    chi0f = reshape( chi0w(ifreq,:),nmtx,nmtx); %# Re chi0 symmetric, Im chi0 changes sign

    chi0I = zeros(nmtx,nmtx);
    for nrow = 1:nmtx;

      %# transform each row of chi0 to real space
      %# multiply by Ixc and transform back
      row_chi0 = zeros(lenS,1);
      row_chi0(indx) = chi0f(nrow,:);
      row_chi0   = cJ( row_chi0);
      chi0I(nrow,:) = cI( row_chi0 .* Ixc)(indx);
    endfor;
    
    epsmat = eye(nmtx) - chi0I;
    eigeps(ifreq,iq) = eigs(epsmat,1,'sm'); 
    epsinv = inv(epsmat);

    %# chi = chi0 + chi0 * Ixc * chi = epsInv * chi0 where eps=1 - chi0 * Ixc
    chi  = epsinv * chi0f;
    chiI = epsinv * chi0I;
    %#------------------------------------------------------
    %# calculate screened interaction: W = Ixc * chi * Ixc = Ixc * (epsInv * chi0 * Ixc)
    Wmat = zeros(nmtx,nmtx);
    for ncol = 1:nmtx;

      %# transform each col of epsinv to real space
      %# multiply by Ixc and transform back
      col_epsI = zeros(lenS,1);
      #col_epsI(indx) = epsinv(:,ncol);
      col_epsI(indx) = chiI(:,ncol);
      col_epsI   = cI( col_epsI);
      Wmat(:,ncol) = cJ( col_epsI .* Ixc)(indx);
    endfor;
    
    Ww(ifreq,:) = reshape(Wmat,nmtx^2,1); 
    %#------------------------------------------------------

    chi0h(ifreq,iq) = chi0f(1,1);
    chih(ifreq,iq)  = chi(1,1);
    eigchi0(ifreq,iq) = eigs(chi0f,1,'lm'); 
    eigchi(ifreq,iq)  = eigs(chi,1,'lm'); 
  endfor; %# freq loop
Ws(ind1*Nfreq+1 : ind2*Nfreq) = reshape(Ww,Nfreq*nmtx^2,1);
endfor; %# iq loop

Ws *= 3/2; %# factor 3: Overhauser, factor 1/2: double chi, half each Ixc

%# calculate freq grid:
ws = [];
tmpfreq = init_frequency;
ifreqcounter = 0;
freqstep = delta_frequency;
while ( tmpfreq <= frequency_high_cutoff);

  ws = [ws tmpfreq];
  ifreqcounter += 1;
  if(tmpfreq < frequency_low_cutoff);
    tmpfreq += delta_frequency;
  else;
    freqstep += delta_frequency_step;
    tmpfreq  += freqstep;
  endif;
endwhile;

%# write output for file:
save output ws chi0h chih eigchi0 eigchi eigeps;
fid = fopen("Wdat","w");
fprintf(fid,"%f %f \n",[real(Ws) imag(Ws)]');
fclose(fid); 

more on;