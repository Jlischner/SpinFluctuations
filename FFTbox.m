function indx = FFTbox(gkvectors,S);

  global gbl_N;

  Ngk = size( gkvectors)(1);
  lenS = prod(S);
  ivec = gkvectors + ones(Ngk,1)*[1 1 1];

  ivec(:,1) += S(1)*(ivec(:,1)<=0) ;
  ivec(:,2) += S(2)*(ivec(:,2)<=0) ;
  ivec(:,3) += S(3)*(ivec(:,3)<=0) ;
  indx = S(2)*S(1)*( ivec(:,3) -1) + S(1)*( ivec(:,2)-1) + ivec(:,1);

  %# test mapping
  if( sum( abs( gbl_N(indx,:)-gkvectors ),2) > 0.01 );
    printf("problem in gvec mapping - FFTbox.m \n");
    exit(1);
  endif;

endfunction;
