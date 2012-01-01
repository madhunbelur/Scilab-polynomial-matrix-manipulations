function [NormCandidatesSq,bezou_mat,size_bezou_mat,size_A_E,pg,qg]=better_method(M_0,M_1)
  [big_mat,max_deg_w_plus1]=two_var_mat(M_0,M_1);
  max_deg_w_plus1=(max_deg_w_plus1+1)/2;

  big_mat_not_even=big_mat(:,[1:2:2*max_deg_w_plus1-1]);

  big_mat_not_even=precondition(big_mat_not_even);
  temp_degree_in_g=size(big_mat);
  degree_in_g=temp_degree_in_g(1)-1;

  pg=real(s^[degree_in_g:-1:0]*big_mat_not_even);//this is polynomial in g^(-1)
  gain_at_origin=roots(pg(1));
  qg=pg(2:max_deg_w_plus1).*[1:1:max_deg_w_plus1-1];

  clear big_mat big_mat_not_even degree_in_g
  clear max_deg_w_plus1 temp_degree_in_g vector_monomials_g_inverse
  
  bezou_mat=bezoutian(pg,qg);
// Sylvester=sylvester_mat(pg,qg);
//   clear pg qg 

// detS=determ_patient(Sylvester);

  [A,E]=make_pencil_for_symm(bezou_mat);
  size_bezou_mat=size(bezou_mat,'r');
  size_A_E=size(A,'r');

  candidates_comp=spec(E,-A); // because A was nonsingular, we get
           // inverse of the required eigenvalues. (wrong answer now 26/11/09)

  tolerance_reality=0.001;
good_cndts=[gain_at_origin gsort(candidates_comp((abs(imag(candidates_comp)) < tolerance_reality) & (real(candidates_comp) > 0)))'];

NormCandidatesSq=[gain_at_origin 10000-gsort(10000-good_cndts)];// gain_at_infinity]

endfunction


