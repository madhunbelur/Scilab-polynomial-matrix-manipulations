// ( [ 
// function to compute the  matrix corresponding to the
// two variable polynomial matrix
function [big_mat,max_deg_w_plus1]=two_var_mat(M_0,M_1)

 M_1_w=(horner(M_1',-s*%i))*(horner(M_1,s*%i));
 M_0_w=-(horner(M_0',-s*%i))*(horner(M_0,s*%i));
 
 mat_m0=coeff(real(det(M_0_w)));
 mat_mf=coeff(real(det(M_1_w))); // leading coefficient (final depending on degree)
 
 deg_mf_plus1=length(mat_mf);
 deg_m0_plus1=length(mat_m0);
 max_deg_w_plus1=max([deg_m0_plus1 deg_mf_plus1]);
 
 mat_m1=[]; mat_m2=[];
 row_col_M1=size(M_1)
 if row_col_M1(1)==2, 
    m_1=real(det([M_0_w(:,1) M_1_w(:,2)])+det([M_1_w(:,1) M_0_w(:,2)]));
    mat_m1=coeff(m_1);
    deg_m1_plus1=length(mat_m1);
    max_deg_w_plus1=max([max_deg_w_plus1 deg_m1_plus1]);
    mat_m1=[mat_m1 zeros(1,max_deg_w_plus1-length(mat_m1))];
   elseif row_col_M1(1)==3, 
     m_1=real(det([M_0_w(:,1:2) M_1_w(:,3)])+det([M_1_w(:,1) M_0_w(:,2:3)])+ ...
         det([M_0_w(:,1) M_1_w(:,2) M_0_w(:,3)]));
     m_2=real(det([M_1_w(:,1:2) M_0_w(:,3)])+det([M_0_w(:,1) M_1_w(:,2:3)])+ ...
         det([M_1_w(:,1) M_0_w(:,2) M_1_w(:,3)]));
     mat_m1=coeff(m_1); mat_m2=coeff(m_2);
     deg_m1_plus1=length(mat_m1); deg_m2_plus1=length(mat_m2);
     max_deg_w_plus1=max([max_deg_w_plus1 deg_m1_plus1 deg_m2_plus1]);
     mat_m1=[mat_m1 zeros(1,max_deg_w_plus1-length(mat_m1))];
     mat_m2=[mat_m2 zeros(1,max_deg_w_plus1-length(mat_m2))];
 end
 
  big_mat=[mat_m0 zeros(1,max_deg_w_plus1-length(mat_m0)); ....
    mat_m1;....// the zero-padding takes place just after defn of mat_m1
    mat_m2; ...//
    mat_mf];
endfunction


