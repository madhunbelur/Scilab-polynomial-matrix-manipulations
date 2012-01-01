// This Scilab program was written by Madhu N. Belur (IIT Bombay) in
// year 2008 in the context of an alternative method to compute the
// H-infinity norm of a transfer matrix.
// Comments/criticisms about the program -> Madhu
// Madhu Belur's website: http://www.ee.iitb.ac.in/~belur
// This is a free program, you can do whatever you want.

// program that takes a polynomial matrix P and constructs 
// a *linearization* E and A such that the spectrum of * A and -E * is exactly
// zeros of P

// from D. Steven Mackey's thesis (found on the net) suggested by Bibhas Adhikari.
// The program below is for *symmetric* coefficients of the polynomial matrix.

function [A,E]=make_pencil_for_symm(P)

  sizeProw_col=size(P);
  sizeP=sizeProw_col(1);
  degP=max(degree(P));

  Pcoeff_rev=coeff(P,[degP:-1:0]);
  Pcoeff_revA=coeff(P,[degP-1:-1:0]);
  Pcoeff_revE=coeff(P,[degP-2:-1:0]);

  A=constr_hankel([Pcoeff_revA zeros(sizeP,(degP-1)*sizeP)]);
  Esubmtx=constr_hankel([Pcoeff_revE zeros(sizeP,(degP-2)*sizeP)]);

  E=[coeff(P,degP) zeros(sizeP,(degP-1)*sizeP);zeros(degP*sizeP-sizeP,sizeP) -Esubmtx]

endfunction

// --------------- reproducing constr_hankel.sci for completeness sake
// test examples later below

// This Scilab program was written by Madhu N. Belur (IIT Bombay) in
// year 2008 in the context of an alternative method to compute the
// H-infinity norm of a transfer matrix.
// Comments/criticisms about the program -> Madhu
// Madhu Belur's website: http://www.ee.iitb.ac.in/~belur
// This is a free program, you can do whatever you want.

// to construct the block hankel matrix

// this function takes a WIDE matrix A containing various square
// blocks (of same size) kept in block row. Row dimension (=n) is
// taken as size of block. Number of columns should be divisible by n.
// Number of *blocks* (n_blocks:=col_size/n, say) should be odd.
// If q is the size of the output hankel matrix, then 2*q/n -1 blocks
// are expected as input inside A. 

// The row-dimension should ideally divide number of columns.
// The number of blocks should be 
// of columns of A. Suppose A has (2*m-1) blocks as rows, 
// and each block a square matrix of size n by n

function [hankel_mat,m]=constr_hankel(A)

  row_col=size(A);
  n=row_col(1);
  m=(row_col(2)/n+1)/2;
  hankel_mat=A(:,1:m*n);
  for i=2:m,
    hankel_mat = [hankel_mat;A(:,(i-1)*n+1:(m+i-1)*n)];
  end 
endfunction

// test examples now
// I vaguely remember some regularity conditions necessary for the linearization
// to work (read Mackey's thesis)

s = %s;
P = [0 1;1 0]+ s*[2 3;3 0]+s^2*[4 -9;-9 8];
[A,E] = make_pencil_for_symm(P)
list1 = roots(det(P))
list2 = spec(A,-E)
