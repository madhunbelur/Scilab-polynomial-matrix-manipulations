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

