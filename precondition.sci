// This Scilab program was written by Madhu N. Belur (IIT Bombay) in
// year 2008 in the context of an alternative method to compute the
// H-infinity norm of a transfer matrix.
// Comments/criticisms about the program -> Madhu
// Madhu Belur's website: http://www.ee.iitb.ac.in/~belur
// This is a free program, you can do whatever you want.

// function to precondition a matrix to such that first and
// last  element have same absolute value. (for polynomial
// matrices, the scaling is on the indeterminate

// function o_mtx=precondition(in_mtx)
function o_mtx=precondition(in_mtx)
  szinp = size(in_mtx); rows=szinp(1); cols=szinp(2);
  lambda = (abs(in_mtx(rows,1))/abs(in_mtx(rows,cols))/((cols-1)^(0.10)))^(1/(cols-1));
  scaling_vect_mat = ones(rows,1)*lambda^[0:cols-1]/abs(in_mtx(rows,1));
  o_mtx = in_mtx.*scaling_vect_mat
endfunction


