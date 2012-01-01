// This Scilab program was written by Madhu N. Belur (IIT Bombay) in
// year 2008 in the context of an alternative method to compute the
// H-infinity norm of a transfer matrix.
// Comments/criticisms about the program -> Madhu
// Madhu Belur's website: http://www.ee.iitb.ac.in/~belur
// This is a free program, you can do whatever you want.


better_method.sci  for overall implementation of H-infinity norm computation
two_var_mat.sci   some two-variable polynomial matrix manipulation,
         (unlikely to be very useful).

bezoutian.sci   Takes two VECTORS (coefficients) and returns bezoutian MATRIX.
constr_hankel.sci Takes a WIDE matrix A (square blocks in a row) and returns
                  a Hankel matrix (A chopped if little too large.)

make_pencil_for_symm.sci Takes a square polynomial matrix P with each
             coefficient a symmetric matrix, and returns
             two large matrices A and E such that zeros(P)=spec(A,-E)
              (Structured linearization)
precondition.sci Some routine preconditioning

sylvester_mat.sci Takes p and q (with deg q = deg p -1) and gives their 
                   Sylvester matrix


