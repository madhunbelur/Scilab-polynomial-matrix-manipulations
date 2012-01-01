// This Scilab program was written by Madhu N. Belur (IIT Bombay) in
// year 2008 in the context of an alternative method to compute the
// H-infinity norm of a transfer matrix.
// Comments/criticisms about the program -> Madhu
// Madhu Belur's website: http://www.ee.iitb.ac.in/~belur
// This is a free program, you can do whatever you want.

// function to build sylvester matrix of polynomial p and another polynomial q
// of degree ONE less (relevant when q is 'derivative' of p)
// Note: the coefficients can THEMSELVES be polynomials (in another variable).
function Sylv=sylvester_mat(p,q);
  degp=length(p)-1;
  // the other degree is one less than this always!
  Sylv=[p zeros(1,degp-2);zeros(1,degp-1) q;q zeros(1,degp-1)];
  for i=1:degp-2;
    Sylv=[Sylv;zeros(1,i) p zeros(1,degp-i-2);zeros(1,i) q zeros(1,degp-i-1)];
  end
endfunction


// degq=length(q)-1;
// Sylv=[p zeros(1,degq-1)];
// for i=1:degq-1;
//   Sylv=[Sylv;zeros(1,i) p zeros(1,degq-i-1)];
// end
// for i=0:degp-1;
//   Sylv=[Sylv;zeros(1,i) q zeros(1,degp-i-1)];
// end
//
// the below one is for general sylvester matrix for polynomials p ang q
// the top case is for the special case when degq is one less than degp

