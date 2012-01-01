// This Scilab program was written by Madhu N. Belur (IIT Bombay) in
// year 2008 in the context of an alternative method to compute the
// H-infinity norm of a transfer matrix.
// Comments/criticisms about the program -> Madhu
// Madhu Belur's website: http://www.ee.iitb.ac.in/~belur
// This is a free program, you can do whatever you want.

// to construct the bezoutian of two polynomials.
// Note: this program can take coefficients that are themselves
// polynomials (in another variable, say z).
// Hence coefficients have to be provided as a VECTOR 

// The output is a symmetric MATRIX, and if coefficients were
// polynomials (in one indeterminate say z), then output will be polynomial
// matrix in indeterminate z

// Output is understood as follows. If B is the output matrix, then
// [1 x x^2 ...] * B *[1 y y^2 ..]' is equal to ( a(x)b(y)-b(x)a(y) ) / (x-y)
// where *coefficients* of polynomials a and b are given as input to this
// program. (If coefficients are always coefficients, then it sure
// is possible to modify this program to take the polynomial (instead of 
// coefficients. -Madhu Belur.)


function bezou_mat=bezoutian(a1,b1)
  // equalize sizes of a and b coeffient vector
  if ~(length(a1)==length(b1)) then
     a=[a1 zeros(1,max(0,length(b1)-length(a1)))]
     b=[b1 zeros(1,max(0,length(a1)-length(b1)))]
  else
     a=a1; b=b1;
  end
  
  cfp=a'*b-b'*a
  
  [rowbz,colbz]=size(cfp);
  
  bezo=[-cfp(1,2:colbz) 0];
  for i=2:rowbz,
      newrow=[bezo(i-1,:)]-cfp(i,1:colbz)
      bezo=[bezo;[newrow(2:colbz) 0]]
  end 
  bezou_mat=bezo(1:rowbz-1,1:colbz-1)
endfunction
  

// first example
pa=poly([2 4 2 -1],'s');
pb=poly([-2 4 1],'s'); // pa and pb have 4 as their common root
acoeff=coeff(pa);
bcoeff=coeff(pb);
B=bezoutian(acoeff,bcoeff)
rank(B) // rank expected to fall by 1 (since only 4 is a common root)
B*[1 4 4^2 4^3]'


// second example (to find parameter z when two polynomials are noncoprime)
z=%z;
acoeff=[5 5] // root is -1
bcoeff=[z 2 1] // root is -1 for SOME z. Now, which z?
Bz=bezoutian(a1,b1)
det(Bz) // roots of det(Bz) are values of z that make Bz singular, hence common root

