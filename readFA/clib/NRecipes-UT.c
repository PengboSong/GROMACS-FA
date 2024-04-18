/* Numerical Recipes in C
  Some programs
  we rewrite it to start the array from 0 but not 1
  we transfer matrix as array, so the matrix should be dimensioned as n*n and transfered as a[0]
    for an alternative, see NRecipes-UT-t.c
  also change float to double
*/

#include "NRecipes-UT.h"

#define SQR(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int myEigen(double *a, int n, double d[], double e[], double *z)
/* my calculation of the eigenvalues and eigenvectors of a symmetric matrix a[0:n-1][0:n-1]
   on input a[0:n-1][0:n-1] is the symmetry matix to be calculated
   on output, a[][] is replaced by the orthogonal matrix Q to results in a tridiagonal matrix,
              d[0:n-1] returns the eigenvalues.
              z[0:n-1][0:n-1]: the kth column of z returns the normalized eigenvector corresponding to d[k]. 
   e[0:n-1] is a tempary array
 return 1 if success, otherwise return 0
*/
{
	int i, j, k;

	tred2(a, n, d, e); // into tridiagonal matrix
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			z[i * n + j] = a[i * n + j];
	/*			z[i][j]=0.0;
	for(i=0; i<n; i++)
		z[i][i]=1.0;
	*/
	if (tqli(d, e, n, z) == 1)
	{
		// eigenvectors: QX
		/*		for(j=0; j<n; j++)
		{
			for(i=0; i<n; i++)
			{
				e[i]=0.0;
				for(k=0; k<n; k++)
					e[i]+=a[i][k]*z[k][j];
			}
			for(i=0; i<n; i++)
				z[i][j]=e[i];
		}
		*/
		return (1);
	}
	else
		return (0);
}

// chapter2-6, page 70
double pythag(double a, double b)
// Computes (a^2 + b^2 )^1/2 without destructive underflow or overflow.
{
	double absa, absb;

	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb)
		return (absa * sqrt(1.0 + SQR(absb / absa)));
	else
		return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

// c11-3, p.480
int tqli(double d[], double e[], int n, double *z)
// return 1 if succeed, otherwise 0
/* QL algorithm with implicit shifts, to determine the eigenvalues 
and eigenvectors of a real, symmetric, tridiagonal matrix, 
or of a real, symmetric matrix previously reduced by tred2 � 11.2. 
On input, d[0..n-1] contains the diagonal elements of the tridiagonal matrix. 
On output, it returns the eigenvalues. 
The vector e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, 
with e[0] arbitrary. On output e is destroyed. When finding only the 
eigenvalues, several lines may be omitted, as noted in the comments. 
If the eigenvectors of a tridiagonal matrix are de-sired, the matrix 
z[0..n-1][0..n-1] is input as the identity matrix. 
If the eigenvectors of a matrix that has been reduced by tred2 are 
required, then z is input as the matrix output by tred2. 
In either case, the kth column of z returns the normalized eigenvector 
corresponding to d[k]. 
*/
{
	int m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;

	for (i = 1; i <= n - 1; i++)
		e[i - 1] = e[i]; // Convenient to renumber the el-ements of e.
	e[n - 1] = 0.0;

	for (l = 0; l <= n - 1; l++)
	{
		iter = 0;
		do
		{
			for (m = l; m <= n - 2; m++)
			{
				// Look for a single small subdi-agonal element to split the matrix.
				dd = fabs(d[m]) + fabs(d[m + 1]);
				if ((double)(fabs(e[m]) + dd) == dd)
					break;
			}
			if (m != l)
			{
				if (iter++ == 30)
				{
					printf("Too many iterations in tqli.\n");
					return (0);
				}
				g = (d[l + 1] - d[l]) / (2.0 * e[l]); // Form shift.
				r = pythag(g, 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g)); // This is d_m-k_s.
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--)
				{
					// A plane rotation as in the origi-nal QL,
					// followed by Givens rotations to restore
					// tridiag-onal form.
					f = s * e[i];
					b = c * e[i];
					e[i + 1] = (r = pythag(f, g));
					if (r == 0.0)
					{
						// Recover from underflow.
						d[i + 1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					d[i + 1] = g + (p = s * r);
					g = c * r - b;
					/* Next loop can be omitted if eigenvectors not wanted*/
					for (k = 0; k <= n - 1; k++)
					{
						// Form eigenvectors.
						f = z[k * n + i + 1];
						z[k * n + i + 1] = s * z[k * n + i] + c * f;
						z[k * n + i] = c * z[k * n + i] - s * f;
					}
				}
				if (r == 0.0 && i >= l)
					continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}

	return (1);
}

// c11-2, p.474
void tred2(double *a, int n, double d[], double e[])
/* Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1].
On output, a is replaced by the orthogonal matrix Q effecting the transformation.
d[0..n-1] returns the diagonal elements of the tridiagonal matrix, 
and e[0..n-1] the off-diagonal elements, with e[0]=0. Several statements, 
as noted in comments,can be omitted if only eigenvalues are to be found, 
in which case a contains no useful information on output.
Otherwise they are to be included.
*/
// If the eigenvectors of the final tridiagonal matrix are found,
// then the digenvectors of A can be obtained by applying Q
{
	int l, k, j, i;
	double scale, hh, h, g, f;

	for (i = n - 1; i >= 1; i--)
	{
		l = i - 1;
		h = scale = 0.0;
		if (l > 0)
		{
			for (k = 0; k <= l; k++)
				scale += fabs(a[i * n + k]);
			if (scale == 0.0) // Skip transformation.
				e[i] = a[i * n + l];
			else
			{
				for (k = 0; k <= l; k++)
				{
					a[i * n + k] /= scale;			  // Use scaled a's for transformation.
					h += a[i * n + k] * a[i * n + k]; // Form \sigma in h
				}
				f = a[i * n + l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale * g;
				h -= f * g;			  // Now h is equation (11.2.4).
				a[i * n + l] = f - g; // Store u in the ith row of a
				f = 0.0;
				for (j = 0; j <= l; j++)
				{
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j * n + i] = a[i * n + j] / h; // Store u/H in ith column of a
					g = 0.0;						 // Form a element of A � u in g
					for (k = 0; k <= j; k++)
						g += a[j * n + k] * a[i * n + k];
					for (k = j + 1; k <= l; k++)
						g += a[k * n + j] * a[i * n + k];
					e[j] = g / h; // Form element of p in temporarily unused element of e
					f += e[j] * a[i * n + j];
				}
				hh = f / (h + h); // Form K ,equation (11.2.11).
				for (j = 0; j <= l; j++)
				{
					// Form q and store in e overwriting p
					f = a[i * n + j];
					e[j] = g = e[j] - hh * f;
					for (k = 0; k <= j; k++) // Reduce a, equation (11.2.13).
						a[j * n + k] -= (f * e[k] + g * a[i * n + k]);
				}
			}
		}
		else
			e[i] = a[i * n + l];
		d[i] = h;
	}

	/* Next statement can be omitted if eigenvectors not wanted */
	d[0] = 0.0;
	e[0] = 0.0;
	/* Contents of this loop can be omitted if eigenvectors not wanted 
	except for statement d[i]=a[i][i]; */
	for (i = 0; i <= n - 1; i++)
	{
		// Begin accumulation of transformation matrices.
		l = i - 1;
		if (d[i])
		{
			// This block skipped when i=1
			for (j = 0; j <= l; j++)
			{
				g = 0.0;
				for (k = 0; k <= l; k++) // Use u and u/H stored in a to form P � Q
					g += a[i * n + k] * a[k * n + j];
				for (k = 0; k <= l; k++)
					a[k * n + j] -= g * a[k * n + i];
			}
		}
		d[i] = a[i * n + i]; // This statement remains.
		a[i * n + i] = 1.0;	 // Reset row and colum of a to identity matrix for next iteration.
		for (j = 0; j <= l; j++)
			a[j * n + i] = a[i * n + j] = 0.0;
	}

	return;
}

// c2-3, p.46
#define TINY 1.0e-20 // A small number.
int ludcmp(double *a, int n, int *indx, double *d)
/* Given a matrix a[0..n-1][0..n-1],this routine replaces it by the LU decomposition 
of a rowwise permutation of itself. 
a and n are input.
a is output, arranged as in equation (2.3.14) (A=L.U, where L is 
a lower triangular with diagonal equal to 1 and lower part stored in output A, 
U is a upper triangular that stored in upper A as output ); 
indx[0..n-1] is an output vector that records the row permutation effected 
by the partial pivoting;
d is output as +/-1 depending on whether the number of row interchanges was 
even or odd, respectively.
This routine is used in combination with lubksb to solve linear equations 
or invert a matrix. 
*/
// return 1 if succeed; otherwise 0;
{
	int i, imax, j, k;
	double big, dum, sum, temp;
	double *vv; // vv stores the implicit scaling of each row.

	vv = (double *)malloc(n * sizeof(double));
	// vv=vector(1,n);

	*d = 1.0; // No row nterchanges yet.
			  // Loop over rows to get the implicit scaling information.
	for (i = 0; i <= n - 1; i++)
	{
		big = 0.0;
		for (j = 0; j <= n - 1; j++)
			if ((temp = fabs(a[i * n + j])) > big)
				big = temp;
		if (big == 0.0) // No nonzero largest element.
		{
			printf("Singular matrix in routine ludcmp.\n");
			return (0);
			// nrerror("Singular matrix in routine ludcmp");
		}
		vv[i] = 1.0 / big; // Save the scaling.
	}

	// This is the loop over columns of Crout's method.
	for (j = 0; j <= n - 1; j++)
	{
		for (i = 0; i < j; i++) // This s equation (2.3.12) except for i = j.
		{
			sum = a[i * n + j];
			for (k = 0; k < i; k++)
				sum -= a[i * n + k] * a[k * n + j];
			a[i * n + j] = sum;
		}
		big = 0.0; // Initialize for the search for largest pivot element.

		// This is i = j of equation (2.3.12) and i = j +1 ..N of equation (2.3.13).
		for (i = j; i <= n - 1; i++)
		{
			sum = a[i * n + j];
			for (k = 0; k < j; k++)
				sum -= a[i * n + k] * a[k * n + j];
			a[i * n + j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big)
			{ // Is the figure of merit for the pivot better than the best so far?
				big = dum;
				imax = i;
			}
		}
		if (j != imax) // Do we need to nterchange rows?
		{
			for (k = 0; k <= n - 1; k++) // Yes,do so...
			{
				dum = a[imax * n + k];
				a[imax * n + k] = a[j * n + k];
				a[j * n + k] = dum;
			}
			*d = -(*d);		  // ...and change the parity of d
			vv[imax] = vv[j]; // Also nterchange the scale factor.
		}
		indx[j] = imax;
		if (a[j * n + j] == 0.0)
			a[j * n + j] = TINY;
		/* If the pivot element is zero the matrix is singular (at least to the 
precision of the algorithm). For some applications on singular matrices, 
it is desirable to substitute TINY for zero. */
		if (j != n - 1) // Now, finally, divide by the pivot element.
		{
			dum = 1.0 / (a[j * n + j]);
			for (i = j + 1; i <= n - 1; i++)
				a[i * n + j] *= dum;
		}
	} // Go back for the next column in the reduction.

	//	free_vector(vv,1,n);
	free(vv);

	return (1);
}

// c2-3, p.47
void lubksb(double *a, int n, int *indx, double *b, int nb)
/* Solves the set of n linear equat ons A � X = B. B[0..n-1][0..nb-1] are nb vectors 
Here a[0..n-1][0..n-1] is input, not as the matrix A but rather 
as its LU decomposition, determined by the routine ludcmp. 
indx[0..n-1] is input as the permutation vector returned by ludcmp. 
b[0..n-1][0..nb-1] is input as the right-hand side vector B,
and returns with the solution vector X.
a, n, and indx are not modified by this routine and can be left in place 
for successive calls with different right-hand sides b. 
This routine takes into account the possibility that b will begin with
 many zero elements, so it is efficient for use in matrix inversion. 
*/
{
	int i, ii, ip, j;
	double sum;
	int ib;

	for (ib = 0; ib < nb; ib++)
	{
		ii = -1;
		/* When ii is set to a positive value, it will become the index of the  
first nonvanishing element of b. We now do the forward substitution,
equation (2.3.6). The only new wrinkle is to unscramble the permutation 
as we go. */
		for (i = 0; i <= n - 1; i++)
		{
			ip = indx[i];
			sum = b[ip * nb + ib];
			b[ip * nb + ib] = b[i * nb + ib];
			if (ii > -1)
				for (j = ii; j <= i - 1; j++)
					sum -= a[i * n + j] * b[j * nb + ib];
			else if (sum) // A nonzero element was encountered,so from now on we will have to do the sums n the loop above.
				ii = i;
			b[i * nb + ib] = sum;
		}
		// Now we do the backsubstitution,equation (2.3.7).
		for (i = n - 1; i >= 0; i--)
		{
			sum = b[i * nb + ib];
			for (j = i + 1; j <= n - 1; j++)
				sum -= a[i * n + j] * b[j * nb + ib];
			b[i * nb + ib] = sum / a[i * n + i]; // Store a component of the solution vector X .
		}
		// All done!
	}

	return;
}

int mySolLinear(double *a, int n, double *b, int nb)
/* solve linear equation AX=b with LU decomposition
   A[0..n-1][0..n-1], b[0..n-1][0..nb-1]
   On ouput, a will be replaced by the LU form, 
   and b replaced by X */
// return 1 if succeed, else 0
{
	double d;
	int *indx;

	indx = (int *)malloc(n * sizeof(int));
	if (ludcmp(a, n, indx, &d) != 1)
	{
		printf("linear equation solve fails.\n");
		free(indx);
		return (0);
	}

	lubksb(a, n, indx, b, nb);

	free(indx);

	return (1);
}

double myDet(double *a, int n)
// Calculate the determinant of matrix a[0..n-1][0..n-1]
{
	double d;
	int *indx;
	int i;

	indx = (int *)malloc(n * sizeof(int));
	if (ludcmp(a, n, indx, &d) != 1)
	{
		free(indx);
		return (0.0);
	}

	for (i = 0; i < n; i++)
		d *= a[i * n + i];

	free(indx);

	return (d);
}
