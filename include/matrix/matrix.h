#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <mpfr.h>

/**
 * A matrix POD structure.
 */
struct mtx
{
	mpfr_t* storage;	///!< internal storage
	size_t nrows;		///!< number of rows
	size_t ncols;		///!< number of columns
};

/**
 * Initializes a matrix.
 * 
 * Allocates a memory for internal storage and matrix elements.
 * 
 * @param m a matrix
 * @param rows a number of the rows
 * @param columns a number of the columns
 * @param prec a precision
 * @return 0 on success, non-zero otherwise
 */
int mtx_init(struct mtx* const m, size_t rows, size_t columns, mpfr_prec_t prec);

/**
 * Clears a matrix resources.
 * 
 * Destroys the matrix elements end frees the internal storage.
 * 
 * @param m a matrix
 * @return 0 on success, non-zero otherwise
 */
int mtx_clear(struct mtx const m);

/**
 * Outputs a matrix to stream.
 * 
 * @param stream a stream
 * @param m a matrix
 * @return The number of characters written, or 0 if an error occurred.
 */
int mtx_fprint(FILE* stream, struct mtx const m);

/**
 * Sets the matrix elements values.
 * 
 * @param m a matrix
 * @param val a value
 * @param diagval a value will be used only for diagonal elements
 * if matrix is square
 * @return 0 on success, non-zero otherwise
 */
int mtx_fill(struct mtx m, double val, double diagval);

/**
 * Multiplies the matrices.
 * 
 * The dimensions of the matrices must be as follows:
 *		rop.nrows == op1.nrows
 *		rop.ncols == op2.ncols
 *		op1.ncols == op2.nrows
 * 
 * @param rop the result matrix
 * @param op1 a left-hand side operand
 * @param op2 a right-hand side operand
 * @return 0 on success, non-zero otherwise
 */
int mtx_mul(struct mtx rop, struct mtx const op1, struct mtx const op2);

/**
 * Multiplies the matrix by value.
 * 
 * @param rop the result matrix
 * @param op1 a matrix
 * @param op2 a value
 * @return 0 on success, non-zero otherwise
 */
int mtx_mulval(struct mtx rop, struct mtx const op1, mpfr_t op2);

/**
 * Adds the matrices.
 * 
 * The dimensions of the matrices must be the same.
 * 
 * @param rop the result matrix
 * @param op1 the first matrix
 * @param op2 the second matrix
 * @return 0 on success, non-zero otherwise
 */
int mtx_add(struct mtx rop, struct mtx const op1, struct mtx const op2);

/**
 * Transposes the matrix.
 * 
 * The dimension of the rop matrix must be as follows:
 *		rop.nrows == op.ncols
 *		rop.ncols == op.nrows
 * 
 * @param rop the result matrix
 * @param op a matrix to be transposed
 * @return 0 on success, non-zero otherwise
 */
int mtx_tr(struct mtx rop, struct mtx const op);

#endif /* MATRIX_H */
