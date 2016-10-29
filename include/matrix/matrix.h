#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <mpfr.h>

/**
 * A matrix POD structure.
 */
struct matrix
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
int mtx_init(struct matrix* const m, size_t rows, size_t columns, mpfr_prec_t prec);

/**
 * Clears a matrix resources.
 * 
 * Destroys the matrix elements end frees the internal storage.
 * 
 * @param m a matrix
 * @return 0 on success, non-zero otherwise
 */
int mtx_clear(struct matrix const* const m);

/**
 * Outputs a matrix to stream.
 * 
 * @param stream a stream
 * @param m a matrix
 * @return The number of characters written, or 0 if an error occurred.
 */
int mtx_print(FILE* stream, struct matrix const* const m);

/**
 * Multiplies the matrices.
 * 
 * The number of columns of the left-hand side operand must be equal
 * to the number of the rows of the right-hand side operand.
 * 
 * @param rop the result matrix
 * @param op1 a left-hand side operand
 * @param op2 a right-hand side operand
 * @return 0 on success, non-zero otherwise
 */
int mtx_mul(struct matrix* const rop, struct matrix const* const op1, struct matrix const* const op2);

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
int mtx_add(struct matrix* const rop, struct matrix const* const op1, struct matrix const* const op2);

#endif /* MATRIX_H */
