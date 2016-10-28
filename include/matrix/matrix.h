#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>

/**
 * Initializes a matrix.
 * 
 * @param mtx a matrix to be initialized
 * @param nrows number of rows
 * @param ncols number of columns
 * @return 0 on success, non-zero otherwise
 */
int mtx_init(mpq_t** mtx, size_t const nrows, size_t const ncols);

/**
 * Clears a matrix resources.
 * 
 * @param mtx a matrix to be cleared
 * @param nrows number of rows
 * @param ncols number of columns
 * @return 0 on success, non-zero otherwise
 */
int mtx_clear(mpq_t** mtx, size_t const nrows, size_t const ncols);

/**
 * Multiplies the matrices.
 * 
 * @param res the multiplication result
 * @param a a left-hand side matrix
 * @param b a right-hand side matrix
 * @param arows a number of rows in matrix a
 * @param bcols a number of columns in matrix b
 * @param cols a number of columns in the matrices a and b
 * @return 0 on success, non-zero otherwise
 */
int mtx_mul(mpq_t* res, mpq_t* a, mpq_t* b, size_t const arows, size_t const bcols, size_t const cols);

size_t mtx_print(FILE* stream, mpq_t* mtx, size_t const nrows, size_t const ncols, int base);

#endif /* MATRIX_H */
