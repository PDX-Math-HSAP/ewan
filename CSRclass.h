#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;

#ifndef __CSRclass_H_INCLUDED__   
#define __CSRclass_H_INCLUDED__ 


//helper function declaration

void accumulate(size_t* array, size_t start, size_t end); //Sets array[i] = array[start]+...+array[i] for i in [start,end).


//Vec class declaration
template<class T>
class Vec
{
public:
  size_t size;
  T* entry;
  bool talk;
  Vec(size_t dim);
  ~Vec();
  void destroy();
  void accumulate();
  void print(string name);
};

#include "Vec.cc"

//CSRMat class declaration
template<class T>
class CSRMat//Data Structure for compressed row format matrix.
{
public:
  bool talk; //If true, CSRMat will cout a ping whenever it allocates or deallocates memory. 
  bool data;//If true, CSRMat contains nonzero entries stored as doubles. If false, no values are stored and each getEntry call has value (T)1.
  bool empty;//If true, CSRMat does not point to any blocks of memory. 
  bool rowEntry;//If true, matrix rowByNnz has been generated.
  bool entryRow;//If true, matrix nnzByRow has been generated.
  size_t rows;//number of rows
  size_t cols;//number of collumns
  size_t nnz;//number of nonzero entries
  size_t* rowPoint;//pointer to array of row pointers.
  size_t* colIndex;//colIndex[rowPoint[i]+k] returns the collumn index of the kth stored nonzero entry in row i (for k<rowPoint[i+1]-rowPoint[i]).
  T* values;//values[rowPoint[i]+k] returns the entry value of the kth stored nonzero entry in row i (for k<rowPoint[i+1]-rowPoint[i]).
  static CSRMat<size_t> rowByNnz;//associated row by nnz relation binary matrix.
  static CSRMat<size_t> nnzByRow;//associated nnz by row relation binary matrix.
  CSRMat();//creates empty data CSRMat structure.
  CSRMat(string filepath, bool hasData, bool importData, string option = "adj");//Loads CSRMat from entry list textfile consisting of lines with format "i j dat" where i and j are indices and dat has type T. If option = "lap", we the constructed matrix will be a graph laplacian.  If option = "adj" the constructed matrix will be an adjacency matrix. If hasData, target must be a text file with only lines of the form "i j dat" with i and j indices and dat a double.  If not hasData, target must be a text file with only lines of the form "i j" with i and j indices. If not importData, data in target file is ignored and each edge is imported with weight 1.
  CSRMat(size_t rowDim,size_t colDim, size_t numVals, bool hasData);//initializes CSRMat object with given dimensions and number of nonZero entries.  All arrays filled with zeros.
  ~CSRMat();
  T getEntry(size_t n) const;//returns value of the n_th stored entry. O(1)
  T getEntry(size_t row, size_t col) const;//returns value of the entry at (row,col). O(rowPoint[row]-rowPoint[row+1])
  size_t getCol(size_t n) const;//returns col index of the n_th stored entry. O(1)
  size_t getRow(size_t n) const;//returns row index of the n_th stored entry. O(rows)
  bool detectZeros();//returns true if CSRMat stores zero in values.
  void removeZeros();//removes any zero entries stored as non-zero entries. O(rows)
};

#include "CSRclass.cc"

//function declarations


template<class T> //T must have a > operator defined.
size_t* rankSort(T* dataArray,  size_t startIndex, size_t endIndex, string direction="ascend");

template<class T> //T must have a > operator defined.
void permute(T* array, size_t* permutation, size_t startIndex, size_t endIndex);

template<class T>
void destroy(CSRMat<T> matrix);//frees all arrays associated with the matrix.

template<class T>
void initValues(CSRMat<T> matrix);//sets data to true. initializes values array and sets all values to zero.

template<class T>
void writeEdgeList(CSRMat<T> matrix, string filepath, string option = "all" , bool includeData=true ,bool sort = false);//writes an edgelist to the given (absolute) filepath. If option = "upper" or "lower", then only upper or lower triangular entries are recorded.  Otherwise, all entries are recorded. If sort, then entries of each row are sorted.

template<class T>
void sortRowsByCol(CSRMat<T> matrix);//puts the entries of colIndex and values corresponding to each row in ascending order by collumn index

template<class T>
void print(CSRMat<T> matrix); //Prints matrix in compressed sparse row format in O(number nonzero entries).

template<class T>
void print_dense(CSRMat<T> matrix);//Pretty prints the matrix in its full form in O(rows*cols) time. Only use "full" option for matrices under 100 by 100.

template<class T>
CSRMat<size_t> generateRowByNnz(CSRMat<T> matrix);

template<class T>
CSRMat<size_t> generateNnzByRow(CSRMat<T> matrix);

template<class T>
CSRMat<T> identity(CSRMat<T> matrix, string side="right");

template<class T, class S>
Vec<S> spMv(CSRMat<T> matrix,Vec<S> vector);

template<class T>
CSRMat<T> spMt(CSRMat<T> matrix);

template<class T>
CSRMat<T> spMM(CSRMat<T> matrix1,CSRMat<T> matrix2);

template<class T>
T dot(Vec<T> vector1,Vec<T> vector2);

template<class T>
Vec<T> edgeWeights(CSRMat<T> graph, CSRMat<size_t> edgeVert, CSRMat<size_t> vertEdge,int weightType);

template<class T>
Vec<size_t> parLuby(CSRMat<T> graph, CSRMat<size_t> edgeVert, CSRMat<size_t> vertEdge);

template<class T>
CSRMat<T> contract(CSRMat<T> graph, CSRMat<size_t> edgeVert, CSRMat<size_t> vertEdge, Vec<size_t> sortedEdgeList);

#include "auxSort.cc"

#include "CSRfun.cc"

#endif
