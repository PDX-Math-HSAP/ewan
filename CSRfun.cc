///////CSRMat function implementations///////

template<class T>
void destroy(CSRMat<T> matrix)//destroys all arrays associated with the matrix.
{
  delete[]matrix.rowPoint;
  delete[]matrix.colIndex;
  if (matrix.data)
  {
    delete[]matrix.values;
  }
  if (matrix.talk)
  {
    if (matrix.data)
    {
      cout << "Freed memory for " << matrix.rows << " by " << matrix.cols << " matrix with " << matrix.nnz << " non-zero double precision entries.\n";
    }
    else
      {
        cout << "Freed memory for " << matrix.rows << " by " << matrix.cols << " binary matrix with " << matrix.nnz << " non-zero entries.\n";
      }
  }
  matrix.empty = true;
  matrix.data  = false;
  matrix.rows  = 0;
  matrix.cols  = 0;
  matrix.nnz   = 0;
};

template<class T>
void initValues(CSRMat<T> matrix)//sets data to true. initializes values array and sets all values to zero.
{
  if (matrix.data)
  {
    delete[]matrix.values;
  }
  matrix.data = true;
  matrix.values   = new T[matrix.nnz];
  fill_n(matrix.values,matrix.nnz,0);
  if (matrix.talk)
  {
    if (matrix.data)
    {
      cout << "Allocated memory for " << matrix.rows << " by " << matrix.cols << " matrix with " << matrix.nnz << " non-zero double precision entries.\n";
    }
    else
    {
      cout << "Allocated memory for " << matrix.rows << " by " << matrix.cols << " binary matrix with " << matrix.nnz << " non-zero entries.\n";
    }
  }
  matrix.empty = false;
};

template<class T>
void writeEdgeList(CSRMat<T> matrix, string filepath, string option, bool includeData,bool sort)//writes an edgelist to the given (absolute) filepath. If option = "upper" or "lower", then only upper or lower triangular entries are recorded.  Otherwise, all entries are recorded. If sort, then entries of each row are sorted.
{
  if (sort)
  {
    sortRowsByCol(matrix);
  }
  ofstream outfile (filepath,ofstream::out);
  for (size_t i=0;i<matrix.rows;i++)
  {
    for (size_t j=matrix.rowPoint[i];j<matrix.rowPoint[i+1];j++)
    {
      if (
           ((option != "upper" ) && (option != "lower"))
        || ((option == "lower")&&(matrix.colIndex[j]<i))
        || ((option == "upper")&&(matrix.colIndex[j]>i))
             )
      {
        if (includeData)
        {
          outfile << i <<" "<<matrix.colIndex[j]<<" "<<matrix.values[j]<<"\r\n";
        }
        else
        {
          outfile << i <<" "<<matrix.colIndex[j]<<"\r\n";
        }
      }
    }
  }
  outfile.close();
};

template<class T>
void sortRowsByCol(CSRMat<T> matrix)//puts the entries of colIndex and values corresponding to each row in ascending order by collumn index
{
  for (size_t i=0; i<matrix.rows; i++)
  {
    size_t* sortedIndices = rankSort(matrix.colIndex,matrix.rowPoint[i],matrix.rowPoint[i+1]);
    permute(matrix.colIndex,sortedIndices,matrix.rowPoint[i],matrix.rowPoint[i+1]);
    if (matrix.data)
    {
      permute(matrix.values,sortedIndices,matrix.rowPoint[i],matrix.rowPoint[i+1]);
    }
  }
};

template<class T>
void print(CSRMat<T> matrix)//Prints matrix in compressed sparse row format in O(number nonzero entries).
{
    cout<<"Printing " << matrix.rows << " by " << matrix.cols
        << " matrix with " << matrix.nnz << " non-zero entries in CSR format.\n";

    cout << "I = ["  ;
    for (size_t j = 0; j <= matrix.rows ; j++ )
    {
      if (j == matrix.rows) {cout << matrix.rowPoint[j]<< "]\n";}
      else {cout << matrix.rowPoint[j] << ",";}
    }

    cout << "J = ["  ;
    for (size_t j = 0; j < matrix.nnz; j++ )
    {
      if (j == matrix.nnz-1) {cout << matrix.colIndex[j]<< "]\n";}
      else {cout << matrix.colIndex[j] << ",";}
    }

    cout << "D = ["  ;
    for (size_t j = 0; j < matrix.nnz; j++ )
    {
      if (j == matrix.nnz-1) {cout << matrix.getEntry(j)<< "]\n";}
      else {cout << matrix.getEntry(j) << ",";}
    }
};

template<class T>
void print_dense(CSRMat<T>& matrix)//Pretty prints the matrix in its full form in O(rows*cols) time. Only use "full" option for matrices under 100 by 100.
{
 cout<<"Printing " << matrix.rows << " by " << matrix.cols
		<< " matrix with " << matrix.nnz << " non-zero entries.\n";
	sortRowsByCol(matrix);
	size_t entry = 0;
	for (size_t i=0;i<matrix.rows;i++)
	{
	  entry=0;
	  for (size_t j=0;j<matrix.cols;j++)
	  {
		if (matrix.colIndex[matrix.rowPoint[i]+entry]==j && entry<(matrix.rowPoint[i+1]-matrix.rowPoint[i]))
		{
		  cout << matrix.getEntry(matrix.rowPoint[i] + entry++);
		}
		else
		{
		  cout <<"0";
		}
		if (j != matrix.cols-1)
		{
		  cout<<",";
		}
		cout << " ";
	  }
	  cout<<"\n";
	}
}

template<class T>
CSRMat<size_t> generateRowByNnz(CSRMat<T> matrix)//generates row by entry relation (only for upper triangular portion of the matrix).
{
	if (matrix.talk) {
		cout << "generating RowByNnz relation\n\r";
	}
	size_t edgeCount = 0; //count edges
	for (size_t i = 0; i < matrix.rows; i++) {
		for (size_t j = matrix.rowPoint[i]; j < matrix.rowPoint[i + 1]; j++) {
			if (matrix.colIndex[j] > i) {
				edgeCount++;
			}
		}
	}

	CSRMat<size_t> rowByNnz(matrix.rows, edgeCount, 2 * edgeCount, false); //rowByNnz is a binary matrix with 2 nonzero entries for each edge.

	if (matrix.talk) {
		cout << "initialized rowByNnz with " << matrix.rows << " rows and " << edgeCount
				<< " cols\n\r";
	}
	for (size_t i = 0; i < matrix.rows; i++) {
		for (size_t j = matrix.rowPoint[i]; j < matrix.rowPoint[i + 1]; j++) {
			if (matrix.colIndex[j] > i) {
				rowByNnz.rowPoint[i + 1]++;
				rowByNnz.rowPoint[matrix.colIndex[j] + 1]++;
			}
		}
	}
	//cout << " accum\n\r";
	accumulate(rowByNnz.rowPoint, 0, rowByNnz.rows + 1);

	//cout << " init counter\n";
	size_t* counter = new size_t[matrix.rows]; //counter[i] will track number of i-edges processed so far.
	fill_n(counter, matrix.rows, 0);

	size_t edge = 0;

	//cout << " assign colIndex\n\r";
	for (size_t i = 0; i < matrix.rows; i++) {
		for (size_t j = matrix.rowPoint[i]; j < matrix.rowPoint[i + 1]; j++) {
			if (matrix.colIndex[j] > i) {
				rowByNnz.colIndex[rowByNnz.rowPoint[i] + counter[i]] = edge;
				rowByNnz.colIndex[rowByNnz.rowPoint[matrix.colIndex[j]]
						+ counter[matrix.colIndex[j]]] = edge;
				counter[i]++;
				counter[matrix.colIndex[j]]++;
				edge++;
			}
		}
	}
	//cout << "deleting counter\n\r";
	delete[] counter;
	//cout << "rowByNnz constructed\n\r";
	return rowByNnz;
}
;

template<class T>
CSRMat<size_t> generateNnzByRow(CSRMat<T> matrix)//generates entry by row relation (only for upper triangular portion of the matrix).
{
	if (matrix.talk) {
			cout << "generating nnzByRow relation\n\r";
		}
		size_t edgeCount = 0; //count edges
		for (size_t i = 0; i < matrix.rows; i++) {
			for (size_t j = matrix.rowPoint[i]; j < matrix.rowPoint[i + 1]; j++) {
				if (matrix.colIndex[j] > i) {
					edgeCount++;
				}
			}
		}

		CSRMat<size_t> nnzByRow(edgeCount,matrix.rows , 2 * edgeCount, false); //rowByNnz is a binary matrix with 2 nonzero entries for each edge.

		if (matrix.talk) {
			cout << "initialized nnzByRow with " << nnzByRow.rows << " rows and "<< nnzByRow.cols << " cols\n\r";
		}
		for (size_t i = 0; i < edgeCount+1; i++) {
			nnzByRow.rowPoint[i] = 2*i;
		}
		edgeCount = 0;
		for (size_t i = 0; i < matrix.rows; i++) {
					for (size_t j = matrix.rowPoint[i]; j < matrix.rowPoint[i + 1]; j++) {
						if (matrix.colIndex[j] > i) {
							nnzByRow.colIndex[nnzByRow.rowPoint[edgeCount]] = i;
							nnzByRow.colIndex[nnzByRow.rowPoint[edgeCount++]+1]=matrix.colIndex[j];
						}
					}
				}
		if (matrix.talk) {
			cout << "nnzByRow constructed\n\r";
		}
		return nnzByRow;
}
;

template<class T>
CSRMat<T> identity(CSRMat<T> matrix, string side)
{
  size_t n=0;
  if (side == "left"){n=matrix.rows;}
  else {n=matrix.cols;}
  
  CSRMat<T> id(n,n,n, false);
  for (size_t i=0; i<n; i++)
  {
    id.rowPoint[i] = i;
    id.colIndex[i] = i;
  }
  id.rowPoint[n] = n;
  return id;
}

template<class T, class S>
Vec<S> spMv(CSRMat<T> matrix,Vec<S> vector)// CSR Sparse Matrix Dense Vecor Product
{
  if (matrix.talk || vector.talk)
  {
    cout<< "calling spMv\n";
  }
  Vec<S> product(matrix.rows);
  fill_n(product.entry,matrix.rows,0);
  for (size_t k=0; k<matrix.rows; k++) {
    for (size_t j=matrix.rowPoint[k]; j<matrix.rowPoint[k+1] ; j++) {
      product.entry[k] += matrix.getEntry(j)*vector.entry[matrix.colIndex[j]];
    }
  }
  return product;
} 

template<class T>
CSRMat<T> spMt(CSRMat<T> matrix)// CSR Sparse Matrix Transpose
{
  if (matrix.talk)
  {
    cout<< "calling spMt\n";
  }
  //initialize csr matrix "transpose"
  CSRMat<T> transpose(matrix.cols, matrix.rows, matrix.nnz, matrix.data);
  
  //initialize counter   
  size_t* counter = new size_t[transpose.rows];
  fill_n(counter,transpose.rows,0);
  
  //count nnz's in transpose by rows
  for (size_t j=0; j<matrix.nnz ; j++)
  {
    transpose.rowPoint[matrix.colIndex[j]+1]++;
  }
  //Accumulate row pointers
  accumulate(transpose.rowPoint,0,transpose.rows+1);
  
  //compute collumn indices and corresponding values
  for (size_t i = 0; i<matrix.rows; i++)
  {
    for(size_t j = matrix.rowPoint[i]; j< matrix.rowPoint[i+1]; j++)
    {
      size_t l = transpose.rowPoint[matrix.colIndex[j]]+counter[matrix.colIndex[j]];
      transpose.colIndex[l] = i;
      if (matrix.data)
      {
        transpose.values[l] = matrix.values[j];
      }
      counter[matrix.colIndex[j]]++;
    }
  }
  delete[]counter;
  return transpose;
} 

template<class T>
CSRMat<T> spMM( CSRMat<T> matrix1,CSRMat<T> matrix2)//CSR sparse matrix product
{
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "calling spMM\n->initRowArrays\n\r";
  }
  
  
   //initialize row pointer array and flag array
  size_t* rowPoint  = new size_t[matrix1.rows+1];
  fill_n(rowPoint,matrix1.rows+1,0);
  size_t* flags     = new size_t[matrix2.cols];
  fill_n(flags,matrix2.cols,0);
  
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "initProduct\n";
  }

  size_t iCounter = 0;
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "rowPoint\n";
  }
  for (size_t i=0;i<matrix1.rows ; i++)
  {
    rowPoint[i] = iCounter;
    for (size_t jPtr = matrix1.rowPoint[i]; jPtr<matrix1.rowPoint[i+1]; jPtr++)
    {
      size_t j = matrix1.colIndex[jPtr];
      for (size_t lPtr = matrix2.rowPoint[j]; lPtr < matrix2.rowPoint[j+1]; lPtr++)
      {
        size_t l = matrix2.colIndex[lPtr];
        if (flags[l] != i+1)
        {
          iCounter++;
          flags[l] = i+1;
        }
      }
    }
  }
  rowPoint[matrix1.rows] = iCounter;
  
  //initialize CSRMat product and reinitialize flag array
  CSRMat<T> product(matrix1.rows, matrix2.cols, iCounter, true);
  product.rowPoint = rowPoint;
  fill_n(flags,product.cols,0);
  
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "init valTemp\n";
  }
  
  iCounter = 0;
  T* valTemp = new T[product.cols];
  for (size_t i=0;i<product.rows ; i++)
  {
    fill_n(valTemp,product.cols,0);
    for (size_t jPtr = matrix1.rowPoint[i]; jPtr < matrix1.rowPoint[i+1]; jPtr++)
    {
      size_t j = matrix1.colIndex[jPtr];
      T valA = matrix1.getEntry(jPtr);
      for (size_t lPtr = matrix2.rowPoint[j]; lPtr < matrix2.rowPoint[j+1]; lPtr++)
      {
        size_t l = matrix2.colIndex[lPtr];
        valTemp[l] += valA * matrix2.getEntry(lPtr);
        if (flags[l] != i+1 )
        {
          product.colIndex[iCounter] = l;
          iCounter++;
          flags[l] = i+1;
        } 
      }
    }
    for (size_t j = product.rowPoint[i]; j < product.rowPoint[i+1]; j++ )
    {
      product.values[j] = valTemp[product.colIndex[j]];
    }
  }
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "deleting flags/valTemp\n";
  }
  delete[]flags;
  delete[]valTemp;
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "finished spMM\n";
  }
  return product;
} 

template<class T>
T dot(Vec<T> vector1,Vec<T> vector2)//dot product for Vec class
{
  size_t dim = std::min(vector1.size,vector2.size);
  double dot = 0; 
  for (size_t i = 0; i<dim; i++){
    dot += vector1.entry[i]*vector2.entry[i];
  }
  return dot;
}

template<class T>
Vec<T> edgeWeights(CSRMat<T> graph, CSRMat<size_t> edgeVert, CSRMat<size_t> vertEdge, int weightType)//returns array of edge weights (upper triangular entries in row major order).  dataType is 0 for "uniform", 1 for "random", 2 for "inverse edge degree",or 3 for "given".  Option 2 requires both rowByNnz and nnzByRow to be generated.
{
  if (graph.talk)
  {
    cout<<"calling edgeWeights->counting edges";
  }
  size_t edgeCounter = 0;
  for (size_t i=0; i<graph.rows; i++)//counting edges
  {
    for (size_t j=graph.rowPoint[i];j<graph.rowPoint[i+1];j++)
    {
      if (graph.colIndex[j] > i)
      {
        edgeCounter++;
      }
    }    
  }
  if (graph.talk)
    {
	  cout<<"->edgeCounter = "<<edgeCounter;
    }


  Vec<T> edgeWeights(edgeCounter);
  
  
  if (graph.talk)
  {
    cout<<"->calculating edgeWeights";
  }
  size_t edge=0;
  switch (weightType){
  	case 0:
  	{
  	for (size_t i=0; i<graph.rows; i++)
  	  {
  	   for (size_t j=graph.rowPoint[i];j<graph.rowPoint[i+1];j++)
  	   {
  		if (graph.colIndex[j] > i)
  	    {
  		edgeWeights.entry[edge++] = 1;

  	    }
  	   }
  	  }
  	break;
  	}
  	case 1:
  	{
  	  	for (size_t i=0; i<graph.rows; i++)
  	  	  {
  	  	   for (size_t j=graph.rowPoint[i];j<graph.rowPoint[i+1];j++)
  	  	   {
  	  		if (graph.colIndex[j] > i)
  	  	    {
  	  		edgeWeights.entry[edge++] = (T)(rand() % 1000) /1000;

  	  	    }
  	  	   }
  	  	  }
  	break;
  	}
  	case 2:
  	{
  		if (graph.talk)
  		{
  		  cout<<"->initing vectors";
  		}
  		Vec<T> ones(vertEdge.cols);
  		fill_n(ones.entry,vertEdge.cols,1);
  		Vec<T> edgeDegrees = spMv(edgeVert, spMv(vertEdge, ones));
  	  	for (size_t i=0; i<graph.rows; i++)
  	  	  {
  	  	   for (size_t j=graph.rowPoint[i];j<graph.rowPoint[i+1];j++)
  	  	   {
  	  		if (graph.colIndex[j] > i)
  	  	    {
  	  		 T edgeDegree = static_cast<T>(edgeDegrees.entry[edge]);
  	  		 edgeWeights.entry[edge++] = 1/(edgeDegree-1); // edgeDegree[e] includes e twice.
  	  	    }
  	  	   }
  	  	  }
  	break;
  	}
  	default:
  	{
  	for (size_t i=0; i<graph.rows; i++)
     {
      for (size_t j=graph.rowPoint[i];j<graph.rowPoint[i+1];j++)
  	  {
  	   if (graph.colIndex[j] > i)
  	   {
  		edgeWeights.entry[edge++] = graph.values[j];
  	   }
  	  }
  	 }
  	}
  }

  if (graph.talk)
  {
    cout<<"->edgeWeights computed\n\r";
    edgeWeights.print("edgeWeights:");
  }
  return edgeWeights;  
}

template<class T>
Vec<size_t> parLuby(CSRMat<T> graph, CSRMat<size_t> edgeVert, CSRMat<size_t> vertEdge)//parallel style Luby's Algorithm
{
  if (graph.talk) {cout<<"calling parLuby\n";}
  Vec<T> edgeData = edgeWeights(graph,edgeVert, vertEdge,2);//edge weights.  Parameter is 0 for "uniform", 1 for "random", 2 for "inverse of edge degree",or 3 for "given"
  //note that edgeData.entry[i] is the weight of edge i in either relation above.
  
  size_t* rankEdge = new size_t[vertEdge.cols];
  
  for (size_t i=0;i<vertEdge.cols;i++)
  {
    rankEdge[i]=i;
  }
  
  if (graph.talk) {cout<<" Sorting...\n";}
  sort(&rankEdge[0],&rankEdge[vertEdge.cols],[&edgeData] (size_t a, size_t b)
  {
    return edgeData.entry[a]<edgeData.entry[b];
  });//rankEdge[i] = (i+1)st heaviest edge
  if (graph.talk) {cout<<" finished Sort\n"<<"heaviest edge is "
                        << rankEdge[0]<<" with weight "<< edgeData.entry[rankEdge[0]]<<"\n";}
  if (graph.talk) {cout<<" init matching...\n";}
  Vec<size_t> matching(edgeData.size); //empty list to hold the matching edge names.
  size_t matchCount = 0; //will count edges in matching.
  if (graph.talk) {cout<<" computing matching...\n";}
  for (size_t i=0; i<edgeData.size; i++)//loop through all edges 
  {
    if (edgeData.entry[rankEdge[i]]!=0)//if (i+1)st largest edge is active
    { 
      matching.entry[matchCount] = rankEdge[i];//add index of this active edge to matching list
      matchCount++;//increment match counter
      size_t endPoint1 = edgeVert.colIndex[edgeVert.rowPoint[rankEdge[i]]];//first endpoint of matched edge
      size_t endPoint2 = edgeVert.colIndex[edgeVert.rowPoint[rankEdge[i]]+1];//second endpoint of matched edge
      
      for (size_t j=vertEdge.rowPoint[endPoint1]; j<vertEdge.rowPoint[endPoint1+1]; j++)//loop through edges sharing endPoint1
      {
        edgeData.entry[vertEdge.colIndex[j]]=0;//set to inactive
      }
      
      for (size_t j=vertEdge.rowPoint[endPoint2]; j<vertEdge.rowPoint[endPoint2+1]; j++)//loop through edges sharing endPoint2
      {
        edgeData.entry[vertEdge.colIndex[j]]=0;
      }
    }
  }
  

  edgeData.destroy();
  delete[]rankEdge;
  
  if (graph.talk) {cout<<" init match...\n";}
  Vec<size_t> match(matchCount);
  for (size_t i=0;i<matchCount;i++)
  {
    match.entry[i]=matching.entry[matchCount-i-1];
  } 
  
  matching.destroy();
  if (graph.talk) { match.print("match: ");
    cout<<" returning match...\n";}
  return match;
}

template<class T>
CSRMat<T> contract(CSRMat<T> graph, CSRMat<size_t> edgeVert, CSRMat<size_t> vertEdge, Vec<size_t> edgeList)//Graph is assumed to be a symmetric matrix.  EdgeList is assumed to be cycle free.  That is edges must comprise a sub-forest of the graph.
{
  if (graph.talk)
  {
    cout<< "calling contract\n";
  }
  size_t newSize = graph.cols - edgeList.size; // will contract one vertex for each edge in edgeList
  CSRMat<T> projection(graph.cols,newSize,graph.cols, false);

  size_t edge=0;

  for (size_t i=0;i<=projection.rows;i++)//projection will have one entry per row.
  {
    projection.rowPoint[i] = i;
  }
  //cout<< "sort\n";
  sort(&edgeList.entry[0],&edgeList.entry[edgeList.size], [&edgeVert] (size_t a, size_t b)
      {
        size_t maxEndPointa = max(edgeVert.colIndex[edgeVert.rowPoint[a]], edgeVert.colIndex[edgeVert.rowPoint[a]+1]);
        size_t maxEndPointb = max(edgeVert.colIndex[edgeVert.rowPoint[b]], edgeVert.colIndex[edgeVert.rowPoint[b]+1]);
        return maxEndPointa < maxEndPointb;
      });//sort edgeList by largest end point.
  Vec<size_t> matchCount(graph.cols); //matchCount.entry[i] counts number of matched edges with a max endpoint less than or equal to i.
  for (size_t i = 0; i < edgeList.size; i++)
  {
    size_t currEdge    = edgeList.entry[i];
    size_t minEndPoint = min(edgeVert.colIndex[edgeVert.rowPoint[currEdge]], edgeVert.colIndex[edgeVert.rowPoint[currEdge]+1]);
    size_t maxEndPoint = max(edgeVert.colIndex[edgeVert.rowPoint[currEdge]], edgeVert.colIndex[edgeVert.rowPoint[currEdge]+1]);
    if (minEndPoint<maxEndPoint)
    {
      matchCount.entry[maxEndPoint]++;
    }
  }
  
  matchCount.accumulate();
  
 // matchCount.print("matchCount");
  size_t matchPos = 0;
  //edgeList.print("edgeList: ");
  for (size_t i = 0; i < graph.cols && matchPos<edgeList.size; i++)
  {
    //cout<<"edgeList has size "<<edgeList.size<<"\n\r";
	//cout<<"matchPos= "<<matchPos<<" and indexes edgeList\n\r";
	size_t currEdge = edgeList.entry[matchPos];
    //cout<< "min2\n";
    //cout<< "currEdge = "<< currEdge <<"\n\r";
    size_t minEndPoint = min(edgeVert.colIndex[edgeVert.rowPoint[currEdge]], edgeVert.colIndex[edgeVert.rowPoint[currEdge]+1]);
    size_t maxEndPoint = max(edgeVert.colIndex[edgeVert.rowPoint[currEdge]], edgeVert.colIndex[edgeVert.rowPoint[currEdge]+1]);
    //cout<< i <<": contracting ("<<minEndPoint<<","<<maxEndPoint<<")\n";
    
    if (maxEndPoint == i)//maxEndPoint is contracted into vertex (minEndPoint - matchCount.entry[minEndPoint])
    {
      projection.colIndex[projection.rowPoint[i]] = minEndPoint-matchCount.entry[minEndPoint];
      //cout<< minEndPoint<<"-"<<matchCount.entry[minEndPoint]<<"-a\n";
      matchPos++;
    }
    else//
    {
      projection.colIndex[projection.rowPoint[i]] = i-matchCount.entry[i];
      //cout<< i<<"-"<<matchCount.entry[i] <<"-b\n";
    }
    //cout<<"done!\n";
  }
  matchCount.destroy();
  return projection;
}
