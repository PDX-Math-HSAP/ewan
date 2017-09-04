//class method definitions

template<class T>
CSRMat<T>::CSRMat() //creates empty data CSRMat structure.
{
	talk = false;
	data = false;
	empty = true;
	rows = 0;
	cols = 0;
	nnz = 0;
	rowEntry = false;
	entryRow = false;
	rowPoint = NULL;
	colIndex = NULL;
	values = NULL;
};

template<class T>
CSRMat<T>::CSRMat(string filepath, bool hasData, bool importData, string option) //Loads CSRMat from entry list textfile consisting of lines with format "i j dat" where i and j are indices and dat has type T. If option = "lap", we the constructed matrix will be a graph laplacian.  If option = "adj" the constructed matrix will be an adjacency matrix. If hasData, target must be a text file with only lines of the form "i j dat" with i and j indices and dat a double.  If not hasData, target must be a text file with only lines of the form "i j" with i and j indices. If not importData, data in target file is ignored and each edge is imported with weight 1.
{

	  data = hasData;
	  rowEntry = false;
	  entryRow = false;
	  talk = false;
	  ifstream iFile(filepath);

	  if (!iFile)
	  {
	    cout<< "Read Failed!\n";
	  }
	  else
	  {
	   empty=false;
	  }

	  size_t maxI = 0;
	  size_t x;
	  size_t y;

	  if (data)
	  {
	    double z;

	    while (iFile >> x >> y >> z)
	    {
	      if (x!=y)
	      {
	    	    nnz+=2 ;
	        maxI = max(maxI,max(x,y));
	      }
	    }

	    rows = 1+maxI;
	    cols = rows;
	    size_t* jCounter = new size_t[rows];
	    rowPoint = new size_t[rows+1];
	    if (option == "lap")
	    {
	    	data = true;
	    	nnz += rows;
	      fill_n(jCounter,rows,1);
	      fill_n(rowPoint,rows+1,1);
	      rowPoint[0]=0;
	    }
	    else
	    {
	      fill_n(jCounter,rows,0);
	      fill_n(rowPoint,rows+1,0);
	    }
	    colIndex = new size_t[nnz];
	    values   = new double[nnz];



	    if (talk)
	    {
	      cout << "Allocated memory for " << rows << " by " << cols << " matrix with " << nnz << " non-zero double precision entries.\n";
	    }

	    iFile.clear();
	    iFile.seekg(0, iFile.beg);


	    while (iFile >> x >> y >> z)
	    {
	      if (x!=y)
	      {
	    	  rowPoint[x+1]++;
	    	  rowPoint[y+1]++;
	      }
	    }

	    accumulate(rowPoint,0,rows+1);

	    iFile.clear();
	    iFile.seekg(0, iFile.beg);

	    if (option == "lap" && !importData)
	    {
	      while (iFile >> x >> y >> z)
	      {
	        if (x!=y)
	        {
	        	colIndex[rowPoint[x]+jCounter[x]] = y;
	        	colIndex[rowPoint[y]+jCounter[y]] = x;
	        values[rowPoint[x]+jCounter[x]]  = -1;
	        values[rowPoint[y]+jCounter[y]]  = -1;
	        jCounter[x]++;
	        jCounter[y]++;
	        }
	      }
	    }
	    else if (option == "lap" && importData)
	    {
	      while (iFile >> x >> y >> z)
	      {
	        if (x!=y)
	        {
	        	colIndex[rowPoint[x]+jCounter[x]] = y;
	        	colIndex[rowPoint[y]+jCounter[y]] = x;
	        	values[rowPoint[x]+jCounter[x]]  = -z;
	        	values[rowPoint[y]+jCounter[y]]  = -z;
	        jCounter[x]++;
	        	jCounter[y]++;
	        }
	      }
	    }
	    else if (!importData)
	    {
	      while (iFile >> x >> y >> z)
	      {
	        if (x!=y)
	        {
	        	colIndex[rowPoint[x]+jCounter[x]] = y;
	        	colIndex[rowPoint[y]+jCounter[y]] = x;
	        	values[rowPoint[x]+jCounter[x]]  = 1;
	        	values[rowPoint[y]+jCounter[y]]  = 1;
	        jCounter[x]++;
	        jCounter[y]++;
	        }
	      }
	    }
	    else
	    {
	      while (iFile >> x >> y >> z)
	      {
	        if (x!=y)
	        {
	        	colIndex[rowPoint[x]+jCounter[x]] = y;
	        	colIndex[rowPoint[y]+jCounter[y]] = x;
	        	values[rowPoint[x]+jCounter[x]]  = z;
	        	values[rowPoint[y]+jCounter[y]]  = z;
	        jCounter[x]++;
	        jCounter[y]++;
	        }
	      }
	    }

	    if (option == "lap")//add in diagonal entry for each row.
	    {
	      for (size_t i = 0; i<rows; i++)
	      {
	    	  colIndex[rowPoint[i]]=i;
	        for (size_t j = rowPoint[i]+1; j<rowPoint[i+1]; j++)
	        {
	        	values[rowPoint[i]] += -values[j];
	        }
	      }
	    }


	    delete[]jCounter;
	    iFile.close();

	  }
	  else
	  {
	    while (iFile >> x >> y)
	    {
	      if (x!=y)
	      {
	    	  nnz+=2 ;
	      maxI = max(maxI,max(x,y));
	      }
	    }

	    rows = 1+maxI;
	    cols = rows;
	    size_t* jCounter = new size_t[rows];
	    rowPoint = new size_t[rows+1];
	    if (option == "lap")
	    {
	    	data = true;
	    	nnz += rows;
	    	values = new T[nnz];
	    	fill_n(jCounter,rows,1);
	    fill_n(values, nnz, 0);
	    fill_n(rowPoint,rows+1,1);
	    rowPoint[0]=0;
	    }
	    else
	    {
	      fill_n(jCounter,rows,0);
	      fill_n(rowPoint,rows+1,0);
	    }
	    colIndex = new size_t[nnz];



	    iFile.clear();
	    iFile.seekg(0, iFile.beg);

	    while (iFile >> x >> y)
	    {
	      if (x!=y)
	      {
	    	  rowPoint[x+1]++;
	    	  rowPoint[y+1]++;
	      }
	    }

	    accumulate(rowPoint,0,rows+1);

	    iFile.clear();
	    iFile.seekg(0, iFile.beg);




	    while (iFile >> x >> y)
	    {
	      if (x<y)
	      {
	    	  colIndex[rowPoint[x]+jCounter[x]] = y;
	    	  colIndex[rowPoint[y]+jCounter[y]] = x;
	        if (option == "lap")
	        {
	        	values[rowPoint[x]+jCounter[x]]  = -1;
	        	values[rowPoint[y]+jCounter[y]]  = -1;
	        }
	        jCounter[x]++;
	        jCounter[y]++;
	      }
	    }



	    if (option == "lap")//add in diagonal entry for each row.
	    {
	      for (size_t i = 0; i<rows; i++)
	      {
	    	  colIndex[rowPoint[i]]=i;
	        for (size_t j = rowPoint[i]+1; j<rowPoint[i+1]; j++)
	        {
	        	values[rowPoint[i]] += -values[j];
	        }
	      }
	    }

	    if (talk && option == "adj")
	    {
	      cout << "Allocated memory for " << rows << " by " << cols << " binary matrix with " << nnz << " non-zero entries.\n";
	    }
	    if (talk && option == "lap")
	    {
	      cout << "Allocated memory for " << rows << " by " << cols << " matrix with " << nnz << " non-zero entries.\n";
	    }

	    iFile.close();
	    delete[]jCounter;
	  }
};

template<class T>
CSRMat<T>::CSRMat(size_t rowDim, size_t colDim, size_t numVals, bool hasData) //initializes CSRMat object with given dimensions and number of nonZero entries.  All arrays filled with zeros.
{
	data = hasData;
	rows = rowDim;
	cols = colDim;
	nnz = numVals;
	talk = false;
	rowPoint = new size_t[rows + 1];
	colIndex = new size_t[nnz];
	rowEntry = false;
	entryRow = false;

	fill_n(rowPoint, rows + 1, 0);
	fill_n(colIndex, nnz, 0);
	if (data) {
		values = new T[nnz];
		fill_n(values, nnz, 0);
	}
	if (talk) {
		if (data) {
			cout << "Allocated memory for " << rows << " by " << cols
					<< " matrix with " << nnz
					<< " non-zero double precision entries.\n";
		} else {
			cout << "Allocated memory for " << rows << " by " << cols
					<< " binary matrix with " << nnz << " non-zero entries.\n";
		}
	}
	empty = false;
};

template<class T>
CSRMat<T>::~CSRMat() //deconstructor.  Does not delete memory pointed to by data structure.
{

}
;

template<class T>
T CSRMat<T>::getEntry(size_t n) const //returns value of the n_th stored entry. O(1)
		{
	if (data) {
		return values[n];
	} else {
		return (T) 1;
	}
}
;

template<class T>
T CSRMat<T>::getEntry(size_t row, size_t col) const //returns value of the n_th stored entry. O(max nnz per row)
{

	for (size_t j=rowPoint[row]; j<rowPoint[row+1]; j++)
	{

	   if (!data && colIndex[j]==col)
	   {
		   return (T) 1;
	   }
	   else if(colIndex[j]==col)
	   {
		return values[j];
	   }
	   else
	   {
		return 0;
	   }
	 }
}
;


template<class T>
size_t CSRMat<T>::getCol(size_t n) const //returns col index of the n_th stored entry. O(1)
		{
	return colIndex[n];
}
;
template<class T>
size_t CSRMat<T>::getRow(size_t n) const//returns row index of the n_th stored entry. O(number of rows) -- avoid using in loops!
		{
	size_t row = 0;
	while (rowPoint[row + 1] < n) {
		row++;
	}
	return row;
}
;


template<class T>
void CSRMat<T>::removeZeros()
{
 size_t* flag = new size_t[rows];
 fill_n(flag, rows, 0);
 for (size_t i=0; i<rows; i++)
 {
  for (size_t j=rowPoint[i]; j<rowPoint[i+1]; j++ )
  {
   if (values[j]==0)
   {
	flag[i]++;
   }
  }
 }

 accumulate(flag,0,rows);

 nnz -= flag[rows-1];
 size_t* newColIndex = new size_t[nnz];
 T* newValues;
 if (data)
 {
  newValues = new T[nnz];
 }
 size_t counter = 0;
 for (size_t i=0; i<rows; i++)
  {
   for (size_t j=rowPoint[i]; j<rowPoint[i+1]; j++ )
   {
    if (values[j]!=0)
    {
    	 if (data)
    	    	  {
    	    	   newValues[counter]   = values[j];
    	    	  }
    	 newColIndex[counter++] = colIndex[j];
    }
   }
  }

 for (size_t i=0; i<rows; i++)
   {
    rowPoint[i+1] -= flag[i];
   }
  delete[] colIndex;
  colIndex = newColIndex;
  if (data)
  {
   delete[] values;
   values = newValues;
  }
}
;
