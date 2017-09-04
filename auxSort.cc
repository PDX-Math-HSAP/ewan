template<class T> //T must have a > method.
size_t* rankSort(T* dataArray,  size_t startIndex, size_t endIndex, string direction)
{
  size_t    arrSize  = endIndex - startIndex;
  size_t* indexArray = new size_t[arrSize];//construct array of indices
  for(size_t k=0; k < arrSize; k++)
  {
    indexArray[k] = k+startIndex;
  }
   
  if (direction == "descend")
  {
    sort(&indexArray[0],&indexArray[arrSize],[&dataArray] (size_t a, size_t b)
      {
        return dataArray[a]>dataArray[b];
      });
  }
  else
  {
    sort(&indexArray[0],&indexArray[arrSize],[&dataArray] (size_t a, size_t b)
      {
        return dataArray[a]<dataArray[b];
      });
  }
  return indexArray;
}

template<class T> //T must have a > method.
void permute(T* array, size_t* permutation, size_t startIndex, size_t endIndex)
{
  size_t size = endIndex - startIndex;
  auto temp = new T[size]; 
  for (size_t j=0; j<size ; j++)
  {
    temp[j] = array[j+startIndex];
    if (permutation[j]<j+startIndex)
    {
      array[j+startIndex] = temp[permutation[j]-startIndex];
    }
    else
    {
      array[j+startIndex] = array[permutation[j]];
    }
    
  }
  delete[]temp;
}


void accumulate(size_t* array, size_t start, size_t end) //Sets array[i] = array[start]+...+array[i] for i in [start,end).
{
  for (size_t i=start; i<end; i++)
  {
    array[i+1] += array[i];
  }
};
