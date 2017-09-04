//Vec class method definitions
template<class T>
Vec<T>::Vec(size_t dim)
{
  size = dim;
  talk = false;
  entry = new T[dim];
  fill_n(entry,size,0);
  if (talk)
  {
    cout << "Allocating memory for size " << size <<" Vec.\n"; 
  }
};
template<class T>
void Vec<T>::destroy()
{
  delete[]entry;
  if (talk)
  {
    cout << "Freed memory for " << size <<" by 1 array.\n"; 
  };
};
template<class T>
Vec<T>::~Vec()
{
  //destroy();
};
template<class T>
void Vec<T>::accumulate() //Sets rowPoint[i] = rowPoint[0]+...+rowPoint[i]. O(number of rows)
{
  for (size_t i=0; i<size-1; i++)
  {
    entry[i+1] += entry[i];
  }
};
template<class T>
void Vec<T>::print(const string name)
{
  cout<< name << "= [";
  for (size_t j = 0; j < size; j++ ) 
  {
    if (j == size-1) 
    {
      cout << entry[j]<< "]\n";
    }
    else 
    {
      cout << entry[j] << ",";
    }
  }
};

