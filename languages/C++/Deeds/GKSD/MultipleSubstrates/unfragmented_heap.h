/* Templated class to implement an "unfragmented heap"
   Based on the suggestion of Jean Krivine!
*/

template <class T> class unfragmented_heap
{
 public:
  unfragmented_heap() ; 
  unfragmented_heap(int) ;

  T& at(int) ;
  T& operator[] (int n) ;
  void add(T) ;
  void erase(int) ;
  int size() ;
  int tot_size() ;
  void gc() ;
 protected:
  int s ;
  vector <T> v ;
};

template <class T> unfragmented_heap<T>::unfragmented_heap()
{
  s = 0 ;
  v.resize(0) ;
}

template <class T> unfragmented_heap<T>::unfragmented_heap(int i)
{
  s = 0 ;
  v.resize(i) ;
}

template <class T> T& unfragmented_heap<T>::operator[](int i)
{
  return v[i] ;
}


template <class T> T& unfragmented_heap<T>::at(int i)
{
  if (i >= s)
    {
      cout << "Attempted to access an element out-of-range in unfragmented heap" << endl ;
      exit(1) ;
    }
  return v[i] ;
}

template <class T> void unfragmented_heap<T>::add(T x)
{
  if (s + 1 > v.size())
    {
      v.resize(2*s +1) ;
      v[s] = x ;
      s ++ ;
    }
  else
    {
      v[s] = x ;
      s ++ ;
    }
}

template <class T> void unfragmented_heap<T>::erase(int i)
{
  v[i] = v[s -1] ;
  s -- ;
}

template <class T> int unfragmented_heap<T>::size()
{
  return s ;
}

template <class T> int unfragmented_heap<T>::tot_size()
{
  return v.size() ;
}

template <class T> void unfragmented_heap<T>::gc()
{
  if (4*s < v.size())
    {
      v.resize(4*s) ;
    }
}
