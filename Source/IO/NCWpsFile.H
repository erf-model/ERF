#ifndef _NCWPSFILE_H_
#define _NCWPSFILE_H_

template <typename DataType>
struct NDArray
{
 using DType = typename std::remove_const<DataType>::type;

 template <typename T>
 explicit NDArray (const std::string vname, const std::vector<size_t>& vshape) : name{vname}, shape{vshape} {
     data = new T [shape.size()]; 
 } 

 template <typename T>
 T* get_data() {
   return static_cast<T*>(data);
 }
 
 std::string name;
 std::vector<size_t> shape;
 DType* data;
};

void ReadNCWpsFile(const std::string  &fname);
#endif
