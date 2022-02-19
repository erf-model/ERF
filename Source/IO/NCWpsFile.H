#ifndef _NCWPSFILE_H_
#define _NCWPSFILE_H_

template <typename DataType>
struct NDArray
{
 using DType = typename DataType;
 std::vector<size_t> shape;
 DType* data;
}

void ReadNCWpsFile(const std::string  &fname);
#endif
