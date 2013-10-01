#ifndef __Defines_h_wanghan__
#define __Defines_h_wanghan__

typedef int	Identity;
typedef double	ValueType;

struct IntVectorType
{
  int x, y, z;
  IntVectorType (const int xx=0, const int yy=0, const int zz=0);
};

struct VectorType
{
  ValueType x, y, z;
  VectorType (const ValueType xx=0, const ValueType yy=0, const ValueType zz=0);
};

inline IntVectorType::
IntVectorType (const int xx, const int yy, const int zz)
    : x(xx), y(yy), z(zz)
{
}

inline VectorType::
VectorType (const ValueType xx, const ValueType yy, const ValueType zz)
    : x(xx), y(yy), z(zz)
{
}

#endif
