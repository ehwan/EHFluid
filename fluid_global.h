#pragma once

#include "../EHMatrix/EHMatrix.h"

typedef EH::Matrix::Matrix< float , 2 , 1 > vector2;
template < typename T >
using vec2 = EH::Matrix::Matrix< T , 2 , 1 >;
