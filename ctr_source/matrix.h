// This file is part of CTR, a kinematics library for concentric tube robots
//
// Copyright (C) 2017 Konrad Leibrandt <konrad.lei@gmx.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef CTMATRIX_H
#define CTMATRIX_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Erl/Utility/static.h>
#include "_forwards.h"

#define CT_ASSERT_ALIGNMENT_THIS         assert(checkAlignment(this->data())  && "Alignment of this  not satisfied");
#define CT_ASSERT_CORRECT_R_C_PARAMETER  assert( ( (_Rows==Dynamic) || (_R==_Rows) ) && \
                                                 ( (_Cols==Dynamic) || (_C==_Cols) ) && \
                                                 ( (_R   >= 0      ) && (_C>= 0)   ) && \
                                                 "Wrong _R or _C parameter for this class");
#define CT_ASSERT_SIZE_AT_COMPILE        assert( (_Rows>0) && (_Cols>0) && \
                                                "Size must be known at compile time" );

namespace  CTR {
using   Eigen::Dynamic;
using   Eigen::ColMajor;
using   Eigen::RowMajor;



namespace details
{
    template<typename Matrix, typename Real, int _Rows, int _Cols, bool _MemAllocFree>
    struct MatrixSetSize
    {
        typedef Eigen::DenseIndex DenseIndex;
        static void setSize(Matrix * _thisptr, DenseIndex _R, DenseIndex _C)
        {
            Erl::static_if<(_Rows==Dynamic)>()([&_R,_thisptr](){
                memcpy(reinterpret_cast<void*>(std::ptrdiff_t(_thisptr)+sizeof(Real*)),
                       &_R,sizeof(DenseIndex));
            });
            Erl::static_if<(_Cols==Dynamic)>()([&_C,_thisptr](){
                Erl::static_if_else<_Rows==Dynamic>()
                (
                    [&_C,_thisptr]()
                    {
                        memcpy(reinterpret_cast<void*>(std::ptrdiff_t(_thisptr)+sizeof(Real*)+sizeof(DenseIndex)),
                               &_C,sizeof(DenseIndex));
                    },
                    [&_C,_thisptr]()
                    {
                        memcpy(reinterpret_cast<void*>(std::ptrdiff_t(_thisptr)+sizeof(Real*)),
                               &_C,sizeof(DenseIndex));
                    }
                );
            });
        }
    };
    template<typename Matrix, typename Real, int _Rows, int _Cols>
    struct MatrixSetSize<Matrix, Real,_Rows,_Cols,false>
    {
        typedef Eigen::DenseIndex DenseIndex;
        static void setSize(Matrix * _thisptr, DenseIndex _R, DenseIndex _C)
        {
            _thisptr->resize(_R,_C);
        }
    };
}

template<typename Real, int _Rows, int _Cols, size_t _Align=16, int _Options = ColMajor, bool _MemAllocFree = true>
class
alignas(
        (  (_Rows>0)*(_Cols>0) *(alignof(Eigen::Matrix<Real,_Rows, _Cols, _Options>) >_Align)*alignof(Eigen::Matrix<Real,_Rows, _Cols, _Options>))+
        (  (_Rows>0)*(_Cols>0) *(alignof(Eigen::Matrix<Real,_Rows, _Cols, _Options>)<=_Align)*_Align)+
        (!((_Rows>0)*(_Cols>0))*_Align)
       )
//alignas(_Align)
Matrix : public Eigen::Matrix<Real,_Rows, _Cols, _Options>
{

public:
    typedef Eigen::Matrix<Real,_Rows, _Cols, _Options> Base;
    typedef Eigen::DenseIndex                          DenseIndex;
    static constexpr size_t align_value =
           (  (_Rows>0)*(_Cols>0) *(alignof(Eigen::Matrix<Real,_Rows, _Cols, _Options>) >_Align)*alignof(Eigen::Matrix<Real,_Rows, _Cols, _Options>))+
           (  (_Rows>0)*(_Cols>0) *(alignof(Eigen::Matrix<Real,_Rows, _Cols, _Options>)<=_Align)*_Align)+
           (!((_Rows>0)*(_Cols>0))*_Align);
//            _Align;

    CTR_INLINE Matrix()
        : Base()
    {
        CT_ASSERT_ALIGNMENT_THIS
    }

    CTR_INLINE Matrix(const Matrix& _o)
        : Base(_o)
    {
        CT_ASSERT_ALIGNMENT_THIS
    }
    template<typename _T>
    CTR_INLINE Matrix(const Eigen::MatrixBase<_T>& _o)
        : Base(_o)
    {
        CT_ASSERT_ALIGNMENT_THIS
    }
    CTR_INLINE ~Matrix()
    {
        CT_ASSERT_ALIGNMENT_THIS
        setDataPtr(nullptr);
    }

    CTR_INLINE Matrix(int _R, int _C)
        : Base()
    {
        CT_ASSERT_ALIGNMENT_THIS
        CT_ASSERT_CORRECT_R_C_PARAMETER
        setSize(_R,_C);
    }
    CTR_INLINE Matrix(int _R, int _C, void * _Data)
        : Base()
    {
        CT_ASSERT_ALIGNMENT_THIS
        CT_ASSERT_CORRECT_R_C_PARAMETER
        setDataPtr(_Data);
        setSize(_R,_C);
    }
    CTR_INLINE void setDataPtr(void * _DataPtr)
    {
        Erl::static_if_else<_MemAllocFree>()
        (
            [&_DataPtr,this]()
            {
                Erl::static_if<( (_Rows==Dynamic) || (_Cols==Dynamic) )>()([&_DataPtr,this](){
                    CT_ASSERT_ALIGNMENT_THIS
                    memcpy(this,&_DataPtr,sizeof(_DataPtr));
                });
            },
            []()
            {
                // Do Nothing
            }
        );
    }
    CTR_INLINE void setDataPtr(void * _DataPtr,DenseIndex _R, DenseIndex _C)
    {
        Erl::static_if_else<_MemAllocFree>()
        (
            [&_DataPtr,this,&_R,&_C]()
            {
                Erl::static_if<( (_Rows==Dynamic) || (_Cols==Dynamic) )>()([&_DataPtr,this](){
                    CT_ASSERT_ALIGNMENT_THIS
                    memcpy(this,&_DataPtr,sizeof(_DataPtr));
                });
                setSize(_R,_C);
            },
            []()
            {
                // Do Nothing
            }
        );
    }
    CTR_INLINE void setSize(DenseIndex _R, DenseIndex _C)
    {
        details::MatrixSetSize<Matrix,Real,_Rows,_Cols,_MemAllocFree>::setSize(this,_R,_C);
    }
    CTR_INLINE Real static_sum() const
    {
        CT_ASSERT_SIZE_AT_COMPILE;
        Real _sum = Base::data()[0];
        Erl::static_for<1,_Rows*_Cols>()([&](auto&& i_ic)
        {
            constexpr auto iel = std::decay<decltype(i_ic)>::type::value;
            _sum += this->Base::data()[iel];
        });
        return _sum;
    }
    CTR_INLINE static size_t requiredSize(int _R, int _C)
    {
        CT_ASSERT_CORRECT_R_C_PARAMETER
        return
        Erl::static_if_else<( (_Rows==Dynamic) || (_Cols==Dynamic) )>()
        (
            [&_R,&_C]()
            {
                return sizeof(Matrix) + sizeof(Real)*(size_t(_R)*size_t(_C));

            },
            []()
            {
                return sizeof(Matrix);
            }

        );
    }
    CTR_INLINE static size_t requiredPtrMemSize(int _R, int _C)
    {
        CT_ASSERT_CORRECT_R_C_PARAMETER
        return
        Erl::static_if_else<( (_Rows==Dynamic) || (_Cols==Dynamic) )>()
        (
            [&_R,&_C]()
            {
                return sizeof(Real)*(size_t(_R)*size_t(_C));

            },
            []()
            {
                return 0;
            }

        );
    }
    CTR_INLINE static size_t requiredSizePadded(int _R, int _C, int& _RowsPadded)
    {
        assert(_Options == ColMajor && "This function gives only valid values for ColumnMajor Matrix");
        assert(_Rows    == Dynamic  && "This function gives only valid values if Rows are Dynamic");

        size_t c_PaddedRBytes = ( (sizeof(Real)*_R) +(align_value-1)) & ~(align_value-1);
        assert( (c_PaddedRBytes%sizeof(Real) == 0) && "Something went wrong");
        _RowsPadded = c_PaddedRBytes/sizeof(Real);
        return requiredSize(_RowsPadded,_C);
    }
    CTR_INLINE static size_t requiredSizePadded(int _R, int _C)
    {
        int c_RowsPadded;
        return requiredSizePadded(_R,_C,c_RowsPadded);
    }

    CTR_INLINE static size_t requiredPtrMemSizePadded(int _R, int _C, int& _RowsPadded)
    {
        assert(_Options == ColMajor && "This function gives only valid values for ColumnMajor Matrix");
        assert(_Rows    == Dynamic  && "This function gives only valid values if Rows are Dynamic");

        size_t c_PaddedRBytes = ( (sizeof(Real)*_R) +(align_value-1)) & ~(align_value-1);
        assert( (c_PaddedRBytes%sizeof(Real) == 0) && "Something went wrong");
        _RowsPadded = c_PaddedRBytes/sizeof(Real);
        return requiredPtrMemSize(_RowsPadded,_C);
    }
    CTR_INLINE static size_t requiredPtrMemSizePadded(int _R, int _C)
    {
        int c_RowsPadded;
        return requiredPtrMemSizePadded(_R,_C,c_RowsPadded);
    }

    CTR_INLINE static bool checkAlignment(void * _Ptr)
    {
        return ((std::uintptr_t(_Ptr) & (align_value-1)) == 0);
    }
};

template<class Real,int N, int Align> using Matrix3N1 = Matrix<Real,           3*N,1,Align>;
template<class Real,int N, int Align> using Matrix1N1 = Matrix<Real,             N,1,Align>;
template<class Real,       int Align> using Matrix1X3 = Matrix<Real,Eigen::Dynamic,3,Align>;
template<class Real,       int Align> using Matrix31X = Matrix<Real,3,Eigen::Dynamic,Align>;
template<class Real,       int Align> using Matrix1X1 = Matrix<Real,Eigen::Dynamic,1,Align>;

}

#endif // CTMATRIX_H
