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

#ifndef _FORWARDS_H
#define _FORWARDS_H

#include <Eigen/src/Core/util/Macros.h>
#ifdef __INTEL_COMPILER
    #define CTR_INLINE inline
#elif __clang__
    #define CTR_INLINE EIGEN_ALWAYS_INLINE
#elif __GNUC__
    #define CTR_INLINE EIGEN_ALWAYS_INLINE
#elif _MSC_VER
    #define CTR_INLINE EIGEN_STRONG_INLINE
#endif

#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x
#define LINE_STRING STRINGIZE(__LINE__)
#define FILE_STRING STRINGIZE(__FILE__)
#ifdef __linux__
    #define FUNC_STRING __PRETTY_FUNCTION__
#elif _WIN32
    #define FUNC_STRING __FUNCSIG__
#endif


#ifdef NDEBUG
    #undef  CT_DEBUG_LEVEL
    #define CT_DEBUG_LEVEL CT_DEBUG_OFF
#endif

#if !defined( CTR_BUILD_SYMBOLS )
    #ifdef WIN32
        #define CTR_IMEXPORT __declspec( dllimport )
    #else
        #define CTR_IMEXPORT
    #endif
#else
    #if defined(WIN32) && !defined( __MINGW32__ )
        #define CTR_IMEXPORT __declspec( dllexport )
    #else
        #define CTR_IMEXPORT
    #endif
#endif

#include <random>
using CT_RNG    = std::mt19937;
typedef std::random_device CT_RD ;
template<class Real>
using CT_RND = std::uniform_real_distribution<Real>;

namespace details{
    template <size_t Size> struct size2type{ typedef int     signed_type; typedef unsigned unsigned_type; };
    template <> struct size2type< 8UL>     { typedef int8_t  signed_type; typedef uint8_t  unsigned_type; };
    template <> struct size2type<16UL>     { typedef int16_t signed_type; typedef uint16_t unsigned_type; };
    template <> struct size2type<32UL>     { typedef int32_t signed_type; typedef uint32_t unsigned_type; };
    template <> struct size2type<64UL>     { typedef int64_t signed_type; typedef uint64_t unsigned_type; };
    template <class RNG , typename... TSeed>
    RNG seed_randomnumbergenerator(TSeed&&... _seeds)
    {

        using std::random_device;
        using std::seed_seq;

        static constexpr size_t ESeedSize = sizeof...(TSeed);
        typedef typename size2type<RNG::word_size>::unsigned_type rng_type;

        std::initializer_list<rng_type> c_SeedIL( {_seeds...} );
        std::array<rng_type, RNG::state_size> c_SeedData;
        std::copy (c_SeedIL.begin(), c_SeedIL.end(), c_SeedData.begin());
        random_device c_RD ;
        std::generate_n(c_SeedData.data()+ESeedSize,
                        c_SeedData.size()-ESeedSize, std::ref(c_RD));

        seed_seq c_seedSeq(std::begin(c_SeedData), std::end(c_SeedData));
        return  RNG(c_seedSeq);
    }
}

#endif // _FORWARDS_H
