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

#ifndef CTR_SECTION_H
#define CTR_SECTION_H
#include <initializer_list>
#include <algorithm>
#include <utility>
#include <iostream>

#include <Erl/Utility/util.h>
#include "_forwards.h"
#include "matrix.h"
#include "transrotscale.h"

namespace CTR
{
    namespace sec
    {
        typedef std::integral_constant<size_t,1> CC;
        typedef std::integral_constant<size_t,2> VC;

        namespace details
        {
            template<class TupleT, long NSec, long iSec, long iInc>
            struct TubeCounterR
            {
                 typedef typename std::tuple_element<iSec,TupleT>::type tuple_t;
                 static constexpr size_t ntubes_r =
                             tuple_t::ntube
                           + TubeCounterR<TupleT, (iInc>0)*(iSec+iInc*(1+((NSec-iSec-1)/ iInc)))
                                                 +(iInc<0)*(iSec+iInc*(1+((iSec-NSec-1)/-iInc))),iSec+iInc,iInc>::ntubes_r;

            };
            template<class TupleT, long NSec, long iInc>
            struct TubeCounterR<TupleT,NSec,NSec,iInc>
            {
                 static constexpr size_t ntubes_r = 0;
            };

            template<class TupleT,long iS, long iE, long iInc>
            struct TubeCounter
            {
                static_assert((iInc!= 0) ,"iInc!=0 required.");
                static_assert((iInc < 0) || ( ( (iS<=iE) || (iE==-2) ) && (iE>=-2)  && (iE!=-1) && (iE<=std::tuple_size<TupleT>::value) && (iS>=0) ) ,"When Inc>0,  iS <= iE AND iE<=NSection AND iS>=0 , required");
                static_assert((iInc > 0) || (   (iE<=iS) && (iE>=-1)   && (iS<std::tuple_size<TupleT>::value) )                                      ,"When Inc>0,  iE <= iS AND iE>=0 AND iS<=NSection, required");

                static constexpr size_t ntubes = TubeCounterR<TupleT,iE,iS,iInc>::ntubes_r;
            };

            template<class TupleT, long iS, long iInc>
            struct TubeCounter<TupleT,iS,long(-2), iInc>
            {
                static_assert((iInc!= 0) ,"iInc!=0 required.");
                static_assert((iInc < 0) || ( (iS<std::tuple_size<TupleT>::value) && (iS>=0) ) ,"When Inc>0,  iS <= iE AND iE<NSection AND iS>=0 , required");
                static_assert((iInc > 0) || ( (iS<std::tuple_size<TupleT>::value) && (iS>=0) ) ,"When Inc>0,  iE <= iS AND iE>=0 AND iS<=NSection, required");

                private:
                    static constexpr size_t ntuple = std::tuple_size<TupleT>::value;
                public:
                    static constexpr size_t ntubes = TubeCounterR<TupleT,ntuple,iS,iInc>::ntubes_r;
            };

        }

        template<class TupleT,long iS = 0, long iE = -2, long iInc = 1 >
        struct count_tubes{ static constexpr size_t value =  details::TubeCounter<TupleT,iS,iE,iInc>::ntubes; };

//        template<class TupleT,long iS , long iE, long iInc >
//        struct count_tubes{ static constexpr size_t value =  details::TubeCounter<TupleT,iS,iE,iInc>::ntubes; };

        template<class tuple_type, long isec>
        struct tube_iter{ static constexpr size_t value = count_tubes<tuple_type,0,isec,1>::value;};

    }
    template <class Real, class NTube>
    class  alignas(CTR_ALIGNMENT) Section
    {
        public:
            typedef typename NTube::value_type   ntube_type;
            typedef          Section<Real,NTube> section_type;
            static constexpr ntube_type ntube =  NTube::value ;
            typedef Matrix1N1<Real,ntube,CTR_ALIGNMENT> Mat1T1;

            template<class EigenDerived>
            using Mat = Eigen::MatrixBase<EigenDerived>;




            CTR_INLINE explicit  Section(){}

            template<class IL=std::initializer_list<Real> >
            CTR_INLINE explicit  Section(IL&& _Stiff, IL&& _Curv, Real&& _Length,
                              IL&& _Radius, IL&& _PoissRatio)
                : m_Stiffness()
                , m_Curvature{ (std::copy(_Stiff.begin(),_Stiff.end(),m_Stiffness.data()), Mat1T1::Zero()) }
                , m_Length   ( (std::copy(_Curv .begin(),_Curv .end(),m_Curvature.data()),_Length) )
                , m_Radius         ()
                , m_PoissonRatio   {(std::copy(_Radius.begin(),_Radius.end(),m_Radius.data()) , Mat1T1::Zero() )}
                , m_CurvedArcStart ((std::copy(_PoissRatio.begin(),_PoissRatio.end(),m_PoissonRatio.data()), Real(0) ))
                , m_CurvedArcEnd   (0)
                , m_Phi            (0)
                , m_AlphaTip       (Mat1T1::Zero())
            {}
            template<class IL=std::initializer_list<Real> >
            CTR_INLINE explicit  Section(const IL& _Stiff, const IL& _Curv, const Real& _Length,
                              const IL& _Radius, const IL& _PoissRatio)
                : m_Stiffness()
                , m_Curvature{ (std::copy(_Stiff.begin(),_Stiff.end(),m_Stiffness.data()), Mat1T1::Zero()) }
                , m_Length   ( (std::copy(_Curv .begin(),_Curv .end(),m_Curvature.data()),_Length) )
                , m_Radius         ()
                , m_PoissonRatio   {(std::copy(_Radius.begin(),_Radius.end(),m_Radius.data()) , Mat1T1::Zero() )}
                , m_CurvedArcStart ((std::copy(_PoissRatio.begin(),_PoissRatio.end(),m_PoissonRatio.data()), Real(0) ))
                , m_CurvedArcEnd   (0)
                , m_Phi            (0)
                , m_AlphaTip       (Mat1T1::Zero())
            {}
            template<class IL=std::initializer_list<Real> >
            CTR_INLINE explicit  Section(IL&& _Stiff, IL&& _Curv, Real&& _Length,
                              IL&& _Radius, IL&& _PoissRatio,
                              Real&& _CurvArcS, Real&& _CurvArcE, Real&& _Phi,IL&& _AlphaTip)
                : m_Stiffness()
                , m_Curvature{ (std::copy(_Stiff.begin(),_Stiff.end(),m_Stiffness), Real(0) ) }
                , m_Length   ( (std::copy(_Curv .begin(),_Curv .end(),m_Curvature),_Length) )
                , m_Radius         ()
                , m_PoissonRatio   {(std::copy(_Radius.begin(),_Radius.end(),m_Radius) , Real(0) )}
                , m_CurvedArcStart ((std::copy(_PoissRatio.begin(),_PoissRatio.end(),m_PoissonRatio), _CurvArcS ))
                , m_CurvedArcEnd   (_CurvArcE)
                , m_Phi            (_Phi)
                , m_AlphaTip       ()
            { std::copy(_AlphaTip.begin(),_AlphaTip.end(),m_AlphaTip); }
            CTR_INLINE ~Section(){}


            // ********* GET FUNCTIONS *****************************************
            CTR_INLINE Real getOuterRadius() const { return m_Radius[0      ];}
            CTR_INLINE Real getInnerRadius() const { return m_Radius[ntube-1];}
            template<ntube_type i>
            CTR_INLINE Real getRadius()                    const { static_assert(i<ntube ,  "wrong index"); return m_Radius[i]; }
            CTR_INLINE Real getRadius(const ntube_type& i) const {        assert(i<ntube && "wrong index"); return m_Radius[i]; }
            template<ntube_type i>
            CTR_INLINE Real getPossionRatio()                    const { static_assert(i<ntube ,  "wrong index"); return m_PoissonRatio[i]; }
            CTR_INLINE Real getPossionRatio(const ntube_type& i) const {        assert(i<ntube && "wrong index"); return m_PoissonRatio[i]; }
            CTR_INLINE Real getPhi()                             const {                                          return m_Phi;             }
            template<ntube_type i>
            CTR_INLINE Real getTubeStiffness()                    const { static_assert(i<ntube ,  "wrong index"); return m_Stiffness[i]; }
            CTR_INLINE Real getTubeStiffness(const ntube_type& i) const {        assert(i<ntube && "wrong index"); return m_Stiffness[i]; }
            template<ntube_type i>
            CTR_INLINE Real getTubeCurvature()                    const { static_assert(i<ntube ,  "wrong index"); return m_Curvature[i]; }
            CTR_INLINE Real getTubeCurvature(const ntube_type& i) const {        assert(i<ntube && "wrong index"); return m_Curvature[i]; }
            template<ntube_type i>
            CTR_INLINE Real getTubeAlphaTip()                    const { static_assert(i<ntube ,  "wrong index"); return m_AlphaTip[i]; }
            CTR_INLINE Real getTubeAlphaTip(const ntube_type& i) const {        assert(i<ntube && "wrong index"); return m_AlphaTip[i]; }

            CTR_INLINE Real getCurvedArcStart()                  const {                                          return m_CurvedArcStart; }
            CTR_INLINE Real getCurvedArcEnd()                    const {                                          return m_CurvedArcEnd;   }
            CTR_INLINE Real getLength()                          const {                                          return m_Length;         }

            // ********* SET FUNCTIONS *****************************************
            template<ntube_type i_tubn>
            CTR_INLINE void print_something() {std::cout<<i_tubn<<std::endl;}

            static CTR_INLINE std::ostream& print_type(std::ostream& os)
            {
                os << Erl::static_if_else<(ntube==sec::CC::value)>()(
                    [](){ return "CC"; },
                    [](){ return "VC"; } );
                return os;
            }

            // ***********************************************************************************************************
            static CTR_INLINE void getCurvatureXY(Real _stiff,Real _curvature, Real _sinAlpha, Real _cosAlpha, Real& _curvatureX, Real& _curvatureY)
            {
                _curvatureX += -_sinAlpha*_curvature*(_stiff/ntube);
                _curvatureY +=  _cosAlpha*_curvature*(_stiff/ntube);
            }
            static CTR_INLINE void getTorsionalCurvatureXY(Real _sinAlpha, Real _cosAlpha,
                                         Real _cumm_curvatureX, Real _cumm_curvatureY, Real _cumm_stiffness,
                                         Real & _torsional_curvatureX, Real & _torsional_curvatureY)
            {
                _torsional_curvatureX = ( _cosAlpha*_cumm_curvatureX + _sinAlpha*_cumm_curvatureY)/_cumm_stiffness;
                _torsional_curvatureY = (-_sinAlpha*_cumm_curvatureX + _cosAlpha*_cumm_curvatureY)/_cumm_stiffness;
            }
            // ***********************************************************************************************************

            template<ntube_type i_secn,ntube_type i_tubn, class EigenDerived>
            CTR_INLINE void setJointsCentreLineStates(const Mat<EigenDerived>& _JValIn,Real * _AlphaTip, Real & _PrevArcEnd)
            {
                m_Phi       = _JValIn[1+i_secn+i_tubn];
                assert( (m_Phi>=0) && (m_Phi<= m_Length) );
                Erl::static_for<0,ntube>()([&](auto&& i_ic)
                {
                    constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                    _AlphaTip[i_tubn+itub] = m_AlphaTip[itub]  = _JValIn[2+i_secn+i_tubn+itub];
                });
                m_CurvedArcEnd   = (_PrevArcEnd+=m_Phi);
                m_CurvedArcStart = std::max<Real>(Real(0), m_CurvedArcEnd- m_Length);
            }

            template<ntube_type i_tubn, class EigenDerived>
            CTR_INLINE void addCurvatureAndStiffness(const Mat <EigenDerived>& _sinAlpha, const Mat <EigenDerived>& _cosAlpha,
                                                 const Real& _stiffness, const Real& _curvature,
                                                 Real& _cumm_curvatureX, Real& _cumm_curvatureY, Real& _cumm_stiffness) const
            {
                _cumm_stiffness += _stiffness;
                Erl::static_for<0,ntube>()([&](auto&& i_ic)
                {
                    constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                    CTR::Section<Real,NTube>::getCurvatureXY(_stiffness,_curvature,_sinAlpha[i_tubn+itub],_cosAlpha[i_tubn+itub],_cumm_curvatureX,_cumm_curvatureY);
                });
            }

            template<ntube_type i_tubn, class EigenDerived>
            CTR_INLINE void calcTorsionalCurvatureXY(const Mat <EigenDerived>& _sinAlpha, const Mat <EigenDerived>& _cosAlpha,
                                    const Real& _cumm_curvatureX, const Real& _cumm_curvatureY, const Real& _cumm_stiffness,
                                    Real * _torsional_curvatureXYZ) const
            {
                Erl::static_for<0,ntube>()([&](auto&& i_ic)
                {
                   constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                   CTR::Section<Real,NTube>::getTorsionalCurvatureXY(_sinAlpha[itub+i_tubn],_cosAlpha[itub+i_tubn],
                                                                     _cumm_curvatureX,_cumm_curvatureY,
                                                                     _cumm_stiffness,
                                                                     _torsional_curvatureXYZ[3*(itub+i_tubn)],
                                                                     _torsional_curvatureXYZ[3*(itub+i_tubn)+1]);
                });
            }
            template<ntube_type i_tubn>
            CTR_INLINE void getTorsionalCurvatureZ(Real _ArcLengthStep, Real _stiffness, Real _curvature,
                                        const Real & _torsional_curvatureX,
                                        Real & _torsional_curvatureZ, Real &_cumm_torsional_curvatureZ) const
            {
                _cumm_torsional_curvatureZ += (_stiffness/ntube)/(Real(1)+m_PoissonRatio[i_tubn])*_torsional_curvatureZ;
                _torsional_curvatureZ      += _ArcLengthStep*(Real(1)+m_PoissonRatio[i_tubn])*_torsional_curvatureX*_curvature;
            }


            template<ntube_type i_tubn>
            CTR_INLINE void calcTorsionalCurvatureZ(Real _ArcLengthStep, Real   _stiffness, Real _curvature, Real* _torsional_curvatureXYZ, Real&  _cumm_torsional_curvatureZ) const
            {
                Erl::static_for<ntube-1,-1,-1>()([&](auto && i_ic)
                {
                    constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                    this->getTorsionalCurvatureZ<itub>(_ArcLengthStep,_stiffness,_curvature,
                                               _torsional_curvatureXYZ[3*(i_tubn+itub)  ],
                                               _torsional_curvatureXYZ[3*(i_tubn+itub)+2],
                                                _cumm_torsional_curvatureZ);
                });
            }
            template<ntube_type i_tubn>
            CTR_INLINE void calcTorsionalCurvatureZLastSection(Real  _ArcLengthStep,
                                                    Real  _stiffness, Real _curvature,
                                                    Real* _torsional_curvatureXYZ    ,
                                                    Real& _cumm_torsional_curvatureZ ) const
            {
                Erl::static_if<(ntube==sec::VC::value)>()([&]()
                {
                    this->getTorsionalCurvatureZ<ntube-1>(_ArcLengthStep,_stiffness,_curvature,
                                                    _torsional_curvatureXYZ[3*(i_tubn)  ],
                                                    _torsional_curvatureXYZ[3*(i_tubn)+2],
                                                    _cumm_torsional_curvatureZ);
                });
            }

            template<ntube_type i_tubn>
            CTR_INLINE void calcTorsionalCurvatureZLastTube( Real   _stiffness,
                                                         Real*  _torsional_curvatureXYZ,
                                                         const Real&  _cumm_torsional_curvatureZ) const
            {
                    _torsional_curvatureXYZ[3*(i_tubn+ntube-1)+2] = (-_cumm_torsional_curvatureZ)*(Real(1)+m_PoissonRatio[ntube-1])/_stiffness;
            }


            template<ntube_type i_tubn>
            CTR_INLINE void calcAlpha(Real _ArcLengthStep,
                           Real * _sample_alpha,
                           Real * _torsional_curvatureXYZ) const
            {
                Erl::static_for<ntube-1,-1,-1>()([&](auto && i_ic)
                {                   
                    constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                    _sample_alpha[itub+i_tubn] = _sample_alpha[itub+i_tubn]+_ArcLengthStep*(_torsional_curvatureXYZ[3*(itub+i_tubn)+2]-_torsional_curvatureXYZ[2]);
                });

            }

            template<ntube_type i_tubn>
            CTR_INLINE void calcAlphaFirstSection(Real _ArcLengthStep,
                                       Real * _sample_alpha,
                                       Real * _torsional_curvatureXYZ) const
            {

                _sample_alpha[i_tubn]=Real(0);
                Erl::static_if<(ntube==sec::VC::value)>()([&](){
                    _sample_alpha[i_tubn+1] = _sample_alpha[i_tubn+1]+_ArcLengthStep*(_torsional_curvatureXYZ[3*(i_tubn+1)+2]-_torsional_curvatureXYZ[2]);
                });
            }

            CTR_INLINE Real getStiffness(const Real& _SamplePoint) const
            {
                return ( (_SamplePoint<=m_CurvedArcEnd) * m_Stiffness.static_sum() );
            }
            CTR_INLINE Real getCurvature(const Real& _SamplePoint) const
            {
                return (((_SamplePoint<=m_CurvedArcEnd) && (_SamplePoint>=m_CurvedArcStart))*m_Curvature[0]);
            }
            CTR_INLINE void getCurvatureStiffness(const Real& _SamplePoint,Real& _Curvature,Real& _Stiffness) const
            {
              _Curvature =    (((_SamplePoint<=m_CurvedArcEnd) && (_SamplePoint>=m_CurvedArcStart))*m_Curvature[0]);
              _Stiffness =    (_SamplePoint<=m_CurvedArcEnd) * (m_Stiffness.static_sum());
            }

            template<ntube_type i_secn,ntube_type i_tubn>
            CTR_INLINE void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Real*_JValReturn) const
            {
                _JValReturn[1+i_secn+i_tubn]= m_Length*_RND(_RNG);
                Erl::static_for<0,ntube>()([&](auto && i_ic)
                {
                     constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                    _JValReturn[2+i_secn+i_tubn+itub]= std::nextafter(Real(2*M_PI),Real(0))*_RND(_RNG);
                });
            }
            template<ntube_type i_secn,ntube_type i_tubn>
            CTR_INLINE void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Real*_JValReturn, const TransRotScale<Real>& _TS) const
            {
                _JValReturn[1+i_secn+i_tubn]= m_Length*(_TS.scaleTrans(_RND(_RNG)));;
                Erl::static_for<0,ntube>()([&](auto && i_ic)
                {
                     constexpr auto itub = std::decay<decltype(i_ic)>::type::value ;
                    _JValReturn[2+i_secn+i_tubn+itub]= std::nextafter(Real(2*M_PI),Real(0))*_RND(_RNG);
                });
            }


            friend  std::ostream&  operator<< (std::ostream& os, const Section &obj)
            {
                os<<"NTube      : "<<obj.ntube                     <<std::endl;
                os<<"Stiffness  : "<<obj.m_Stiffness   .transpose()<<std::endl;
                os<<"Curvature  : "<<obj.m_Curvature   .transpose()<<std::endl;
                os<<"Length     : "<<obj.m_Length                  <<std::endl;
                os<<"Radius     : "<<obj.m_Radius      .transpose()<<std::endl;
                os<<"PoissionRat: "<<obj.m_PoissonRatio.transpose()<<std::endl;
                os<<"CurvAS     : "<<obj.m_CurvedArcStart          <<std::endl;
                os<<"CurvAE     : "<<obj.m_CurvedArcEnd            <<std::endl;
                os<<"Phi        : "<<obj.m_Phi                     <<std::endl;
                os<<"AlphaTip   : "<<obj.m_AlphaTip    .transpose()<<std::flush;
                return os;
            }
        protected:
            alignas(CTR_ALIGNMENT) Mat1T1 m_Stiffness    ;
            alignas(CTR_ALIGNMENT) Mat1T1 m_Curvature    ;
            alignas(CTR_ALIGNMENT) Real   m_Length         ;
            alignas(CTR_ALIGNMENT) Mat1T1 m_Radius       ;
            alignas(CTR_ALIGNMENT) Mat1T1 m_PoissonRatio ;
            alignas(CTR_ALIGNMENT) Real   m_CurvedArcStart ;
            alignas(CTR_ALIGNMENT) Real   m_CurvedArcEnd   ;
            alignas(CTR_ALIGNMENT) Real   m_Phi            ;
            alignas(CTR_ALIGNMENT) Mat1T1 m_AlphaTip     ;
    };

    template<class Real>
    using SectionCC = Section<Real,sec::CC>;
    template<class Real>
    using SectionVC = Section<Real,sec::VC>;
}
#endif // CTR_SECTION_H
