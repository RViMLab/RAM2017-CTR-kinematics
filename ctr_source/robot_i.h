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

#ifndef CTR_ROBOT_I_H
#define CTR_ROBOT_I_H

#include <tuple>
#include <Erl/Core/transform.h>
#include <Erl/Core/rotmat.h>
#include <Erl/Core/vector3.h>
#include <Erl/Core/vector6.h>

#include <Erl/Utility/util.h>
#include "_forwards.h"
#include "section_i.h"
#include "matrix.h"

#include "transrotscale.h"
#include <Erl/Utility/singleton.h>

#ifdef WITH_CTROPENCL
#include "robot_cl.h"
#endif

namespace CTR
{
    template<class Real>
    void JacobianColumnFiniteDifference(const Erl::Transform<Real> & _Tq_lo ,
                                        const Erl::Transform<Real> & _Tq_hi ,
                                        const Real& _dq,
                                        Real* _JacobianColumn)
    {
        /* */
        Erl::Transform<Real> c_Tq_dq = _Tq_hi - _Tq_lo;
        c_Tq_dq.setRotation(c_Tq_dq.getRotationReference()*_Tq_hi.getRotationReference().transpose());
        _JacobianColumn[0] = c_Tq_dq.getTranslationReference().x()/_dq;
        _JacobianColumn[1] = c_Tq_dq.getTranslationReference().y()/_dq;
        _JacobianColumn[2] = c_Tq_dq.getTranslationReference().z()/_dq;
        _JacobianColumn[3] = ((c_Tq_dq(2,1)-c_Tq_dq(1,2))*Real(.5))/_dq ;
        _JacobianColumn[4] = ((c_Tq_dq(0,2)-c_Tq_dq(2,0))*Real(.5))/_dq ;
        _JacobianColumn[5] = ((c_Tq_dq(1,0)-c_Tq_dq(0,1))*Real(.5))/_dq ;
        /* */

        /*
        _JacobianColumn[0] = ( _Tq_hi.getTranslationReference().x()
                              -_Tq_lo.getTranslationReference().x())/_dq;
        _JacobianColumn[1] = ( _Tq_hi.getTranslationReference().y()
                              -_Tq_lo.getTranslationReference().y())/_dq;
        _JacobianColumn[2] = ( _Tq_hi.getTranslationReference().z()
                              -_Tq_lo.getTranslationReference().z())/_dq;

        Erl::Vector3<Real> c_Rq_dq =  Erl::Quaternion<Real>(_Tq_hi.getRotationReference()*
                                                            _Tq_lo.getRotationReference().inv())
                                      .getRotationVector();

        _JacobianColumn[3] = c_Rq_dq.x()/_dq;
        _JacobianColumn[4] = c_Rq_dq.y()/_dq;
        _JacobianColumn[5] = c_Rq_dq.z()/_dq;
        */
    }

    namespace details
    {
        template<bool geq, long iSec>
        struct TubNGEqTubMax
        { static constexpr auto nsec = iSec; };
        template<long iSec>
        struct TubNGEqTubMax<false,iSec>
        { static constexpr auto nsec = 0;    };

        template<class TupleT, long iTube, long iTubeF>
        struct SectionCounterR
        {
            static constexpr auto nsec   = std::tuple_size<TupleT>::value;
            static constexpr auto ntub = CTR::sec::tube_iter<TupleT, iTube>::value;
            static constexpr auto isec = TubNGEqTubMax<ntub>=iTubeF,iTube> ::nsec
                                       + (ntub<iTubeF)*
                                         SectionCounterR<TupleT
                                                        ,iTube+1
                                                        ,(ntub< iTubeF)*iTubeF+
                                                         (ntub>=iTubeF)*(iTube+1)>::isec;
        };
        template<class TupleT, long iTube>
        struct SectionCounterR<TupleT,iTube,iTube>
        {
             static constexpr auto nsec   = std::tuple_size<TupleT>::value;
             static constexpr auto ntub   = CTR::sec::tube_iter<TupleT, (iTube<=nsec)*iTube>::value;
             static constexpr size_t isec = TubNGEqTubMax<ntub>=iTube,iTube> ::nsec;
        };
    }
    template<class TupleT, long iTube>
    struct SectionCounter
    { static constexpr size_t isec = details::SectionCounterR<TupleT,0,iTube>::isec; };

    template<class TupleT, long i_alpha>
    struct JointValue_ALPHA
    {
        static_assert(i_alpha< CTR::sec::count_tubes<TupleT>::value, "not fullfilled: i_phi <  number of tubes");
        static_assert(i_alpha>=0                                   , "not fullfilled: i_phi >= 0");
        static constexpr size_t ijv = SectionCounter<TupleT,i_alpha+1>::isec
                                      +i_alpha+1;

    };
    template<class TupleT, long i_phi>
    struct JointValue_PHI
    {
        static_assert(i_phi< std::tuple_size<TupleT>::value, "not fullfilled: i_phi <  number of sections");
        static_assert(i_phi>=0                              ,"not fullfilled: i_phi >= 0");
        static constexpr size_t ijv = CTR::sec::tube_iter<TupleT, i_phi>::value+
                                    i_phi+1;
    };

}


namespace CTR
{
    template<class Real, class PSecTup, class ...Args>
    class RobotI ;

    template<class Real, class PSecTup, class ...Args>
    class RobotI_ ;

    template<class Real, class PSecTup>
    class RobotBase;

    namespace details {
        template<long NumSection>
        struct append_helper
        {
            template<class Real, class SecApp, class PSecTup, class ...Args >
            CTR_INLINE static RobotBase<Real,PSecTup>* appendSec(const RobotI<Real,PSecTup,Args...>& _RobPrev, const SecApp& _Sapp) noexcept
            { //return Erl::alloc_aligned_init_array<RobotI<Real,PSecTup,Args...,SecApp> >(size_t(1),RobotI<Real,PSecTup,Args...,SecApp>::align_value,_RobPrev,_Sapp);
              return new RobotI<Real,PSecTup,Args...,SecApp> (_RobPrev,_Sapp);
            }

            template<class Real, class SecApp, class PSecTup, class ...Args >
            CTR_INLINE static RobotBase<Real,PSecTup>* appendSec(const RobotI_<Real,PSecTup,Args...>& _RobPrev, const SecApp& _Sapp) noexcept
            {
                typedef RobotI <Real,PSecTup,Args...,SecApp> R ;
                //typedef RobotI_<Real,PSecTup,Args...       > P0;
                //typedef SecApp                               P1;
                //return Erl::alloc_aligned_init_array< R, P0, P1 >(size_t(1),R::align_value,_RobPrev,_Sapp);
                return new R(_RobPrev,_Sapp);
            }
        };

        template<>
        struct append_helper<CTR_MAXSECTION>
        {
            template<class Real,  class SecApp, class PSecTup, class ...Args>
            CTR_INLINE static RobotBase<Real,PSecTup>* appendSec(const RobotI<Real,PSecTup,Args...>& _Ts,const SecApp& _Sapp) noexcept
            { return nullptr; }

            template<class Real,  class SecApp, class PSecTup, class ...Args>
            CTR_INLINE static RobotBase<Real,PSecTup>* appendSec(const RobotI_<Real,PSecTup,Args...>& _Ts,const SecApp& _Sapp) noexcept
            { return nullptr; }
        };
    }
template<class Real, class PSecTup>
class alignas(CTR_ALIGNMENT) RobotBase
{
public:
    typedef PSecTup                  PSecTupRBase;
    typedef Real                     RealRBase   ;
    typedef RobotBase<Real, PSecTup> this_type   ;

    typedef Erl::Transform<Real>     Transform;
    typedef Erl::Rotmat   <Real>     Rotmat   ;
    typedef Erl::Vector6  <Real>     Vector6  ;
    typedef Erl::Vector3  <Real>     Vector3  ;
    typedef Eigen::Matrix<Real,Eigen::Dynamic,1,Eigen::ColMajor> VectorJ    ;
    typedef Eigen::Matrix<Real,6,Eigen::Dynamic,Eigen::ColMajor> MatJacobian;

    template<long iE>
    using ST = typename std::decay<typename std::tuple_element<iE,PSecTup>::type>::type;

    static constexpr size_t NPSecTupType = std::tuple_size<PSecTup >::value;
    virtual CTR_INLINE ~RobotBase<Real,PSecTup>() noexcept {/*static_assert(NPSecTupType == NARGTYPE , "wrong number of PSecTup and NARGTYPE");*/ }

    virtual this_type* appendSec(const ST<0> & _section) const noexcept   = 0;
    virtual this_type* appendSec(const ST<1> & _section) const noexcept   = 0;
    virtual size_t                   requiredMemory()                  const = 0;
    virtual size_t                   requiredMemory(const size_t & _max_samp,size_t& _alignment) const = 0;
    virtual size_t                   requiredMemory(size_t& _alignment)const = 0;
    virtual this_type* move_into(const size_t & _max_samp, const size_t & _memory_size, void * _memory) const = 0;
    virtual this_type* move_into(const size_t & _memory_size, void * _memory) const = 0;

    virtual size_t  getNSections       () const noexcept = 0;
    virtual size_t  getNTubes          () const noexcept = 0;
    virtual size_t  getDoF             () const noexcept = 0;
    virtual size_t  getNJoints         () const noexcept = 0;
    virtual size_t  getMaxSamples      () const noexcept = 0;
    virtual size_t  getNSamples        () const noexcept = 0;
    virtual bool    getStable          () const noexcept = 0;
    virtual Real    getMaxLength       () const noexcept = 0;
    virtual size_t  getMemoryAlignment () const noexcept = 0;
    virtual VectorJ getJointValues     () const noexcept = 0;
    virtual void    getJointValues     (VectorJ& _JV) const noexcept = 0;

    virtual Real       getRadiusSample     (unsigned _Sample) const noexcept = 0;
    virtual Transform  getTransformSample  (unsigned _Sample) const noexcept = 0;
    virtual Vector3    getTranslationSample(unsigned _Sample) const noexcept = 0;

    virtual Real       getTipRadius     () const noexcept = 0;
    virtual Transform  getTipTransform  () const noexcept = 0;
    virtual Vector3    getTipTranslation() const noexcept = 0;

    virtual Transform calcKinematic(const VectorJ& _JVal,
                                    const Real& _ArcLengthStep,Real& _ArcLengthStepReturn) noexcept= 0;

    virtual Transform calcKinematic(const VectorJ& _JVal, const size_t& _NSamples, Real& _ArcLengthStepReturn) noexcept = 0;

            Transform calcKinematic(const VectorJ& _JVal, const Real& _ArcLengthStep) noexcept
            { Real c_ArcR; return calcKinematic(_JVal,_ArcLengthStep,c_ArcR);}

            Transform calcKinematic(const VectorJ& _JVal, const size_t& _NSamples) noexcept
            { Real c_ArcR; return calcKinematic(_JVal,_NSamples,c_ArcR);}

    virtual bool   calcStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep) = 0;
    virtual Real   calcDegreeStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep, const Real& _MinDegreeStability) = 0;

    virtual void   calcJacobian(const VectorJ& _JVal, const Real& _ArcLengthStep, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip) noexcept = 0;
    virtual void   calcJacobian(const VectorJ& _JVal, const size_t& _NSamples,    const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample) noexcept = 0;
    virtual void   calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const size_t& _NSamples, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip) noexcept = 0;
    virtual void   calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const Transform& _TiSample, const size_t& _NSamples,    const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample) noexcept = 0;


    virtual void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn) const = 0;
    virtual void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn, const TransRotScale<Real>& _TS) const = 0;

    virtual void getSection(unsigned _iSection, const Section<Real,sec::CC> ** _Section) const = 0;
    virtual void getSection(unsigned _iSection, const Section<Real,sec::VC> ** _Section) const = 0;



    #ifdef WITH_CTROPENCL
    virtual size_t copy2opencl(void *_robot_cl, int _mode) const noexcept = 0;
    #endif

    virtual std::ostream&  print(std::ostream& os) const noexcept = 0;
    friend  std::ostream&  operator<< (std::ostream& os, const RobotBase &obj)
    { return obj.print(os); }

};



template<class Real, class PSecTup, class ...Args>
class
alignas(CTR_ALIGNMENT)
RobotI_ : public RobotBase<Real,PSecTup>
{
  public:
    typedef Real                                  value_type;
    typedef std::tuple<Args... >                  tuple_type;
    typedef RobotI_<Real, PSecTup, Args ... >     this_type ;
    typedef RobotBase<Real,PSecTup>               base_type ;
    static constexpr size_t tuple_size = std::tuple_size<tuple_type>::value;
    static constexpr size_t n_sections = tuple_size;
    static constexpr size_t n_tubes    = sec::count_tubes<tuple_type>::value;
    static constexpr size_t n_jointval = n_tubes+n_sections+1;

    typedef typename base_type::Transform         Transform  ;
    typedef typename base_type::Rotmat            Rotmat     ;
    typedef typename base_type::Vector6           Vector6    ;
    typedef typename base_type::Vector3           Vector3    ;
    typedef typename base_type::VectorJ           VectorJ    ;
    typedef typename base_type::MatJacobian       MatJacobian;

    typedef Eigen::Matrix<Real,n_jointval,1,Eigen::ColMajor> VectorJs;

    typedef Matrix3N1<Real,n_tubes   ,CTR_ALIGNMENT> Mat3T1;
    typedef Matrix1N1<Real,n_tubes   ,CTR_ALIGNMENT> Mat1T1;
    typedef Matrix1N1<Real,n_sections,CTR_ALIGNMENT> Mat1S1;
    typedef Matrix1X3<Real,           CTR_ALIGNMENT> Mat1X3;
    typedef Matrix31X<Real,           CTR_ALIGNMENT> Mat31X;
    typedef Matrix1X1<Real,           CTR_ALIGNMENT> Mat1X1;
    typedef Matrix1X1<Transform,      CTR_ALIGNMENT> MatTX1;

    template<long i>
    using ST = typename RobotBase<Real,PSecTup>::template ST<i>;

    template<class EigenDerived>
    using MatT = Eigen::MatrixBase<EigenDerived>;
private:
    friend struct details::append_helper<tuple_size>;
    friend class  RobotI_<Real,PSecTup,Args..., ST<0> >;
    friend class  RobotI_<Real,PSecTup,Args..., ST<1> >;

public:
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, const Args&... _joints) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(std::make_tuple(_joints...))
    {}
    CTR_INLINE explicit RobotI_(const Args&... _joints) noexcept
        : m_MaxSamples(0)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(std::make_tuple(_joints...))
    {}
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, Args&&... _joints) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(std::forward_as_tuple(_joints...))
    {}
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, const std::tuple<Args... >& _robot) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(_robot)
    {}
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, std::tuple<Args... >&& _robot) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(std::move(_robot))
    {}
    CTR_INLINE explicit RobotI_(const std::tuple<Args... >& _robot) noexcept
        : m_MaxSamples(0)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(_robot)
    {}
    CTR_INLINE explicit RobotI_(std::tuple<Args... >&& _robot) noexcept
        : m_MaxSamples(0)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(Transform::Identity())
        , m_Sections(std::move(_robot))
    {}
    // *******************************
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, const Transform& _InitTrafo, const Args&... _joints) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(std::make_tuple(_joints...))
    {}
    CTR_INLINE explicit RobotI_(const Transform& _InitTrafo, const Args&... _joints) noexcept
        : m_MaxSamples(0)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(std::make_tuple(_joints...))
    {}
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, const Transform& _InitTrafo, Args&&... _joints) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(std::forward_as_tuple(_joints...))
    {}
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, const Transform& _InitTrafo, const std::tuple<Args... >& _robot) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(_robot)

    {}
    CTR_INLINE explicit RobotI_(const size_t& _max_samp, const Transform& _InitTrafo, std::tuple<Args... >&& _robot) noexcept
        : m_MaxSamples(_max_samp)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(std::move(_robot))
    {}
    CTR_INLINE explicit RobotI_(const Transform& _InitTrafo, const std::tuple<Args... >& _robot) noexcept
        : m_MaxSamples(0)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(_robot)
    {}
    CTR_INLINE explicit RobotI_(const Transform& _InitTrafo, std::tuple<Args... >&& _robot) noexcept
        : m_MaxSamples(0)
        , m_Samples   (0)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_InitTrafo(_InitTrafo)
        , m_Sections(std::move(_robot))
    {}
    // *******************************
    template<class _R, class _ST>
    CTR_INLINE explicit RobotI_(_R&& _robot, _ST&& _section) noexcept
        : m_MaxSamples(_robot.m_MaxSamples)
        , m_Samples   (_robot.m_Samples)
        , m_MemorySize(requiredMemory(m_MaxSamples))
        , m_RigidTheta     (_robot.m_RigidTheta     )
        , m_ConfigStable   (_robot.m_ConfigStable   )
        , m_InitTrafo      (_robot.m_InitTrafo      )
        , m_Sections(std::tuple_cat(_robot.m_Sections,std::forward_as_tuple((_section))))
    {}
    CTR_INLINE RobotI_(const RobotI_& _robot) noexcept
        : m_MaxSamples     (_robot.m_MaxSamples     )
        , m_Samples        (_robot.m_Samples        )
        , m_MemorySize     (_robot.m_MemorySize     )
        , m_RigidTheta     (_robot.m_RigidTheta     )
        , m_ConfigStable   (_robot.m_ConfigStable   )
        , m_InitTrafo      (_robot.m_InitTrafo      )
        , m_SampleCurvature(_robot.m_SampleCurvature)
        , m_SampleTubeAlpha(_robot.m_SampleTubeAlpha)
        , m_Sections       (_robot.m_Sections       )
    {}
    CTR_INLINE RobotI_(RobotI_&& _robot) noexcept
        : m_MaxSamples     (std::move(_robot.m_MaxSamples     ))
        , m_Samples        (std::move(_robot.m_Samples        ))
        , m_MemorySize     (std::move(_robot.m_MemorySize     ))
        , m_RigidTheta     (std::move(_robot.m_RigidTheta     ))
        , m_ConfigStable   (std::move(_robot.m_ConfigStable   ))
        , m_InitTrafo      (std::move(_robot.m_InitTrafo      ))
        , m_SampleCurvature(std::move(_robot.m_SampleCurvature))
        , m_SampleTubeAlpha(std::move(_robot.m_SampleTubeAlpha))
        , m_Sections       (std::move(_robot.m_Sections       ))
    {}
    CTR_INLINE RobotI_(const RobotI_& _robot, void * _memory) noexcept
        : m_MaxSamples     (_robot.m_MaxSamples     )
        , m_Samples        (_robot.m_Samples        )
        , m_MemorySize     (_robot.m_MemorySize     )
        , m_RigidTheta     (_robot.m_RigidTheta     )
        , m_ConfigStable   (_robot.m_ConfigStable   )
        , m_InitTrafo      (_robot.m_InitTrafo      )
        , m_SampleCurvature(_robot.m_SampleCurvature)
        , m_SampleTubeAlpha(_robot.m_SampleTubeAlpha)
        , m_Sections       (_robot.m_Sections       )
    {
        std::ptrdiff_t c_MemoryPtr = std::ptrdiff_t(_memory)+ Erl::align_padded_size(sizeof(this_type),CTR_ALIGNMENT);
        m_CentreLineTransform   .setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples  ,1);
        c_MemoryPtr += MatTX1::requiredPtrMemSizePadded(m_MaxSamples,1);
        m_CentreLineRadius      .setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples  ,1);
        c_MemoryPtr += Mat1X1::requiredPtrMemSizePadded(m_MaxSamples,1);
        m_AlphaInnermostTube    .setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples  ,3);
        c_MemoryPtr += Mat1X3::requiredPtrMemSizePadded(m_MaxSamples,3);
        m_CurvatureInnermostTube.setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples*3,1);
    }
    CTR_INLINE RobotI_(const RobotI_& _robot, const size_t& _max_samp, void * _memory) noexcept
        : m_MaxSamples     (_max_samp               )
        , m_Samples        (0                       )
        , m_MemorySize     ( requiredMemory()       )
        , m_RigidTheta     (_robot.m_RigidTheta     )
        , m_ConfigStable   (_robot.m_ConfigStable   )
        , m_InitTrafo      (_robot.m_InitTrafo      )
        , m_SampleCurvature(_robot.m_SampleCurvature)
        , m_SampleTubeAlpha(_robot.m_SampleTubeAlpha)
        , m_Sections       (_robot.m_Sections       )
    {
        std::ptrdiff_t c_MemoryPtr = std::ptrdiff_t(_memory)+ Erl::align_padded_size(sizeof(this_type),CTR_ALIGNMENT);
        m_CentreLineTransform   .setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples,1);
        c_MemoryPtr += MatTX1::requiredPtrMemSizePadded(int(m_MaxSamples),1);
        m_CentreLineRadius      .setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples,1);
        c_MemoryPtr += Mat1X1::requiredPtrMemSizePadded(int(m_MaxSamples),1);
        m_AlphaInnermostTube    .setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples,3);
        c_MemoryPtr += Mat1X3::requiredPtrMemSizePadded(int(m_MaxSamples),3);
        m_CurvatureInnermostTube.setDataPtr(reinterpret_cast<void*>(c_MemoryPtr),m_MaxSamples*3,1);
    }
    static void* operator new(std::size_t _n)
    { return Erl::aligned_malloc(sizeof(RobotI_)*_n,CTR_ALIGNMENT); }
    static void* operator new[](std::size_t _n)
    { return Erl::aligned_malloc(sizeof(RobotI_)*_n,CTR_ALIGNMENT);}
    static void* operator new(std::size_t, void* _ptr)
    { return _ptr; }
    static void* operator new[](std::size_t, void* _ptr)
    { return _ptr;}
    static void operator delete  ( void* _ptr)
    { return Erl::aligned_free(_ptr);}
    static void operator delete[]( void* _ptr)
    { return Erl::aligned_free(_ptr); }
    CTR_INLINE virtual base_type* appendSec(const ST<0> & _section) const noexcept
    {  typedef typename std::decay<decltype(_section)>::type _T;
       return details::append_helper<tuple_size>::template appendSec<Real,_T,PSecTup,Args... >(*this,_section);
    }
    CTR_INLINE virtual base_type* appendSec(const ST<1> & _section) const noexcept
    {  typedef typename std::decay<decltype(_section)>::type _T;
       return details::append_helper<tuple_size>::template appendSec<Real,_T,PSecTup,Args... >(*this,_section);
    }

    CTR_INLINE const std::tuple<Args... >& getRobot() const noexcept
    { return m_Sections; }

    CTR_INLINE virtual std::ostream&  print(std::ostream& os) const noexcept
    {
        os<<"ConfigStable   : " <<m_ConfigStable   <<std::endl;
        os<<"MaxSamples     : " <<m_MaxSamples     <<std::endl;
        os<<"Samples        : " <<m_Samples        <<std::endl;
        os<<"MemorySize     : " <<m_MemorySize     <<std::endl;
        os<<"RigidTheta     : " <<m_RigidTheta     <<std::endl;
        os<<"InitTrafo      :\n"<<m_InitTrafo      <<std::flush;
        os<<"SampleCurvature: " <<m_SampleCurvature.transpose()<<std::endl;
        os<<"SampleTubeAlpha: " <<m_SampleTubeAlpha.transpose()<<std::endl;
        Erl::static_for<0,tuple_size-1>()([&](auto&& tup_el, auto&& i){ os<<"Section: "<<i<<"\n"<< tup_el <<std::endl; },m_Sections);
        os<<"Section: "<<tuple_size-1<<"\n"<<std::get<tuple_size-1>(m_Sections);
        return os;
    }
    CTR_INLINE virtual size_t getNSections      () const noexcept {return this_type::NSections(); }
    CTR_INLINE virtual size_t getNTubes         () const noexcept {return this_type::NTubes   (); }
    CTR_INLINE virtual size_t getMaxSamples     () const noexcept {return m_MaxSamples          ; }
    CTR_INLINE virtual size_t getNSamples       () const noexcept {return m_Samples             ; }
    CTR_INLINE virtual size_t getDoF            () const noexcept {return this_type::NSections()+this_type::NTubes   ();}
    CTR_INLINE virtual size_t getNJoints        () const noexcept {return 1+getDoF();}
    CTR_INLINE virtual size_t getMemoryAlignment() const noexcept {return alignof(this_type);};


    CTR_INLINE virtual Real       getRadiusSample      (unsigned _Sample) const noexcept { return m_CentreLineRadius   [_Sample]; }
    CTR_INLINE virtual Transform  getTransformSample   (unsigned _Sample) const noexcept { return m_CentreLineTransform[_Sample]; }
    CTR_INLINE const Transform&   getTransformSampleRef(unsigned _Sample) const noexcept { return m_CentreLineTransform[_Sample]; }
    CTR_INLINE virtual Vector3    getTranslationSample (unsigned _Sample) const noexcept { return m_CentreLineTransform[_Sample].getTranslation(); }

    CTR_INLINE virtual Real       getTipRadius     () const noexcept { return m_CentreLineRadius   [m_Samples-1];}
    CTR_INLINE virtual Transform  getTipTransform  () const noexcept { return m_CentreLineTransform[m_Samples-1]; }
    CTR_INLINE virtual Vector3    getTipTranslation() const noexcept { return m_CentreLineTransform[m_Samples-1].getTranslation(); }

    CTR_INLINE bool               getStable()    const noexcept { return m_ConfigStable; }

    CTR_INLINE virtual Real       getMaxLength () const noexcept
    {
        Real c_Length(0);
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        { c_Length += iSec.getLength();}, m_Sections);
        return c_Length;
    }

    CTR_INLINE static constexpr size_t NSections() { return n_sections; }
    CTR_INLINE static constexpr size_t NTubes   () { return n_tubes   ; }
    CTR_INLINE static           size_t requiredMemory(const size_t& _max_samp)
    {
        if(_max_samp == 0)
        { return Erl::align_padded_size(sizeof(this_type) , CTR_ALIGNMENT); }

        return ( Erl::align_padded_size(sizeof(this_type) , CTR_ALIGNMENT)
                 + MatTX1::requiredPtrMemSizePadded(int(_max_samp  ),1)
                 + Mat1X1::requiredPtrMemSizePadded(int(_max_samp  ),1)
                 + Mat1X3::requiredPtrMemSizePadded(int(_max_samp  ),3)
                 + Mat1X1::requiredPtrMemSizePadded(int(3*_max_samp),1));
    }
    CTR_INLINE virtual size_t requiredMemory() const
    { return requiredMemory(m_MaxSamples); }
    CTR_INLINE virtual size_t requiredMemory(size_t& _alignment)const
    { _alignment = CTR_ALIGNMENT; return requiredMemory(); }
    CTR_INLINE virtual size_t requiredMemory(const size_t & _max_samp,size_t& _alignment) const
    { _alignment = CTR_ALIGNMENT; return requiredMemory(_max_samp);}
    CTR_INLINE virtual this_type* move_into(const size_t & _memory_size, void * _memory) const
    { return move_into(m_MaxSamples,_memory_size,_memory);}
    CTR_INLINE virtual this_type* move_into(const size_t & _max_samp, const size_t & _memory_size, void * _memory) const
    {   if(requiredMemory(_max_samp) != _memory_size) { return nullptr; }
        return new (_memory) this_type(*this,_max_samp,_memory);
    }

    CTR_INLINE virtual void getSection(unsigned _iSection, const Section<Real,sec::VC> ** _Section) const
    {
        *_Section = nullptr;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr bool same_type = std::is_same<::CTR::Section<Real,::CTR::sec::VC>,typename std::decay<decltype(iSec)>::type>::value;
            Erl::static_if<same_type>()([&]{if(_iSection == isec) { *_Section = reinterpret_cast<const ::CTR::Section<Real,::CTR::sec::VC>*>(&iSec);}});
        },m_Sections);
    }
    CTR_INLINE virtual void getSection(unsigned _iSection, const Section<Real,sec::CC> ** _Section) const
    {
        *_Section = nullptr;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr bool same_type = std::is_same<::CTR::Section<Real,::CTR::sec::CC>,typename std::decay<decltype(iSec)>::type>::value;
            Erl::static_if<same_type>()([&]{if(_iSection == isec) { *_Section = reinterpret_cast<const ::CTR::Section<Real,::CTR::sec::CC>*>(&iSec);}});
        },m_Sections);
    }
    CTR_INLINE virtual VectorJ getJointValues     () const noexcept
    {
        constexpr auto n_jv = n_jointval;
        VectorJ c_JV(n_jv,1);
        this->getJointValues(c_JV);
        return c_JV;
    }
    CTR_INLINE virtual void    getJointValues     (VectorJ& _JV) const noexcept
    {
        _JV[0] = m_RigidTheta;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec    = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto isec_jv = ::CTR::JointValue_PHI<tuple_type,isec>::ijv;
            _JV[isec_jv] = iSec.getPhi();

            Erl::static_for<0,std::decay<decltype(iSec)>::type::ntube>()([&](auto&& j_ic)
            {
                constexpr auto itub_p  = ::CTR::sec::tube_iter<tuple_type,isec>::value;
                constexpr auto itub_s  = std::decay<decltype(j_ic)>::type::value ;
                constexpr auto itub_jv = ::CTR::JointValue_ALPHA<tuple_type,itub_p+itub_s>::ijv;

                _JV[itub_jv] = iSec.template getTubeAlphaTip<itub_s>();
            });
        }, m_Sections);
    }
    #ifdef WITH_CTROPENCL
    CTR_INLINE virtual size_t copy2opencl(void *_robot_cl,int _mode) const noexcept
    {
        if(_mode == Erl::as_int(cl_ctr_mem_mode::private_mem))
        {
            if(_robot_cl == nullptr)
            { return sizeof(RobotCPP<Real,n_sections,n_tubes>);}
            RobotCPP<Real,n_sections,n_tubes>* c_RobotCL
                    = reinterpret_cast<RobotCPP<Real,n_sections,n_tubes>*>(_robot_cl);

            memset(c_RobotCL,0,sizeof(RobotCPP<Real,n_sections,n_tubes>));

            c_RobotCL->m_MaxSamples   = m_MaxSamples;
            c_RobotCL->m_Samples      = m_Samples   ;
            c_RobotCL->m_RigidTheta   = m_RigidTheta;
            c_RobotCL->m_ConfigStable = m_ConfigStable;
            m_InitTrafo.getData(c_RobotCL->m_InitTrafo.data,Erl::ColMajor);
            memcpy(c_RobotCL->m_SampleCurvature.data,m_SampleCurvature.data(),sizeof(Real)*n_tubes*3);
            memcpy(c_RobotCL->m_SampleTubeAlpha.data,m_SampleTubeAlpha.data(),sizeof(Real)*n_tubes*1);

            Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
            {
                constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
                SectionCPP<Real>& clSec = c_RobotCL->m_Sections[isec];
                clSec.m_NTube                 = iSec.ntube;
                clSec.m_Length                = iSec.getLength();
                clSec.m_CurvedArcStart        = iSec.getCurvedArcStart();
                clSec.m_CurvedArcEnd          = iSec.getCurvedArcEnd();
                clSec.m_Phi                   = iSec.getPhi();
                Erl::static_for<0,std::decay<decltype(iSec)>::type::ntube>()([&](auto&& j_ic)
                {
                    constexpr auto itub = std::decay<decltype(j_ic)>::type::value ;
                    clSec.m_Stiffness     [itub] = iSec.template getTubeStiffness<itub>();
                    clSec.m_Curvature     [itub] = iSec.template getTubeCurvature<itub>();
                    clSec.m_Radius        [itub] = iSec.template getRadius       <itub>();
                    clSec.m_PoissonRatio  [itub] = iSec.template getPossionRatio <itub>();
                    clSec.m_AlphaTip      [itub] = iSec.template getTubeAlphaTip <itub>();
                });
            }, m_Sections);
            return sizeof(RobotCPP<Real,n_sections,n_tubes>);
        }
        else if(_mode == Erl::as_int(cl_ctr_mem_mode::local_constant_mem))
        {
            if(_robot_cl == nullptr)
            { return sizeof(RobotCCPP<Real,n_sections>);}
            RobotCCPP<Real,n_sections>* c_RobotCL
                    = reinterpret_cast<RobotCCPP<Real,n_sections>*>(_robot_cl);

            memset(c_RobotCL,0,sizeof(RobotCCPP<Real,n_sections>));

            c_RobotCL->m_MaxSamples   = m_MaxSamples;
            m_InitTrafo.getData(c_RobotCL->m_InitTrafo.data,Erl::ColMajor);

            Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
            {
                constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
                SectionCCPP<Real>& clSec = c_RobotCL->m_Sections[isec];
                clSec.m_NTube                 = iSec.ntube;
                clSec.m_Length                = iSec.getLength();
                Erl::static_for<0,std::decay<decltype(iSec)>::type::ntube>()([&](auto&& j_ic)
                {
                    constexpr auto itub = std::decay<decltype(j_ic)>::type::value ;
                    clSec.m_Stiffness     [itub] = iSec.template getTubeStiffness<itub>();
                    clSec.m_Curvature     [itub] = iSec.template getTubeCurvature<itub>();
                    clSec.m_Radius        [itub] = iSec.template getRadius       <itub>();
                    clSec.m_PoissonRatio  [itub] = iSec.template getPossionRatio <itub>();
                });
            }, m_Sections);
            return sizeof(RobotCCPP<Real,n_sections>);
        }
        return 0;
    }
    #endif
    CTR_INLINE virtual Transform calcKinematic(const VectorJ & _JVal,
                                           const Real& _ArcLengthStep,Real& _ArcLengthStepReturn) noexcept
    {
        _ArcLengthStepReturn=_ArcLengthStep;
        Real   c_CurrArcPosition;
        Real   c_RobotArcEnd;
        size_t c_CurrSample(0);

        calcAnglesCurvature(_JVal,_ArcLengthStepReturn,c_CurrArcPosition,c_RobotArcEnd,c_CurrSample);
        Transform c_InterMediateFrame = Transform::Identity();
        size_t c_SampleMax = std::max<size_t>( size_t(1), size_t(c_RobotArcEnd/_ArcLengthStepReturn));

        m_AlphaInnermostTube.topRows(c_SampleMax).col(1).array() = (m_AlphaInnermostTube.topRows(c_SampleMax).template leftCols<1>().array()+m_RigidTheta).sin();
        m_AlphaInnermostTube.topRows(c_SampleMax).col(2).array() = (m_AlphaInnermostTube.topRows(c_SampleMax).template leftCols<1>().array()+m_RigidTheta).cos();

        for(;
            (    (c_CurrArcPosition <= c_RobotArcEnd)
              && (c_CurrSample<m_MaxSamples) );
            c_CurrArcPosition += _ArcLengthStepReturn, ++c_CurrSample)

        {
            assert( c_CurrSample<m_MaxSamples );
            c_InterMediateFrame *= curvature2Transform(Vector3((m_CurvatureInnermostTube.data()+3*c_CurrSample)),_ArcLengthStepReturn);

            m_CentreLineTransform[c_CurrSample] = m_InitTrafo* curvature2Transform(c_InterMediateFrame,
                                                                                  m_AlphaInnermostTube(c_CurrSample,1),
                                                                                  m_AlphaInnermostTube(c_CurrSample,2));
        }
        calcTubeRadius(_ArcLengthStepReturn,c_CurrSample);
        m_Samples = c_CurrSample;
        return m_CentreLineTransform[c_CurrSample-1];
    }
    CTR_INLINE virtual Transform calcKinematic(const VectorJ& _JVal, const size_t& _NSamples, Real& _ArcLengthStepReturn) noexcept
    {
        Real   c_CurrArcPosition;
        Real   c_RobotArcEnd;
        size_t c_CurrSample(0);

        calcAnglesCurvature(_JVal,_NSamples,_ArcLengthStepReturn,c_CurrArcPosition,c_RobotArcEnd,c_CurrSample);
        Transform c_InterMediateFrame = Transform::Identity();

        m_AlphaInnermostTube.topRows(_NSamples).col(1).array() = (m_AlphaInnermostTube.topRows(_NSamples).template leftCols<1>().array()+m_RigidTheta).sin();
        m_AlphaInnermostTube.topRows(_NSamples).col(2).array() = (m_AlphaInnermostTube.topRows(_NSamples).template leftCols<1>().array()+m_RigidTheta).cos();

        for(;
            ( //   (c_CurrArcPosition <= c_RobotArcEnd)
                 (c_CurrSample<_NSamples   )
              && (c_CurrSample<m_MaxSamples) );
              //c_CurrArcPosition += _ArcLengthStepReturn, ++c_CurrSample)
              ++c_CurrSample)

        {
            assert( c_CurrSample<m_MaxSamples );
            c_InterMediateFrame *= curvature2Transform(Vector3((m_CurvatureInnermostTube.data()+3*c_CurrSample)),_ArcLengthStepReturn);

            m_CentreLineTransform[c_CurrSample] = m_InitTrafo* curvature2Transform(c_InterMediateFrame,
                                                                                  m_AlphaInnermostTube(c_CurrSample,1),
                                                                                  m_AlphaInnermostTube(c_CurrSample,2));
        }
        //std::stringstream c_ss;
        //c_ss<<int(c_CurrSample)<<" | "<<int(_NSamples)<<std::endl;
        //std::cout<<c_ss.str()<<std::endl;
        calcTubeRadius(_ArcLengthStepReturn,c_CurrSample);
        m_Samples = c_CurrSample;
        return m_CentreLineTransform[c_CurrSample-1];
    }
    CTR_INLINE virtual bool calcStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep)
    {
        m_ConfigStable = true;
        Eigen::Matrix<Real,1+n_sections+n_tubes,1,Eigen::ColMajor> c_JVal = _JVal;

        //RigidTheta = Zero;
        c_JVal(0)      = Real(0);

        Real c_ArcLengthStep(_ArcLengthStep);

        Real c_CurrArcPosition;
        Real c_RobotArcEnd;
        size_t  c_CurrSample;

        Real c_AlphaBasesPrev[n_tubes];

        // LOOP OVER TUBES
        Erl::static_for<1,n_tubes>()([&](auto&& i_ic)
        {
            if(!m_ConfigStable){return m_ConfigStable;}
            constexpr auto idxTube = std::decay<decltype(i_ic)>::type::value ;

            Real c_AlphaTip(c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv));
            Real c_ScanAngle;
            Real c_ScanAngleEnd;

            if(c_AlphaTip <= M_PI)
            {
                c_ScanAngle     = Real(0);
                c_ScanAngleEnd  = c_AlphaTip;
            }
            else
            {
                c_ScanAngle    = c_AlphaTip;
                c_ScanAngleEnd = Real(2.)*Real(M_PI);
            }


            // FIRST ITERATION
            c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv) = c_ScanAngle;
            c_ArcLengthStep = _ArcLengthStep;
            this->calcAnglesCurvature(c_JVal,c_ArcLengthStep,c_CurrArcPosition,c_RobotArcEnd,c_CurrSample);
            memcpy(c_AlphaBasesPrev,m_SampleTubeAlpha.data(),sizeof(Real)*n_tubes);
            c_ScanAngle += _AngleDiscretizationStep;

            // FURTHER ITERATIONS
            for(;
                (c_ScanAngle <= c_ScanAngleEnd);
                 c_ScanAngle += _AngleDiscretizationStep)
            {
                c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv) = c_ScanAngle;
                c_ArcLengthStep = _ArcLengthStep;
                this->calcAnglesCurvature(c_JVal,c_ArcLengthStep,c_CurrArcPosition,c_RobotArcEnd,c_CurrSample);
                // AlphaBase in m_SampleTubeAlpha
                if((m_SampleTubeAlpha[idxTube]-c_AlphaBasesPrev[idxTube])<0.){m_ConfigStable = false; return m_ConfigStable;}
                // Remember AlphaBases
                memcpy(c_AlphaBasesPrev,m_SampleTubeAlpha.data(),sizeof(Real)*n_tubes);
            }
            c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv) = _JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv);
            return m_ConfigStable;
        });
        return m_ConfigStable;
    }
    CTR_INLINE virtual Real      calcDegreeStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep, const Real& _MinDegreeStability)
    {
        m_ConfigStable = true;
        Eigen::Matrix<Real,1+n_sections+n_tubes,1,Eigen::ColMajor> c_JVal = _JVal;

        //RigidTheta = Zero;
        c_JVal(0)      = Real(0);

        Real c_ArcLengthStep(_ArcLengthStep);

        Real    c_CurrArcPosition;
        Real    c_RobotArcEnd;
        size_t  c_CurrSample;

        Real c_AlphaBasesPrev[n_tubes];
        bool c_RelBaseRotvsTipRot_checked= (_MinDegreeStability > Erl::Constants<Real>::Zero_Tolerance ? false : true);

        Real c_SmallestRelBaseRotvsTipRot= std::numeric_limits<Real>::max();
        Real c_CurrentRelBaseRotvsTipRot;
        Real c_StabilityDegree;

        // LOOP OVER TUBES
        Erl::static_for<1,n_tubes>()([&](auto&& i_ic)
        {
            if(!m_ConfigStable){return c_StabilityDegree;}
            constexpr auto idxTube = std::decay<decltype(i_ic)>::type::value ;

            Real c_AlphaTip(c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv));
            Real c_ScanAngle;
            Real c_ScanAngleEnd;

            if(c_AlphaTip <= ERL_PI)
            {
                c_ScanAngle     = Real(0);
                c_ScanAngleEnd  = c_AlphaTip;
            }
            else
            {
                c_ScanAngle    = c_AlphaTip;
                c_ScanAngleEnd = Real(2.)*Real(M_PI);
            }


            // FIRST ITERATION
            c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv) = c_ScanAngle;
            c_ArcLengthStep = _ArcLengthStep;
            this->calcAnglesCurvature(c_JVal,c_ArcLengthStep,c_CurrArcPosition,c_RobotArcEnd,c_CurrSample);
            memcpy(c_AlphaBasesPrev,m_SampleTubeAlpha.data(),sizeof(Real)*n_tubes);
            c_ScanAngle += _AngleDiscretizationStep;

            // FURTHER ITERATIONS
            for(;
                (c_ScanAngle <= c_ScanAngleEnd);
                 c_ScanAngle += _AngleDiscretizationStep)
            {
                c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv) = c_ScanAngle;
                c_ArcLengthStep = _ArcLengthStep;
                this->calcAnglesCurvature(c_JVal,c_ArcLengthStep,c_CurrArcPosition,c_RobotArcEnd,c_CurrSample);
                // AlphaBase in m_SampleTubeAlpha
                c_CurrentRelBaseRotvsTipRot = (m_SampleTubeAlpha[idxTube]-c_AlphaBasesPrev[idxTube]);
                if(c_CurrentRelBaseRotvsTipRot<c_SmallestRelBaseRotvsTipRot)
                {
                    c_RelBaseRotvsTipRot_checked = true;
                    c_SmallestRelBaseRotvsTipRot = c_CurrentRelBaseRotvsTipRot;
                    if(c_CurrentRelBaseRotvsTipRot<0)
                    {
                        m_ConfigStable = false;
                        return (c_StabilityDegree = std::atan2(_AngleDiscretizationStep,c_CurrentRelBaseRotvsTipRot));
                    }
                }
                // Remember AlphaBases
                memcpy(c_AlphaBasesPrev,m_SampleTubeAlpha.data(),sizeof(Real)*n_tubes);
            }
            c_JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv) = _JVal(::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv);
            return c_StabilityDegree;
        });

        if(c_RelBaseRotvsTipRot_checked)
        { return (c_StabilityDegree = std::atan2(_AngleDiscretizationStep,c_SmallestRelBaseRotvsTipRot)); }
        return (_MinDegreeStability+Erl::Constants<Real>::Zero_Tolerance);
    }

    CTR_INLINE virtual void   calcJacobian(const VectorJ& _JVal, const Real& _ArcLengthStep, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot,
                                           MatJacobian& _JacobianTip) noexcept
    {
        Real c_ArcLengthStepRet;
        Transform c_Tip = calcKinematic(_JVal,_ArcLengthStep,c_ArcLengthStepRet);
        calcJacobian(_JVal,c_Tip,m_Samples,_FiniteDiff_Trans,_FiniteDiff_Rot,_JacobianTip);
    }

    CTR_INLINE virtual void   calcJacobian(const VectorJ& _JVal, const size_t& _NSamples, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot,
                                           const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample) noexcept
    {
        Real c_ArcLengthStepRet;
        Transform c_Tip      = calcKinematic(_JVal,_NSamples,c_ArcLengthStepRet);
        Transform c_TiSample = getTransformSample(_iSample);
        calcJacobian(_JVal,c_Tip,c_TiSample,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_iSample,_JacobianTip,_JacobianiSample);
    }
    CTR_INLINE virtual void   calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const size_t& _NSamples,
                                           const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip) noexcept
    {
        /**/
        Real c_ArcLengthStepRet;
        //c_T0 = calcKinematic(_JVal,m_Samples,c_ArcLengthStepRet);
        Transform c_T_tip;
        VectorJs  c_JVal = _JVal;
        Real      c_JValBackup ;
        Real      c_FiniteDiff ;

         _JacobianTip.col(2).setZero();
        Erl::static_for<1,n_tubes>()([&](auto&& i_ic)
        {
            constexpr auto idxTube = std::decay<decltype(i_ic)>::type::value;
            constexpr auto idxJV   = ::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv;

            c_JValBackup =  _JVal[idxJV];
            c_JVal[idxJV] = c_JValBackup+_FiniteDiff_Rot;
            c_T_tip = this->calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);
            ::CTR::JacobianColumnFiniteDifference(_TTip,c_T_tip,_FiniteDiff_Rot,_JacobianTip.col(idxJV).data());
            c_JVal[idxJV] = c_JValBackup;
        });

        Erl::static_for<0,n_sections>()([&](auto&& i_ic)
        {
            constexpr auto idxSect = std::decay<decltype(i_ic)>::type::value;
            constexpr auto idxJV   = ::CTR::JointValue_PHI<tuple_type,idxSect>::ijv;

            c_JValBackup =  _JVal[idxJV];
            c_JVal[idxJV] = c_JValBackup+_FiniteDiff_Trans;
            if(c_JVal[idxJV] >= std::get<idxSect>(m_Sections).getLength())
            { c_JVal[idxJV] = c_JValBackup-_FiniteDiff_Trans; c_FiniteDiff = -_FiniteDiff_Trans; }
            else{ c_FiniteDiff = _FiniteDiff_Trans; }
            c_T_tip = this->calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);
            ::CTR::JacobianColumnFiniteDifference(_TTip,c_T_tip,c_FiniteDiff,_JacobianTip.col(idxJV).data());
            c_JVal[idxJV] = c_JValBackup;
        });

        _JacobianTip.template bottomLeftCorner<3,1>() = m_InitTrafo.getColumn2();
        _JacobianTip.template topLeftCorner   <3,1>() = m_InitTrafo.getColumn2().cross(
                                                                  _TTip.getTranslationReference()
                                                           -m_InitTrafo.getTranslationReference() );
         /**/
        /*
        Real c_ArcLengthStepRet;

        Transform c_T_lo,c_T_hi;
        VectorJs  c_JVal = _JVal;
        Real      c_JValBackup ;
        Real      c_FinDiff_lo, c_FinDiff_hi;

         _JacobianTip.col(2).setZero();

        Erl::static_for<1,n_tubes>()([&](auto&& i_ic)
        {
            constexpr auto idxTube = std::decay<decltype(i_ic)>::type::value;
            constexpr auto idxJV   = ::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv;

            c_JValBackup =  _JVal[idxJV];
            c_JVal[idxJV] = c_JValBackup+_FiniteDiff_Rot;
            c_T_hi = calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);
            c_JVal[idxJV] = c_JValBackup-_FiniteDiff_Rot;
            c_T_lo = calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);

            ::CTR::JacobianColumnFiniteDifference(c_T_lo,c_T_hi,Real(2.)*_FiniteDiff_Rot,_JacobianTip.col(idxJV).data());
            c_JVal[idxJV] = c_JValBackup;
        });

        Erl::static_for<0,n_sections>()([&](auto&& i_ic)
        {
            constexpr auto idxSect = std::decay<decltype(i_ic)>::type::value;
            constexpr auto idxJV   = ::CTR::JointValue_PHI<tuple_type,idxSect>::ijv;
            c_JValBackup =  _JVal[idxJV];

            c_JVal[idxJV] = c_JValBackup+_FiniteDiff_Trans;
            if(c_JVal[idxJV] > std::get<idxSect>(m_Sections).getLength())
            { c_JVal[idxJV] = std::get<idxSect>(m_Sections).getLength(); }
            c_FinDiff_hi = c_JVal[idxJV]-c_JValBackup;
            c_T_hi = calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);

            c_JVal[idxJV] = c_JValBackup-_FiniteDiff_Trans;
            if(c_JVal[idxJV] < 0 )
            { c_JVal[idxJV] = 0;  }
            c_FinDiff_lo = c_JVal[idxJV]-c_JValBackup;
            c_T_lo = calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);

            ::CTR::JacobianColumnFiniteDifference(c_T_lo,c_T_hi,c_FinDiff_hi-c_FinDiff_lo,_JacobianTip.col(idxJV).data());

            c_JVal[idxJV] = c_JValBackup;
        });

        _JacobianTip.template bottomLeftCorner<3,1>() = m_InitTrafo.getColumn2();
        _JacobianTip.template topLeftCorner<3,1>() = m_InitTrafo.getColumn2().cross(_TTip.getTranslationReference()-m_InitTrafo.getTranslationReference()) ;
        */

    }
    CTR_INLINE virtual void   calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const Transform& _TiSample, const size_t& _NSamples,
                                           const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample) noexcept
    {
            Real c_ArcLengthStepRet;

            Transform c_T_tip;
            VectorJs  c_JVal = _JVal;
            Real      c_JValBackup ;
            Real      c_FiniteDiff ;

             _JacobianTip    .col(2).setZero();
             _JacobianiSample.col(2).setZero();
            Erl::static_for<1,n_tubes>()([&](auto&& i_ic)
            {
                constexpr auto idxTube = std::decay<decltype(i_ic)>::type::value;
                constexpr auto idxJV   = ::CTR::JointValue_ALPHA<tuple_type,idxTube>::ijv;

                c_JValBackup =  _JVal[idxJV];
                c_JVal[idxJV] = c_JValBackup+_FiniteDiff_Rot;
                c_T_tip = this->calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);
                ::CTR::JacobianColumnFiniteDifference(_TTip    ,c_T_tip,
                                                      _FiniteDiff_Rot,_JacobianTip    .col(idxJV).data());
                ::CTR::JacobianColumnFiniteDifference(_TiSample,this->getTransformSampleRef(_iSample),
                                                      _FiniteDiff_Rot,_JacobianiSample.col(idxJV).data());

                c_JVal[idxJV] = c_JValBackup;
            });

            Erl::static_for<0,n_sections>()([&](auto&& i_ic)
            {
                constexpr auto idxSect = std::decay<decltype(i_ic)>::type::value;
                constexpr auto idxJV   = ::CTR::JointValue_PHI<tuple_type,idxSect>::ijv;

                c_JValBackup =  _JVal[idxJV];
                c_JVal[idxJV] = c_JValBackup+_FiniteDiff_Trans;

                if(c_JVal[idxJV] >= std::get<idxSect>(m_Sections).getLength())
                { c_JVal[idxJV] = c_JValBackup-_FiniteDiff_Trans; c_FiniteDiff = -_FiniteDiff_Trans; }
                else{ c_FiniteDiff = _FiniteDiff_Trans; }

                c_T_tip = this->calcKinematic(c_JVal,_NSamples,c_ArcLengthStepRet);
                ::CTR::JacobianColumnFiniteDifference(_TTip    ,c_T_tip,
                                                      c_FiniteDiff,_JacobianTip       .col(idxJV).data());
                ::CTR::JacobianColumnFiniteDifference(_TiSample,this->getTransformSampleRef(_iSample),
                                                      c_FiniteDiff,_JacobianiSample.col(idxJV).data());

                c_JVal[idxJV] = c_JValBackup;
            });

            _JacobianTip.template bottomLeftCorner<3,1>() = m_InitTrafo.getColumn2();
            _JacobianTip.template topLeftCorner   <3,1>() = m_InitTrafo.getColumn2().cross(
                                                                   _TTip.getTranslationReference()
                                                            -m_InitTrafo.getTranslationReference() ) ;

            _JacobianiSample.template bottomLeftCorner<3,1>() = m_InitTrafo.getColumn2();
            _JacobianiSample.template topLeftCorner   <3,1>() = m_InitTrafo.getColumn2().cross(
                                                                    _TiSample.getTranslationReference()
                                                                 -m_InitTrafo.getTranslationReference() ) ;
    }

    CTR_INLINE void calcAnglesCurvature(const VectorJ & _JVal, Real& c_ArcLengthStep, Real& c_CurrArcPosition,Real& c_RobotArcEnd, size_t & c_CurrSample)
    {
        memset(&m_SampleCurvature[0],0,sizeof(m_SampleCurvature)+sizeof(m_SampleTubeAlpha));
        c_RobotArcEnd = setJointsCentreLineStates(_JVal);

        c_CurrSample = std::max<size_t>( size_t(1), size_t(c_RobotArcEnd/c_ArcLengthStep));
        if(c_CurrSample>m_MaxSamples)
        { c_CurrSample  = m_MaxSamples; c_ArcLengthStep = c_RobotArcEnd/Real(c_CurrSample); }

        c_CurrArcPosition =  c_RobotArcEnd;
        c_CurrSample--;

        calcAnglesCurvature(c_ArcLengthStep,c_CurrArcPosition,c_CurrSample);
    }
    CTR_INLINE void calcAnglesCurvature(const VectorJ & _JVal, const size_t& _NSamples, Real& c_ArcLengthStep, Real& c_CurrArcPosition,Real& c_RobotArcEnd, size_t & c_CurrSample)
    {
        memset(&m_SampleCurvature[0],0,sizeof(m_SampleCurvature)+sizeof(m_SampleTubeAlpha));
        c_RobotArcEnd = setJointsCentreLineStates(_JVal);

        c_CurrSample      = std::min<size_t>(_NSamples,m_MaxSamples);
        c_ArcLengthStep   = c_RobotArcEnd/Real(c_CurrSample);

        c_CurrArcPosition = c_RobotArcEnd;
        c_CurrSample-- ;

        calcAnglesCurvature(c_ArcLengthStep,c_CurrArcPosition,c_CurrSample);
    }
    CTR_INLINE void calcAnglesCurvature(Real& c_ArcLengthStep, Real& c_CurrArcPosition, size_t & c_CurrSample)
    {

        Real c_CurvatureX(0);
        Real c_CurvatureY(0);
        Real c_SampleCummStiffness(0);

        m_CurvatureInnermostTube[3*c_CurrSample+2] = Real(0);

        for(;
            (c_CurrArcPosition >= c_ArcLengthStep) && (c_CurrSample >= 1);
             c_CurrArcPosition -= c_ArcLengthStep, --c_CurrSample)
        {
            calcSampleCurvatureStiffnessAlpha(c_CurrArcPosition,c_CurrSample,c_ArcLengthStep,c_CurvatureX,c_CurvatureY,c_SampleCummStiffness);

            c_CurvatureX         = 0;
            c_CurvatureY         = 0;
            c_SampleCummStiffness= 0;
        }
        calcSampleCurvatureStiffnessAlphaLastIteration(c_CurrArcPosition,c_CurrSample,c_CurvatureX,c_CurvatureY,c_SampleCummStiffness);

    }
    template<class EigenDerived>
    CTR_INLINE Real setJointsCentreLineStates(const MatT<EigenDerived>& _JVal)
    {
        m_RigidTheta =   _JVal[0];
        // Apply JointValues
        Real c_CurrArcEnd = 0;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec  = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itube = CTR::sec::tube_iter<tuple_type, isec>::value;
            iSec.template setJointsCentreLineStates<isec,itube>(_JVal,m_SampleTubeAlpha.data(),c_CurrArcEnd);
        }, m_Sections );

        return c_CurrArcEnd;

    }
    CTR_INLINE  void calcSampleCurvatureStiffnessAlpha(const Real& _sample, const size_t& _sampleIdx, const Real& _ArcLengthStep, Real &_curvatureX, Real &_curvatureY, Real& _stiffness)
    {

        alignas(Erl::details::max_alignment<Mat1T1,CTR_ALIGNMENT>::value) Mat1T1 c_sinAlpha;
        alignas(Erl::details::max_alignment<Mat1T1,CTR_ALIGNMENT>::value) Mat1T1 c_cosAlpha;
        alignas(CTR_ALIGNMENT) Mat1S1 c_curvature;
        alignas(CTR_ALIGNMENT) Mat1S1 c_stiffness;

        c_sinAlpha.array()=m_SampleTubeAlpha.array().sin();
        c_cosAlpha.array()=m_SampleTubeAlpha.array().cos();

        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            iSec.getCurvatureStiffness(_sample,c_curvature[isec],c_stiffness[isec]);
        }, m_Sections);
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itub = CTR::sec::tube_iter<tuple_type,isec>::value;
            iSec.template addCurvatureAndStiffness<itub>(c_sinAlpha,c_cosAlpha,c_stiffness[isec],c_curvature[isec],
                                                         _curvatureX,_curvatureY,_stiffness);
        }, m_Sections);
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itub = CTR::sec::tube_iter<tuple_type,isec>::value;
            iSec.template calcTorsionalCurvatureXY<itub>(c_sinAlpha,c_cosAlpha,
                                  _curvatureX,_curvatureY,_stiffness,m_SampleCurvature.data());
        }, m_Sections);

        m_CurvatureInnermostTube[3*(_sampleIdx)]  =m_SampleCurvature[3*(n_tubes-1)];
        m_CurvatureInnermostTube[3*(_sampleIdx)+1]=m_SampleCurvature[3*(n_tubes-1)+1];

         // CALCULATE ALPHAS **********************
         std::get<0>(m_Sections). template calcAlphaFirstSection<0>(_ArcLengthStep,m_SampleTubeAlpha.data(),m_SampleCurvature.data());
         Erl::static_for<1,n_sections>()([&](auto&& iSec, auto&& i_ic)
         {
             constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
             constexpr auto itub = CTR::sec::tube_iter<tuple_type,isec>::value;
             iSec.template calcAlpha<itub>(_ArcLengthStep,m_SampleTubeAlpha.data(),m_SampleCurvature.data());
         }, m_Sections);
         m_AlphaInnermostTube(_sampleIdx,0)=m_SampleTubeAlpha[n_tubes-1];

        // *********** CALCULATE CURVATURE Z
        Real c_CummulatedTorsionalCurvatureZ(0);

        constexpr auto itub_lsec = CTR::sec::tube_iter<tuple_type,n_sections-1>::value;
        std::get<n_sections-1>(m_Sections).template calcTorsionalCurvatureZLastSection<itub_lsec>(
                _ArcLengthStep,c_stiffness[(n_sections-1)],c_curvature[(n_sections-1)],m_SampleCurvature.data(),c_CummulatedTorsionalCurvatureZ);
        Erl::static_for<long(n_sections)-2,-1,-1>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itub = CTR::sec::tube_iter<tuple_type,isec>::value;
            iSec.template calcTorsionalCurvatureZ<itub>(_ArcLengthStep,c_stiffness[isec],c_curvature[isec],m_SampleCurvature.data(),c_CummulatedTorsionalCurvatureZ);
        }, m_Sections);
        std::get<n_sections-1>(m_Sections).template calcTorsionalCurvatureZLastTube<itub_lsec>(
                    c_stiffness[(n_sections-1)],m_SampleCurvature.data(),c_CummulatedTorsionalCurvatureZ);
        // **************************************************
        // SAVE THE Z-Value
        m_CurvatureInnermostTube[3*(_sampleIdx)-1]=m_SampleCurvature[3*(n_tubes-1)+2];
        // *********** *********** *********** ***********
    }

    CTR_INLINE void calcSampleCurvatureStiffnessAlphaLastIteration(const Real& _sample,const int& _sampleIdx,Real &_curvatureX, Real &_curvatureY, Real& _stiffness)
    {

        alignas(Erl::details::max_alignment<Mat1T1,CTR_ALIGNMENT>::value) Mat1T1 c_sinAlpha;
        alignas(Erl::details::max_alignment<Mat1T1,CTR_ALIGNMENT>::value) Mat1T1 c_cosAlpha;
        alignas(CTR_ALIGNMENT) Mat1S1 c_curvature;
        alignas(CTR_ALIGNMENT) Mat1S1 c_stiffness;

        c_sinAlpha.array()=m_SampleTubeAlpha.array().sin();
        c_cosAlpha.array()=m_SampleTubeAlpha.array().cos();

        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            iSec.getCurvatureStiffness(_sample,c_curvature[isec],c_stiffness[isec]);
        }, m_Sections);
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {           
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itub = CTR::sec::tube_iter<tuple_type,isec>::value;
            iSec.template addCurvatureAndStiffness<itub>(c_sinAlpha,c_cosAlpha,c_stiffness[isec],c_curvature[isec],
                                                         _curvatureX,_curvatureY,_stiffness);

        }, m_Sections);
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itub = CTR::sec::tube_iter<tuple_type,isec>::value;
            iSec.template calcTorsionalCurvatureXY<itub>(c_sinAlpha,c_cosAlpha,
                                  _curvatureX,_curvatureY,_stiffness,m_SampleCurvature.data());
        }, m_Sections);

        m_CurvatureInnermostTube[3*(_sampleIdx)  ] = m_SampleCurvature[3*(n_tubes-1)  ];
        m_CurvatureInnermostTube[3*(_sampleIdx)+1] = m_SampleCurvature[3*(n_tubes-1)+1];
        m_AlphaInnermostTube    (_sampleIdx,   0 ) = m_SampleTubeAlpha[n_tubes-1];
    }


    CTR_INLINE Transform curvature2Transform(Vector3&& _curvature, Real _arcStep)
    {
        Real c_CurvatureNorm = _curvature.norm();

        if(!Erl::equal<Real>(c_CurvatureNorm,Real(0)))
        {
            Vector3 w = _curvature/c_CurvatureNorm;

            Real sn = std::sin(_arcStep*c_CurvatureNorm);
            Real cn = std::cos(_arcStep*c_CurvatureNorm);

            Rotmat W (
            ((w.data()[1]*w.data()[1] + w.data()[2]*w.data()[2])*(cn - 1) + 1 ),(   - sn*w.data()[2] - w.data()[0]*w.data()[1]*(cn - 1))           ,(     sn*w.data()[1] - w.data()[0]*w.data()[2]*(cn - 1))          ,
            (     sn*w.data()[2] - w.data()[0]*w.data()[1]*(cn - 1))           ,( (w.data()[0]*w.data()[0] + w.data()[2]*w.data()[2])*(cn - 1) + 1),(   - sn*w.data()[0] - w.data()[1]*w.data()[2]*(cn - 1))          ,
            (   - sn*w.data()[1] - w.data()[0]*w.data()[2]*(cn - 1))           ,(     sn*w.data()[0] - w.data()[1]*w.data()[2]*(cn - 1))           ,( (w.data()[0]*w.data()[0] + w.data()[1]*w.data()[1])*(cn - 1) + 1));

            w = Vector3(_arcStep*w.data()[0]*w.data()[2] + (-(w.data()[1]*(W.data()[0] - 1)) + (W.data()[3]*w.data()[0]))/c_CurvatureNorm,
                        _arcStep*w.data()[1]*w.data()[2] + ( (w.data()[0]*(W.data()[4] - 1)) - (W.data()[1]*w.data()[1]))/c_CurvatureNorm,
                        _arcStep*w.data()[2]*w.data()[2] + ( (W.data()[5]*w.data()[0])       - (W.data()[2]*w.data()[1]))/c_CurvatureNorm
                                  );

            return Transform(W,w);
        }
        return Transform(Real(1),Real(0),Real(0),Real(0),
                         Real(0),Real(1),Real(0),Real(0),
                         Real(0),Real(0),Real(1),_arcStep*c_CurvatureNorm);


    }


    CTR_INLINE Transform curvature2Transform(const Transform& _iT, Real _angle)
    {
        Real cn = std::cos(_angle);
        Real sn = std::sin(_angle);
        return Transform( _iT.getRotationReference().data()[0]*cn - _iT.getRotationReference().data()[1]*sn, _iT.getRotationReference().data()[3]*cn - _iT.getRotationReference().data()[4]*sn, _iT.getRotationReference().data()[6]*cn - _iT.getRotationReference().data()[7]*sn, _iT.getTranslationReference().data()[0]*cn - _iT.getTranslationReference().data()[1]*sn,
                          _iT.getRotationReference().data()[1]*cn + _iT.getRotationReference().data()[0]*sn, _iT.getRotationReference().data()[4]*cn + _iT.getRotationReference().data()[3]*sn, _iT.getRotationReference().data()[7]*cn + _iT.getRotationReference().data()[6]*sn, _iT.getTranslationReference().data()[1]*cn + _iT.getTranslationReference().data()[0]*sn,
                                                                       _iT.getRotationReference().data()[2],                                              _iT.getRotationReference().data()[5],                                              _iT.getRotationReference().data()[8],                                                  _iT.getTranslationReference().data()[2]);
    }
    CTR_INLINE Transform curvature2Transform(const Transform & _iT, Real _sin_angle, Real _cos_angle)
    {
        return Transform( _iT.getRotationReference().data()[0]*_cos_angle - _iT.getRotationReference().data()[1]*_sin_angle, _iT.getRotationReference().data()[3]*_cos_angle - _iT.getRotationReference().data()[4]*_sin_angle, _iT.getRotationReference().data()[6]*_cos_angle - _iT.getRotationReference().data()[7]*_sin_angle, _iT.getTranslationReference().data()[0]*_cos_angle - _iT.getTranslationReference().data()[1]*_sin_angle,
                          _iT.getRotationReference().data()[1]*_cos_angle + _iT.getRotationReference().data()[0]*_sin_angle, _iT.getRotationReference().data()[4]*_cos_angle + _iT.getRotationReference().data()[3]*_sin_angle, _iT.getRotationReference().data()[7]*_cos_angle + _iT.getRotationReference().data()[6]*_sin_angle, _iT.getTranslationReference().data()[1]*_cos_angle + _iT.getTranslationReference().data()[0]*_sin_angle,
                                                                       _iT.getRotationReference().data()[2],                                              _iT.getRotationReference().data()[5],                                              _iT.getRotationReference().data()[8],                                                  _iT.getTranslationReference().data()[2]);
    }
    CTR_INLINE void calcTubeRadius(Real _arcLength, size_t _SamplePoints)
    {
        size_t c_ConstantRadiusStartSample(0);
        size_t c_ConstantRadiusEndSample;
        Real   c_SectionRadius, c_CurvedArcEnd;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            if(c_ConstantRadiusStartSample>=(_SamplePoints)){ return; }

            c_SectionRadius = iSec.getOuterRadius();
            c_CurvedArcEnd  = iSec.getCurvedArcEnd();
            c_ConstantRadiusEndSample = std::min<size_t>(size_t(c_CurvedArcEnd/_arcLength), _SamplePoints);
            assert(c_ConstantRadiusEndSample>=c_ConstantRadiusStartSample);
            std::fill(m_CentreLineRadius.data()+c_ConstantRadiusStartSample,
                      m_CentreLineRadius.data()+c_ConstantRadiusEndSample,
                      c_SectionRadius);

            c_ConstantRadiusStartSample = c_ConstantRadiusEndSample;
        }, m_Sections);
        ERL_RMV_UNINIT_SINGLE;
        m_CentreLineRadius[_SamplePoints-1] = c_SectionRadius;
    }
    CTR_INLINE void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn) const
    {
        _JValReturn[0]=std::nextafter(Real(2*M_PI),Real(0))*_RND(_RNG); // RigidTheta;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec  = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itube = CTR::sec::tube_iter<tuple_type, isec>::value;
            iSec.template getRandomJointValues<isec,itube>(_RNG,_RND,_JValReturn.data());
        }, m_Sections );
    }
    CTR_INLINE void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn, const TransRotScale<Real>& _TS) const
    {
        _JValReturn[0]=std::nextafter(Real(2*M_PI),Real(0))*_RND(_RNG); // RigidTheta;
        Erl::static_for<0,n_sections>()([&](auto&& iSec, auto&& i_ic)
        {
            constexpr auto isec  = std::decay<decltype(i_ic)>::type::value ;
            constexpr auto itube = CTR::sec::tube_iter<tuple_type, isec>::value;
            iSec.template getRandomJointValues<isec,itube>(_RNG,_RND,_JValReturn.data(),_TS);
        }, m_Sections );
    }

protected:
    alignas(CTR_ALIGNMENT) const size_t m_MaxSamples     ;
    alignas(CTR_ALIGNMENT)       size_t m_Samples        ;
    alignas(CTR_ALIGNMENT) const size_t m_MemorySize     ;
    alignas(CTR_ALIGNMENT)       Real   m_RigidTheta     ;
    alignas(CTR_ALIGNMENT)       bool   m_ConfigStable   ;
    alignas(CTR_ALIGNMENT) Transform    m_InitTrafo      ;
    alignas(Erl::details::max_alignment<Mat3T1,CTR_ALIGNMENT>::value) Mat3T1 m_SampleCurvature;
    alignas(Erl::details::max_alignment<Mat1T1,CTR_ALIGNMENT>::value) Mat1T1 m_SampleTubeAlpha;
    //alignas(CTR_ALIGNMENT) Mat3T1       m_SampleCurvature;
    //alignas(CTR_ALIGNMENT) Mat1T1       m_SampleTubeAlpha;

    alignas(CTR_ALIGNMENT) tuple_type   m_Sections       ;


    alignas(CTR_ALIGNMENT) MatTX1   m_CentreLineTransform    ;  // PointerToArray -> Depending on MaxSamples
    alignas(CTR_ALIGNMENT) Mat1X1   m_CentreLineRadius       ;  // PointerToArray -> Depending on MaxSamples
    alignas(CTR_ALIGNMENT) Mat1X3   m_AlphaInnermostTube     ;
    alignas(CTR_ALIGNMENT) Mat1X1   m_CurvatureInnermostTube ;  // PointerToArray -> Depending on MaxSamples [actualy X3 [fake impression here] ]



};

template<class Real, class PSecTup, class ...Args>
class
alignas(Erl::details::max_alignment<RobotI_<Real,PSecTup, Args ... >,CTR_ALIGNMENT>::value)
RobotI : public RobotI_<Real,PSecTup,Args ...>
{
public:
    typedef RobotI_<Real,PSecTup,Args ...> RobotIT   ;
    typedef RobotIT                        base_type ;
    typedef typename base_type::Transform  Transform ;
    static constexpr size_t align_value = Erl::details::max_alignment<RobotI_<Real,PSecTup, Args ... >,CTR_ALIGNMENT>::value;


    CTR_INLINE explicit RobotI(const size_t& _max_samp, const Args&... _joints) noexcept
        : RobotIT(_max_samp,_joints...)
    {}
    CTR_INLINE explicit RobotI(const Args&... _joints) noexcept
        : RobotIT(_joints ...)
    {}
    CTR_INLINE explicit RobotI(const size_t& _max_samp, Args&&... _joints) noexcept
        : RobotIT(_max_samp, _joints...)
    {}
    CTR_INLINE explicit RobotI(const size_t& _max_samp, const std::tuple<Args... >& _robot) noexcept
        : RobotIT(_max_samp,_robot)
    {}
    CTR_INLINE explicit RobotI(const size_t& _max_samp, std::tuple<Args... >&& _robot) noexcept
        : RobotIT(_max_samp,_robot)
    {}
    CTR_INLINE explicit RobotI(const std::tuple<Args... >& _robot) noexcept
        : RobotIT(_robot)
    {}
    CTR_INLINE explicit RobotI(std::tuple<Args... >&& _robot) noexcept
        : RobotIT(_robot)
    {}
    // *******************************
    CTR_INLINE explicit RobotI(const size_t& _max_samp, const Transform& _InitTrafo, const Args&... _joints) noexcept
        : RobotIT(_max_samp,_InitTrafo,_joints...)
    {}
    CTR_INLINE explicit RobotI(const Transform& _InitTrafo, const Args&... _joints) noexcept
        : RobotIT(_InitTrafo,_joints...)
    {}
    CTR_INLINE explicit RobotI(const size_t& _max_samp, const Transform& _InitTrafo, Args&&... _joints) noexcept
        : RobotIT(_max_samp,_InitTrafo,_joints...)
    {}
    CTR_INLINE explicit RobotI(const size_t& _max_samp, const Transform& _InitTrafo, const std::tuple<Args... >& _robot) noexcept
        : RobotIT(_max_samp,_InitTrafo,_robot)
    {}
    CTR_INLINE explicit RobotI(const size_t& _max_samp, const Transform& _InitTrafo, std::tuple<Args... >&& _robot) noexcept
        : RobotIT(_max_samp,_InitTrafo,_robot)
    {}
    CTR_INLINE explicit RobotI(const Transform& _InitTrafo, const std::tuple<Args... >& _robot) noexcept
        : RobotIT(_InitTrafo,_robot)
    {}
    CTR_INLINE explicit RobotI(const Transform& _InitTrafo, std::tuple<Args... >&& _robot) noexcept
        : RobotIT(_InitTrafo,_robot)
    {}
    // *******************************
    template<class _R, class _ST>
    CTR_INLINE explicit RobotI(_R&& _robot, _ST&& _section) noexcept
        : RobotIT(_robot,_section)
    {}
    CTR_INLINE RobotI(const RobotIT& _robot) noexcept
        : RobotIT(_robot)
    {}
    CTR_INLINE RobotI(RobotIT&& _robot) noexcept
        : RobotIT(_robot)
    {}
    CTR_INLINE RobotI(const RobotIT& _robot, void * _memory) noexcept
        : RobotIT(_robot,_memory)
    {}
    CTR_INLINE RobotI(const RobotIT& _robot, const size_t& _max_samp, void * _memory) noexcept
        : RobotIT(_robot,_max_samp,_memory)
    {}
    // *******************************
    CTR_INLINE RobotI(const RobotI& _robot) noexcept
        : RobotIT(_robot)
    {}
    CTR_INLINE RobotI(RobotI&& _robot) noexcept
        : RobotIT(_robot)
    {}
    CTR_INLINE RobotI(const RobotI& _robot, void * _memory) noexcept
        : RobotIT(_robot,_memory)
    {}
    CTR_INLINE RobotI(const RobotI& _robot, const size_t& _max_samp, void * _memory) noexcept
        : RobotIT(_robot,_max_samp,_memory)
    {}
    CTR_INLINE size_t getMemoryAlignment() const noexcept {return align_value;}
    static void* operator new(std::size_t _n)
    { return Erl::aligned_malloc(sizeof(RobotI)*_n,align_value); }
    static void* operator new[](std::size_t _n)
    { return Erl::aligned_malloc(sizeof(RobotI)*_n,align_value);}
    static void* operator new(std::size_t, void* _ptr)
    { return _ptr; }
    static void* operator new[](std::size_t, void* _ptr)
    { return _ptr;}
    static void operator delete  ( void* _ptr)
    { return Erl::aligned_free(_ptr);}
    static void operator delete[]( void* _ptr)
    { return Erl::aligned_free(_ptr); }
};

template<class RBase, class ...ContArgs>
RBase * createRobot(const ContArgs&... _Cont)
{ return createRobot<RBase,ContArgs...>(Erl::Transform<typename RBase::RealRBase>::Identity(),_Cont...); }

template<class RBase, class ...ContArgs>
RBase * createRobot(const Erl::Transform<typename RBase::RealRBase>& _InitTrafo, const ContArgs&... _Cont)
{
    //typedef RobotBase<Real,PSecTup>      RBase       ;
    typedef typename RBase::RealRBase    RealRBase   ;
    typedef typename RBase::PSecTupRBase PSecTupRBase;

    typedef typename std::tuple< ContArgs... >                                             ContTuple    ;
    typedef typename std::tuple< typename std::decay<ContArgs>::type::const_iterator ... > ContIterTuple;
    typedef typename std::tuple_size<ContTuple>::value_type size_type;
    constexpr size_type NCont = std::tuple_size<ContTuple>::value ;



    size_type     c_NumEl(0);
    ContIterTuple c_ContPos;
    ContTuple     c_ContTup = std::make_tuple(_Cont...);
    Erl::static_for<0,NCont>()([&c_NumEl](auto&& tup_c,auto&& tup_i, auto&& i){ tup_i = tup_c.begin(); c_NumEl+=tup_c.size();},c_ContTup,c_ContPos);

    RBase* c_RobotPtr(nullptr);

    for(size_type iEl(0);iEl<c_NumEl;iEl++)
    {
        Erl::static_for<0,NCont>()([&](auto&& tup_i,auto&& tup_c, int i)
        {
            typedef typename std::decay<decltype(tup_i->first)>::type E0;
            //E0::print_type(std::cout)<<std::endl;
            if((tup_i!=tup_c.end()) && (size_type(tup_i->second) == iEl))
            {
                if(c_RobotPtr == nullptr)
                {
                    std::tuple<E0> c_AppSection = std::make_tuple(tup_i->first);
                    c_RobotPtr = new CTR::RobotI_<RealRBase, PSecTupRBase ,E0>(size_t(0),_InitTrafo,c_AppSection);
                }
                else
                {
                    RBase* c_RobotPtr_tmp = c_RobotPtr->appendSec(tup_i->first);
                    delete c_RobotPtr; c_RobotPtr = c_RobotPtr_tmp;
                }
                tup_i++;
            }
        },
        c_ContPos,c_ContTup);
    }
    return c_RobotPtr;
}

template<class RBase>
RBase * createRobot2Ptr(const SectionCC<typename RBase::RealRBase>* _SecCC, const std::vector<size_t>& _SecCC_I, const SectionVC<typename RBase::RealRBase>* _SecVC, const std::vector<size_t>& _SecVC_I)
{ return createRobot2Ptr<RBase>(Erl::Transform<typename RBase::RealRBase>::Identity(),_SecCC,_SecCC_I,_SecVC,_SecVC_I); }

template<class RBase>
RBase * createRobot2Ptr(const Erl::Transform<typename RBase::RealRBase>& _InitTrafo,
                        const SectionCC<typename RBase::RealRBase>* _SecCC, const std::vector<size_t>& _SecCC_I,
                        const SectionVC<typename RBase::RealRBase>* _SecVC, const std::vector<size_t>& _SecVC_I)
{
    typedef typename RBase::RealRBase    RealRBase   ;
    typedef typename RBase::PSecTupRBase PSecTupRBase;


    size_t    c_NumEl(0);
    c_NumEl = _SecCC_I.size()+_SecVC_I.size();

    RBase* c_RobotPtr(nullptr);
    size_t c_IterCC = 0;
    size_t c_IterVC = 0;

    for(size_t iEl(0);iEl<c_NumEl;iEl++)
    {
        // SECTION CC
        if((c_IterCC!=_SecCC_I.size() ) && (size_t(_SecCC_I[c_IterCC]) == iEl))
        {
            if(c_RobotPtr == nullptr)
            {
                //typedef CTR::RobotI<RealRBase, PSecTupRBase ,SectionCC<RealRBase> > R_t;
                std::tuple<SectionCC<RealRBase>> c_AppSection = std::make_tuple(_SecCC[c_IterCC]);
                //c_RobotPtr = Erl::alloc_aligned_init_array< R_t>(size_t(1),R_t::align_value, size_t(0),_InitTrafo,c_AppSection);
                c_RobotPtr = new CTR::RobotI<RealRBase, PSecTupRBase ,SectionCC<RealRBase> >(size_t(0),_InitTrafo,c_AppSection);
            }
            else
            {
                RBase* c_RobotPtr_tmp = c_RobotPtr->appendSec(_SecCC[c_IterCC]);
                //Erl::free_aligned_init_array(c_RobotPtr,size_t(1));
                delete c_RobotPtr;
                c_RobotPtr = c_RobotPtr_tmp;
            }
            c_IterCC++;
        }
        // SECTION VC
        if((c_IterVC!=_SecVC_I.size()) && (size_t(_SecVC_I[c_IterVC]) == iEl))
        {
            if(c_RobotPtr == nullptr)
            {
                //typedef CTR::RobotI<RealRBase, PSecTupRBase ,SectionVC<RealRBase> > R_t;
                std::tuple<SectionVC<RealRBase>> c_AppSection = std::make_tuple(_SecVC[c_IterVC]);
                //c_RobotPtr = Erl::alloc_aligned_init_array< R_t>(size_t(1),R_t::align_value, size_t(0),_InitTrafo,c_AppSection);
                c_RobotPtr = new CTR::RobotI<RealRBase, PSecTupRBase ,SectionVC<RealRBase> >(size_t(0),_InitTrafo,c_AppSection);
            }
            else
            {
                RBase* c_RobotPtr_tmp = c_RobotPtr->appendSec(_SecVC[c_IterVC]);
                //Erl::free_aligned_init_array(c_RobotPtr,size_t(1));
                delete c_RobotPtr;
                c_RobotPtr = c_RobotPtr_tmp;
            }
            c_IterVC++;
        }

    }
    return c_RobotPtr;
}

}
#endif // CTR_ROBOT_I_H
