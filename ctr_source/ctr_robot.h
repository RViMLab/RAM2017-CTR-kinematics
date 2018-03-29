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

#ifndef CTR_ROBOT_H
#define CTR_ROBOT_H

#include <string>
#include <vector>
#include <deque>
#include <list>
#include <tuple>
#include <type_traits>

#include <Erl/Core.h>

#include "_forwards.h"

#ifdef WITH_CTROPENCL
template <class T>
class CTR_IMEXPORT BufferElementAbstract;

template <class T>
class CTR_IMEXPORT Buffer;
#endif

namespace CTR
{
    template<class Real, class PSecTup>
    class CTR_IMEXPORT RobotBase;

    template <class Real, class NTube>
    class  CTR_IMEXPORT Section;

    template <class Real>
    class  CTR_IMEXPORT TransRotScale;


    #ifdef WITH_CTROPENCL
    template <class Real>
    class CTR_IMEXPORT ManSampleContainer;

    template <class T>
    class CTR_IMEXPORT ManSampleAbstract;

    #endif

    namespace types {
        typedef std::integral_constant<size_t,1> CCic    ;
        typedef std::integral_constant<size_t,2> VCic    ;

        template<class Real>
        using  SecCC = Section<Real,types::CCic> ;
        template<class Real>
        using  SecVC = Section<Real,types::VCic> ;
        template<class Real>
        using PSTuple = std::tuple<SecCC<Real>,SecVC<Real>>;
        template<class Real>
        using RobotA = RobotBase<Real,PSTuple<Real>>;
    }


    template<class Real>
    class CTR_IMEXPORT Robot
    {

        typedef types::SecCC<Real>        SecCC ;
        typedef types::SecVC<Real>        SecVC ;
        typedef types::PSTuple<Real>      PSTuple;


        public:
            typedef Real                     RReal;
            typedef Erl::Transform<Real>     Transform;
            typedef Erl::Vector6  <Real>     Vector6  ;
            typedef Erl::Vector3  <Real>     Vector3  ;
            typedef Eigen::Matrix<Real,Eigen::Dynamic,1,Eigen::ColMajor> VectorJ    ;
            typedef Eigen::Matrix<Real,6,Eigen::Dynamic,Eigen::ColMajor> MatJacobian;
            typedef types::RobotA<Real>      RobotA;

            Robot();
            Robot(const Robot& _R);
            Robot(const size_t& _MaxSample, const Transform& _InitTrafo,
                  const std::vector<std::pair<SecCC,size_t> > & _C0,
                  const std::vector<std::pair<SecVC,size_t> > & _C1);
            Robot(const size_t& _MaxSample, const Transform& _InitTrafo,
                  const std::deque<std::pair<SecCC,size_t> > & _C0,
                  const std::deque<std::pair<SecVC,size_t> > & _C1);
            Robot(const size_t& _MaxSample, const Transform& _InitTrafo,
                  const std::list<std::pair<SecCC,size_t> > & _C0,
                  const std::list<std::pair<SecVC,size_t> > & _C1);
            Robot(const size_t& _MaxSample, const std::string& _XmlFilename);
            ~Robot();

            bool init() const;
            Robot& operator = (const Robot& _R);

            Transform calcKinematic(const VectorJ& _JVal, const Real& _ArcLengthStep);
            Transform calcKinematic(const VectorJ& _JVal, const Real& _ArcLengthStep,Real& _ArcLengthStepReturn);
            Transform calcKinematic(const VectorJ& _JVal, const size_t& _NSamples);
            Transform calcKinematic(const VectorJ& _JVal, const size_t& _NSamples, Real& _ArcLengthStepReturn);


            void      calcKinematic(Real * const _JVal, const Real& _ArcLengthStep, Real* _Tr);
            bool      calcStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep);
            Real      calcDegreeStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep, const Real& _MinDegreeStability = Real(-1.));

            void      calcJacobian(const VectorJ& _JVal, const Real& _ArcLengthStep, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip);
            void      calcJacobian(const VectorJ& _JVal, const size_t& _NSamples,    const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample);

            void      calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const size_t& _NSamples, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip);
            void      calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const Transform& _TiSample, const size_t& _NSamples,    const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample);

            size_t  getNSections      (            ) const ;
            size_t  getNTubes         (            ) const ;
            size_t  getNJointValues   (            ) const ;
            bool    getStable         (            ) const noexcept;
            Real    getMaxLength      (            ) const noexcept;
            size_t  getMemoryAlignment(            ) const noexcept;
            VectorJ getJointValues    (            ) const noexcept;
            void    getJointValues    (VectorJ& _JV) const noexcept;

            size_t     getMaxSamples       ()                 const noexcept ;
            size_t     getNSamples         ()                 const noexcept ;
            Real       getRadiusSample     (unsigned _Sample) const noexcept ;
            Transform  getTransformSample  (unsigned _Sample) const noexcept ;
            Vector3    getTranslationSample(unsigned _Sample) const noexcept ;
            Real       getTipRadius        ()                 const noexcept ;
            Transform  getTipTransform     ()                 const noexcept ;
            Vector3    getTipTranslation   ()                 const noexcept ;


            void getSection(unsigned _iSection, const SecCC ** _Section) const noexcept;
            void getSection(unsigned _iSection, const SecVC ** _Section) const noexcept;

            bool    setMaxSample(const size_t _MaxSample);
            template<class _Real>
            friend  std::ostream&  operator<< (std::ostream& os, const Robot<_Real> &R);
            operator bool() const ;


            void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn) const;
            void   getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn, const TransRotScale<Real>& _TS) const;




            #ifdef WITH_CTROPENCL
                size_t copy2opencl(void *_robot_cl, int _mode) const;
            #endif



        protected:
            void     deleteRobot();
            RobotA * m_Robot;

        private:

    };
}



#endif // ROBOT_H
