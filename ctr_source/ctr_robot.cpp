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

#include "ctr_robot.h"
#include "robot_i.h"
#include "section_i.h"
#include "ctr_static.h"


namespace CTR
{
    template<class Real>
    Robot<Real>::Robot()
        : m_Robot(nullptr)
    { }
    template<class Real>
    Robot<Real>::Robot(const Robot& _R)
        : m_Robot(details::constructRobot(_R.m_Robot))
    { }
    template<class Real>
    Robot<Real>::~Robot()
    { deleteRobot(); }
    template<class Real>void
    Robot<Real>:: deleteRobot()
    {
        if(m_Robot != nullptr)
        {
            if(m_Robot->getMaxSamples() > 0)
            { Erl::aligned_free(m_Robot); }
            else{ delete m_Robot; }

            m_Robot = nullptr;
        }
    }

    template<class Real> Robot<Real>&
    Robot<Real>::operator = (const Robot& _R)
    {
        if(_R.m_Robot==m_Robot)
        { return (*this); }
        if(_R.m_Robot==nullptr)
        { deleteRobot(); return (*this); }
        size_t c_Alignment;
        size_t c_MemReqOther =  (_R.m_Robot==nullptr) ? size_t(0) : _R.m_Robot->requiredMemory(c_Alignment);
        size_t c_MemReqThis  =  (   m_Robot==nullptr) ? size_t(0) :    m_Robot->requiredMemory();
        if( c_MemReqOther != c_MemReqThis)
        {
            deleteRobot();
            m_Robot = reinterpret_cast<RobotA*>(Erl::aligned_malloc(c_MemReqOther,c_Alignment));
        }
        m_Robot = _R.m_Robot->move_into(c_MemReqOther,m_Robot);
        return (*this);
    }
    template<class Real>
    Robot<Real>::Robot(const size_t& _MaxSample, const Transform& _InitTrafo,
                       const std::vector<std::pair<SecCC,size_t> > & _C0,
                       const std::vector<std::pair<SecVC,size_t> > & _C1)
        : m_Robot(details::constructRobot<RobotA>(_MaxSample,_InitTrafo,_C0,_C1))
    { }
    template<class Real>
    Robot<Real>::Robot(const size_t& _MaxSample, const Transform& _InitTrafo,
                       const std::deque<std::pair<SecCC,size_t> > & _C0,
                       const std::deque<std::pair<SecVC,size_t> > & _C1)
        : m_Robot(details::constructRobot<RobotA>(_MaxSample,_InitTrafo,_C0,_C1))
    { }
    template<class Real>
    Robot<Real>::Robot(const size_t& _MaxSample, const Transform& _InitTrafo,
                       const std::list<std::pair<SecCC,size_t> > & _C0,
                       const std::list<std::pair<SecVC,size_t> > & _C1)
        : m_Robot(details::constructRobot<RobotA>(_MaxSample,_InitTrafo,_C0,_C1))
    { }
    template<class Real>
    bool Robot<Real>::init() const
    { return (m_Robot != nullptr); }

    template<class Real> typename Robot<Real>::Transform
    Robot<Real>::calcKinematic(const VectorJ& _JVal, const Real& _ArcLengthStep, Real& _ArcLengthStepReturn)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->calcKinematic(_JVal,_ArcLengthStep,_ArcLengthStepReturn);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->calcKinematic(_JVal,_ArcLengthStep,_ArcLengthStepReturn);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Erl::Transform<Real>::Identity();
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Transform
    Robot<Real>::calcKinematic(const VectorJ& _JVal, const Real& _ArcLengthStep)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->calcKinematic(_JVal,_ArcLengthStep);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->calcKinematic(_JVal,_ArcLengthStep);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Erl::Transform<Real>::Identity();
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Transform
    Robot<Real>::calcKinematic(const VectorJ& _JVal, const size_t& _NSamples,Real& _ArcLengthStepReturn)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->calcKinematic(_JVal,_NSamples,_ArcLengthStepReturn);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->calcKinematic(_JVal,_NSamples,_ArcLengthStepReturn);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Erl::Transform<Real>::Identity();
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Transform
    Robot<Real>::calcKinematic(const VectorJ& _JVal, const size_t& _NSamples)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->calcKinematic(_JVal,_NSamples);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->calcKinematic(_JVal,_NSamples);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Erl::Transform<Real>::Identity();
        }
        #endif
    }
    template<class Real> void
    Robot<Real>::calcKinematic(Real* const _JVal, const Real& _ArcLengthStep, Real* _Tr)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            Matrix<Real,Eigen::Dynamic,1> c_JVal(m_Robot->getNJoints(),1,_JVal);
            m_Robot->calcKinematic(c_JVal,_ArcLengthStep).getData(_Tr,Erl::ColMajor);
        #else
        if(m_Robot!=nullptr)
        {
            Matrix<Real,Eigen::Dynamic,1> c_JVal(m_Robot->getNJoints(),1,_JVal);
            m_Robot->calcKinematic(c_JVal,_ArcLengthStep).getData(_Tr,Erl::ColMajor);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return ;
        }
        #endif
    }

    template<class Real> void
    Robot<Real>::calcJacobian(const VectorJ& _JVal, const Real& _ArcLengthStep, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            m_Robot->calcJacobian(_JVal,_ArcLengthStep,_FiniteDiff_Trans,_FiniteDiff_Rot,_JacobianTip);
        #else
        if(m_Robot!=nullptr)
        {
            m_Robot->calcJacobian(_JVal,_ArcLengthStep,_FiniteDiff_Trans,_FiniteDiff_Rot,_JacobianTip);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return ;
        }
        #endif
    }

    template<class Real> void
    Robot<Real>::calcJacobian(const VectorJ& _JVal, const size_t& _NSamples, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            m_Robot->calcJacobian(_JVal,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_iSample,_JacobianTip,_JacobianiSample);
        #else
        if(m_Robot!=nullptr)
        {
            m_Robot->calcJacobian(_JVal,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_iSample,_JacobianTip,_JacobianiSample);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return ;
        }
        #endif
    }

    template<class Real> void
    Robot<Real>::calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const size_t& _NSamples, const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, MatJacobian& _JacobianTip)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            m_Robot->calcJacobian(_JVal,_TTip,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_JacobianTip);
        #else
        if(m_Robot!=nullptr)
        {
            m_Robot->calcJacobian(_JVal,_TTip,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_JacobianTip);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return ;
        }
        #endif
    }

    template<class Real> void
    Robot<Real>::calcJacobian(const VectorJ& _JVal, const Transform& _TTip, const Transform& _TiSample, const size_t& _NSamples,    const Real& _FiniteDiff_Trans, const Real& _FiniteDiff_Rot, const size_t& _iSample, MatJacobian& _JacobianTip, MatJacobian& _JacobianiSample)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            m_Robot->calcJacobian(_JVal,_TTip,_TiSample,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_iSample,_JacobianTip,_JacobianiSample);
        #else
        if(m_Robot!=nullptr)
        {
            m_Robot->calcJacobian(_JVal,_TTip,_TiSample,_NSamples,_FiniteDiff_Trans,_FiniteDiff_Rot,_iSample,_JacobianTip,_JacobianiSample);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return ;
        }
        #endif
    }

    template<class Real> bool
    Robot<Real>::setMaxSample(const size_t _MaxSample)
    {
        if((_MaxSample>0) && (m_Robot != nullptr) )
        {
            RobotA * c_Robot  ;
            size_t c_Alignment;
            auto   c_MemorySize = m_Robot->requiredMemory(_MaxSample,c_Alignment);
            void * c_Memory     = Erl::aligned_malloc(c_MemorySize,c_Alignment);
            if(c_Memory != nullptr)
            {
                c_Robot = m_Robot->move_into(_MaxSample,c_MemorySize,c_Memory);
                Erl::aligned_free(m_Robot);
                m_Robot = c_Robot;
                return true;
            }
        }
        return false;
    }
    template<class Real>
    Robot<Real>::Robot(const size_t& _MaxSample, const std::string& _XmlFilename)
        : m_Robot(details::readConcentricTubefromXML<Real>(_XmlFilename,_MaxSample,false))
    { }
    template<class Real>
    size_t Robot<Real>::getNSections()   const
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getNSections();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getNSections();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real>
    size_t Robot<Real>::getNTubes()      const
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getNTubes();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getNTubes();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real>
    size_t Robot<Real>::getNJointValues()      const
    { return getNTubes()+getNSections()+1; }

    template<class Real>              size_t
    Robot<Real>::getMaxSamples       ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getMaxSamples();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getMaxSamples();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real>              size_t
    Robot<Real>::getNSamples         ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getNSamples();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getNSamples();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real>              Real
    Robot<Real>::getMaxLength         ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getMaxLength();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getMaxLength();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real>              size_t
    Robot<Real>::getMemoryAlignment         ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getMemoryAlignment();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getMemoryAlignment();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real> typename Robot<Real>::VectorJ
    Robot<Real>::getJointValues    (            ) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getJointValues();
        #else
            if(m_Robot!=nullptr)
            { return m_Robot->getJointValues();}
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Robot<Real>::VectorJ();
        #endif
    }
    template<class Real> void
    Robot<Real>::getJointValues    (VectorJ& _JV) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            m_Robot->getJointValues(_JV);
            return ;
        #else
            if(m_Robot!=nullptr)
            { m_Robot->getJointValues(_JV);
              return;
            }
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return;
        #endif
    }

    template<class Real>              bool
    Robot<Real>::getStable         ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getStable();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getStable();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
        }
        return false;
        #endif
    }
    template<class Real>              Real
    Robot<Real>::getRadiusSample     (unsigned _Sample) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getRadiusSample(_Sample);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getRadiusSample(_Sample);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Transform
    Robot<Real>::getTransformSample  (unsigned _Sample) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getTransformSample(_Sample);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getTransformSample(_Sample);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Transform ::Zero();;
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Vector3
    Robot<Real>::getTranslationSample(unsigned _Sample) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getTranslationSample(_Sample);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getTranslationSample(_Sample);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Vector3 ::Zero();;
        }
        #endif
    }
    template<class Real>              Real
    Robot<Real>::getTipRadius        ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getTipRadius();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getTipRadius();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Transform
    Robot<Real>::getTipTransform     ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getTipTransform();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getTipTransform();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Transform::Zero();
        }
        #endif
    }
    template<class Real> typename Robot<Real>::Vector3
    Robot<Real>::getTipTranslation   ()                 const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getTipTranslation();
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getTipTranslation();
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return Vector3 ::Zero();
        }
        #endif
    }
    template<class Real>
    void Robot<Real>::getSection(unsigned _iSection, const Robot<Real>::SecCC ** _Section) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getSection(_iSection,_Section);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getSection(_iSection,_Section);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            *_Section = nullptr;
            return;
        }
        #endif
    }
    template<class Real>
    void Robot<Real>::getSection(unsigned _iSection, const Robot<Real>::SecVC ** _Section) const noexcept
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getSection(_iSection,_Section);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->getSection(_iSection,_Section);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            *_Section = nullptr;
            return;
        }
        #endif
    }

    template<class Real>
    Robot<Real>::operator bool() const
    { return (m_Robot!=nullptr); }

    template<class _Real>
    std::ostream&  operator<< (std::ostream& os, const Robot<_Real> &R)
    {
        if(R.m_Robot!=nullptr)
        {
            os<<(*R.m_Robot);
        }
        else
        {
            #if CTR_DEBUG_LEVEL >= CTR_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
        }
        return os;
    }
    template<class Real>
    void   Robot<Real>::getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn) const
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->getRandomJointValues(_RNG,_RND,_JValReturn);
        #else
        if(m_Robot!=nullptr)
        {
            m_Robot->getRandomJointValues(_RNG,_RND,_JValReturn);
            return;
        }
        #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
        Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
        #endif
        return;
        #endif
    }
    template<class Real>
    void   Robot<Real>::getRandomJointValues(CT_RNG& _RNG, CT_RND<Real>& _RND  , Eigen::Ref<VectorJ > _JValReturn, const TransRotScale<Real>& _TS) const
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            m_Robot->getRandomJointValues(_RNG,_RND,_JValReturn,_TS);
            return;
        #else
        if(m_Robot!=nullptr)
        {
            m_Robot->getRandomJointValues(_RNG,_RND,_JValReturn,_TS);
            return;
        }
        #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
        Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
        #endif
        return;
        #endif
    }
    template<class Real>
    bool   Robot<Real>::calcStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep)
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->calcStability(_JVal,_AngleDiscretizationStep,_ArcLengthStep);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->calcStability(_JVal,_AngleDiscretizationStep,_ArcLengthStep);
        }
        #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
        Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
        #endif
        return false;
        #endif
    }
    template<class Real>
    Real   Robot<Real>::calcDegreeStability(const VectorJ& _JVal, const Real& _AngleDiscretizationStep, const Real& _ArcLengthStep, const Real& _MinDegreeStability )
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->calcDegreeStability(_JVal,_AngleDiscretizationStep,_ArcLengthStep,_MinDegreeStability);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->calcDegreeStability(_JVal,_AngleDiscretizationStep,_ArcLengthStep,_MinDegreeStability);
        }
        #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
        Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
        #endif
        return std::numeric_limits<Real>::quiet_NaN();
        #endif
    }


    #ifdef WITH_CTROPENCL
    template<class Real>
    size_t Robot<Real>::copy2opencl(void *_robot_cl, int _mode) const
    {
        #if CT_DEBUG_LEVEL < CT_DEBUG_NORMAL
            return m_Robot->copy2opencl(_robot_cl,_mode);
        #else
        if(m_Robot!=nullptr)
        {
            return m_Robot->copy2opencl(_robot_cl,_mode);
        }
        else
        {
            #if CT_DEBUG_LEVEL >= CT_DEBUG_NOISY
            Erl::debug_cout_e()<<"FOOL: no robot, FUNC:"<<FUNC_STRING<<std::endl;
            #endif
            return 0;
        }
        #endif
    }
    #endif
    #ifdef CT_ROBOT_FLOAT
        #ifdef  CTR_REAL
            #undef  CTR_REAL
        #endif
        #define CTR_REAL float
        template class Robot<CTR_REAL >;
        template std::ostream& operator<<(std::ostream&, const Robot<CTR_REAL>& );
    #endif
    #ifdef CT_ROBOT_DOUBLE
        #ifdef  CTR_REAL
            #undef  CTR_REAL
        #endif
        #define CTR_REAL double
        template class Robot<CTR_REAL>;
        template std::ostream& operator<<(std::ostream&, const Robot<CTR_REAL>& );
    #endif
    #ifdef CT_ROBOT_LONG_DOUBLE
        #ifdef  CTR_REAL
            #undef  CTR_REAL
        #endif
        #define CTR_REAL long double
        template class Robot<CTR_REAL>;
        template std::ostream& operator<<(std::ostream&, const Robot<CTR_REAL>& );
    #endif


}
