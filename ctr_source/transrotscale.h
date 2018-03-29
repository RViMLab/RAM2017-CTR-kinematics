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

#ifndef TRANSROTSCALE_H
#define TRANSROTSCALE_H

#include <Erl/Realtime/spinlock.h>
#include <cstdlib>
#include "ctr_robot.h"
#include "section_i.h"

namespace CTR
{
template<class Real>
class TransRotScale
{
public:
    virtual ~TransRotScale(){}
    virtual void     set_scale (const Real& _yfactor)         = 0;
    virtual void     add_scale (const Real& _yfactor)         = 0;
    virtual void     next_scale()                             = 0;
    virtual size_t   size_scalefactor()                 const = 0;
    virtual Real     get_scale   ()                     const = 0;
    virtual long     get_iter_position()                const = 0;
    virtual Real     scaleTrans  (const Real&   _value) const = 0;
    virtual Real     scaleRot    (const Real&   _value) const = 0;

};
template <class Real, int LookUpSize>
class GammaCorrection : public TransRotScale<Real>
{
    typedef Eigen::Array <Real,LookUpSize,1>           LU;
    typedef std::vector<Real>                          Vec;
    typedef typename std::vector<Real>::const_iterator VecIter;
    public:
        inline GammaCorrection()
            : m_CurrentScaleFactor(1.)
            , m_ScaleFactor()
            , m_ScaleFactorIter(m_ScaleFactor.begin())
            , m_LookUp     ()
        { Erl::static_if<(LookUpSize>0)>()([&]{m_LookUp.matrix().setLinSpaced(LookUpSize,Real(0.),Real(1.));});}
        inline void set_scale(const Real& _scalefactor)
        {
            //std::lock_guard<Erl::RT::Spinlock> c_lock(m_sl_Lock);
            m_CurrentScaleFactor = _scalefactor;
            Erl::static_if<(LookUpSize>0)>()([&]{
                for(int i(0);i<LookUpSize;i++)
                { m_LookUp[i]= std::pow(m_LookUp[i],_scalefactor); };
            });
        }
        inline void add_scale(const Real& _scalefactor)
        {
            //std::lock_guard<Erl::RT::Spinlock> c_lock(m_sl_Lock);
            if(m_ScaleFactor.empty())
            {
                m_CurrentScaleFactor = _scalefactor;
                m_ScaleFactor.push_back(m_CurrentScaleFactor);
                m_ScaleFactorIter = m_ScaleFactor.begin();
            }
            else
            {
                auto c_BeginPre = std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter);
                m_ScaleFactor.push_back(_scalefactor);
                if(c_BeginPre != std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter))
                {
                    m_ScaleFactorIter = m_ScaleFactor.begin()+c_BeginPre;
                }

            }
        }
        inline void next_scale()
        {
            //std::lock_guard<Erl::RT::Spinlock> c_lock(m_sl_Lock);
            m_ScaleFactorIter++;
            set_scale(*m_ScaleFactorIter);
        }
        inline long get_iter_position() const
        {
            //std::lock_guard<Erl::RT::Spinlock> c_lock(m_sl_Lock);
            return std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter);
        }
        inline size_t size_scalefactor() const
        {
            //std::lock_guard<Erl::RT::Spinlock> c_lock(m_sl_Lock);
            return m_ScaleFactor.size();
        }
        inline Real   get_scale() const
        {
            //std::lock_guard<Erl::RT::Spinlock> c_lock(m_sl_Lock);
            return m_CurrentScaleFactor;
        }
        inline Real scaleTrans(const Real& _value) const
        {
            ERL_ASSERT(_value>=Real(0) && _value<=Real(1));
            return
            Erl::static_if_else<(LookUpSize>0)>()([&]
            {
                //const GammaCorrection& c_const_Ref = GammaCorrection::const_Ref();
                int c_I0 = std::floor(_value*(LookUpSize-1));
                Real c_Mod = _value*(LookUpSize-1)-c_I0;
                Real c_F0 = m_LookUp[c_I0];
                Real c_F1 = m_LookUp[std::ceil (_value*(LookUpSize-1))];
                return c_F1*c_Mod+c_F0*(Real(1)-c_Mod);
            },[&]
            {
                return std::pow(_value,m_CurrentScaleFactor);
            });
        }
        inline Real scaleRot(const Real& _value) const
        {
            ERL_ASSERT(_value>=Real(0) && _value<=Real(1));
            return (_value-Real(0.5));
        }
    private:
        Real                      m_CurrentScaleFactor    ;
        Vec                       m_ScaleFactor    ;
        VecIter                   m_ScaleFactorIter;
        LU                        m_LookUp         ;
    public:
        mutable Erl::RT::Spinlock m_sl_Lock;
};

template <class Real>
class ProportionalScale : public TransRotScale<Real>
{
    typedef Eigen::Matrix<Real,Eigen::Dynamic,1> SLVector;
    typedef Eigen::Matrix<Real,1,2>              PScales;
    typedef std::vector<PScales>                 Vec;
    typedef typename Vec::const_iterator         VecIter;
    public:
        inline ProportionalScale()
            : m_CurrentScaleFactor()
            , m_ScaleFactor()
            , m_ScaleFactorIter(m_ScaleFactor.begin())
        {m_CurrentScaleFactor<<Real(1.),Real(1.);}
        inline bool initialize(const Robot<Real>& _Robot)
        {
            if(!_Robot.init())
            {return false;}
            size_t c_NSection  = _Robot.getNSections();
            m_IterSection      = 0;
            m_OverallMaxLength = Real(0);
            m_SegmentLengths.resize(c_NSection,1);

            const ::CTR::types::SecVC<Real> *c_SecVC;
            const ::CTR::types::SecCC<Real> *c_SecCC;
            for(size_t iS(0);iS<c_NSection;iS++)
            {
                _Robot.getSection(iS,&c_SecCC);
                _Robot.getSection(iS,&c_SecVC);
                ERL_ASSERT( ((c_SecVC!=nullptr) || (c_SecCC!=nullptr)) );
                if(c_SecCC!=nullptr)
                { m_SegmentLengths(iS,0) = c_SecCC->getLength();}
                else
                { m_SegmentLengths(iS,0) = c_SecVC->getLength();}
                m_OverallMaxLength += m_SegmentLengths(iS,0);
            }
            m_CurrentScaleFactorInst = ((m_CurrentScaleFactor[1]-m_CurrentScaleFactor[0])*
                    Real(std::rand())/Real(RAND_MAX)+m_CurrentScaleFactor[0]);
            m_GoalLengthRemain = m_OverallMaxLength*m_CurrentScaleFactorInst;
            m_MaxLengthRemain  = m_OverallMaxLength;
            return true;
        }
        inline void set_scale(const PScales& _scalefactor)
        {
            ERL_ASSERT(m_IterSection==0);
            m_CurrentScaleFactor = _scalefactor;
            m_CurrentScaleFactorInst = ((m_CurrentScaleFactor[1]-m_CurrentScaleFactor[0])*
                    Real(std::rand())/Real(RAND_MAX)+m_CurrentScaleFactor[0]);
            m_GoalLengthRemain = m_OverallMaxLength*m_CurrentScaleFactorInst;
        }
        inline void set_scale(const Real& _scalefactor)
        {
            ERL_ASSERT(m_IterSection==0);
            m_CurrentScaleFactor[0] =  _scalefactor;
            m_CurrentScaleFactor[1] =  _scalefactor;
            m_GoalLengthRemain = m_OverallMaxLength*_scalefactor;
        }
        inline void set_scale(const Real& _scalefactor_min,const Real& _scalefactor_max)
        {
            ERL_ASSERT(m_IterSection==0);
            Real c_min = std::min(_scalefactor_min,_scalefactor_max);
            Real c_max = std::max(_scalefactor_min,_scalefactor_max);
            m_CurrentScaleFactor[0] = c_min;
            m_CurrentScaleFactor[1] = c_max;

            m_CurrentScaleFactorInst = ((m_CurrentScaleFactor[1]-m_CurrentScaleFactor[0])*
                    Real(std::rand())/Real(RAND_MAX)+m_CurrentScaleFactor[0]);
            m_GoalLengthRemain = m_OverallMaxLength*m_CurrentScaleFactorInst;

        }
        inline void add_scale(const Real& _scalefactor)
        {
            ERL_ASSERT(m_IterSection==0);
            if(m_ScaleFactor.empty())
            {
                m_CurrentScaleFactor << _scalefactor,_scalefactor;
                m_ScaleFactor.push_back(m_CurrentScaleFactor);
                m_ScaleFactorIter = m_ScaleFactor.begin();
                m_CurrentScaleFactorInst = ((m_CurrentScaleFactor[1]-m_CurrentScaleFactor[0])*
                        Real(std::rand())/Real(RAND_MAX)+m_CurrentScaleFactor[0]);
                m_GoalLengthRemain = m_OverallMaxLength*m_CurrentScaleFactorInst;
            }
            else
            {
                PScales c_CurrentScaleFactor;
                c_CurrentScaleFactor << _scalefactor,_scalefactor;
                auto c_BeginPre = std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter);
                m_ScaleFactor.push_back(c_CurrentScaleFactor);
                if(c_BeginPre != std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter))
                { m_ScaleFactorIter = m_ScaleFactor.begin()+c_BeginPre; }

            }
        }
        inline void add_scale(const Real& _scalefactor_min,const Real& _scalefactor_max)
        {
            ERL_ASSERT(m_IterSection==0);
            if(m_ScaleFactor.empty())
            {
                m_CurrentScaleFactor << std::min(_scalefactor_min,_scalefactor_max),
                                    std::max(_scalefactor_min,_scalefactor_max);
                m_ScaleFactor.push_back(m_CurrentScaleFactor);
                m_ScaleFactorIter = m_ScaleFactor.begin();
                m_CurrentScaleFactorInst = ((m_CurrentScaleFactor[1]-m_CurrentScaleFactor[0])*
                        Real(std::rand())/Real(RAND_MAX)+m_CurrentScaleFactor[0]);
                m_GoalLengthRemain = m_OverallMaxLength*m_CurrentScaleFactorInst;
            }
            else
            {
                PScales c_CurrentScaleFactor;
                c_CurrentScaleFactor << std::min(_scalefactor_min,_scalefactor_max),
                                        std::max(_scalefactor_min,_scalefactor_max);
                auto c_BeginPre = std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter);
                m_ScaleFactor.push_back(c_CurrentScaleFactor);
                if(c_BeginPre != std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter))
                { m_ScaleFactorIter = m_ScaleFactor.begin()+c_BeginPre; }

            }
        }
        inline void next_scale()
        {
            ERL_ASSERT(m_IterSection==0);
            m_ScaleFactorIter++;
            ERL_ASSERT(m_ScaleFactorIter != m_ScaleFactor.end());
            set_scale(*m_ScaleFactorIter);
        }
        inline long get_iter_position() const
        {
            return std::distance<VecIter>(m_ScaleFactor.begin(),m_ScaleFactorIter);
        }
        inline size_t size_scalefactor() const
        {
            return m_ScaleFactor.size();
        }
        inline Real   get_scale() const
        {
            return m_CurrentScaleFactorInst;
        }
        inline Real scaleTrans(const Real& _value) const
        {
            ERL_ASSERT(_value>=Real(0) && _value<=Real(1));
            Real c_val_min = std::max((m_GoalLengthRemain - m_MaxLengthRemain + m_SegmentLengths(m_IterSection,0))/ m_SegmentLengths(m_IterSection,0),Real(0));
            Real c_val_max = std::min((m_GoalLengthRemain / m_SegmentLengths(m_IterSection,0)),Real(1));
            Real c_val     =  (c_val_max-c_val_min)*_value + c_val_min;
            m_GoalLengthRemain -= c_val*m_SegmentLengths(m_IterSection,0);
            m_MaxLengthRemain  -= m_SegmentLengths(m_IterSection,0);
            m_IterSection++;
            if(m_IterSection >= size_t(m_SegmentLengths.rows()) )
            {
                m_IterSection      = 0;
                m_MaxLengthRemain  = m_OverallMaxLength;
                m_CurrentScaleFactorInst = ((m_CurrentScaleFactor[1]-m_CurrentScaleFactor[0])*
                         _value+m_CurrentScaleFactor[0]);
                m_GoalLengthRemain = m_OverallMaxLength*m_CurrentScaleFactorInst;
            }
            return Erl::clamp<Real>(c_val,Real(0.),Real(1.));
        }
        inline Real scaleRot(const Real& _value) const
        {
            ERL_ASSERT(_value>=Real(0) && _value<=Real(1));
            return (_value-Real(0.5));
        }
    private:
        PScales                   m_CurrentScaleFactor;
        mutable Real              m_CurrentScaleFactorInst;
        Real                      m_OverallMaxLength  ;
        mutable Real              m_GoalLengthRemain  ;
        mutable Real              m_MaxLengthRemain   ;
        mutable size_t            m_IterSection       ;

        Vec                       m_ScaleFactor       ;
        VecIter                   m_ScaleFactorIter   ;
        SLVector                  m_SegmentLengths    ;

};

}


#endif // TRANSROTSCALE_H
