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

#ifndef CTR_ROBOT_STATIC_H
#define CTR_ROBOT_STATIC_H

#include <vector>
#include <array>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <Erl/Utility/debug.h>

#include <random>

#include "ctr_robot.h"
#include "section_i.h"
#include "robot_i.h"

namespace CTR
{
namespace details
{
template<class RobotA, class ...ContArgs>
RobotA * constructRobot(const size_t& _MaxSample, const Erl::Transform<typename RobotA::RealRBase>& _InitTrafo, const ContArgs&... _Cont)
{
    RobotA * c_Robot(nullptr);
    RobotA * c_InitRobot = createRobot<RobotA>(_InitTrafo,_Cont...);
    if((c_InitRobot!=nullptr) && (_MaxSample>0))
    {
        size_t c_Alignment;
        auto   c_MemorySize = c_InitRobot->requiredMemory(_MaxSample,c_Alignment);
        void * c_Memory     = Erl::aligned_malloc(c_MemorySize,c_Alignment);
        if(c_Memory != nullptr){  c_Robot = c_InitRobot->move_into(_MaxSample,c_MemorySize,c_Memory); }
        delete c_InitRobot;
        //Erl::free_aligned_init_array(c_InitRobot,size_t(1));
    }
    return c_Robot;
}
template<class RobotA>
RobotA * constructRobot2Ptr(const size_t& _MaxSample, const Erl::Transform<typename RobotA::RealRBase>& _InitTrafo,
                            const SectionCC<typename RobotA::RealRBase>* _SecCC, const std::vector<size_t>& _SecCC_I,
                            const SectionVC<typename RobotA::RealRBase>* _SecVC, const std::vector<size_t>& _SecVC_It)
{
    RobotA * c_Robot(nullptr);
    RobotA * c_InitRobot = createRobot2Ptr<RobotA>(_InitTrafo,_SecCC,_SecCC_I,_SecVC,_SecVC_It);
    if((c_InitRobot!=nullptr) && (_MaxSample>0))
    {
        size_t c_Alignment;
        auto   c_MemorySize = c_InitRobot->requiredMemory(_MaxSample,c_Alignment);
        void * c_Memory     = Erl::aligned_malloc(c_MemorySize,c_Alignment);
        if(c_Memory != nullptr){  c_Robot = c_InitRobot->move_into(_MaxSample,c_MemorySize,c_Memory); }
        delete c_InitRobot;
        //Erl::free_aligned_init_array(c_InitRobot,size_t(1));
    }
    return c_Robot;
}
template<class RobotA>
RobotA * constructRobot(const RobotA* const _Robot)
{
    if( (_Robot!=nullptr) )
    {
        size_t c_Alignment;
        size_t c_MaxSamples = _Robot->getMaxSamples();
        auto   c_MemorySize = _Robot->requiredMemory(c_MaxSamples,c_Alignment);
        void * c_Memory     = Erl::aligned_malloc(c_MemorySize,c_Alignment);
        if(c_Memory != nullptr){  return _Robot->move_into(c_MaxSamples,c_MemorySize,c_Memory); }
    }
    return nullptr;
}
#define CTR_USE_PTR_FOR_INIT
#ifdef  CTR_USE_PTR_FOR_INIT
template<class Real>
types::RobotA<Real>* readConcentricTubefromXMLV1(const std::string& _filename, const size_t& _maxSamples, bool _verbose)
{
    size_t c_ISections(0);
    size_t c_ITube(0);
    size_t c_IJv(1);
    types::RobotA<Real> * R(nullptr);

    std::ifstream c_XMLFile(_filename.c_str(),std::ifstream::in);
    if (!c_XMLFile.good())
    {
        std::cerr<<"ERROR: concentricTubefromXML::XMLParse, could not read file: "<<_filename<<std::endl;
        return nullptr;
    }

    using boost::property_tree::ptree;
    ptree c_PropertyTree;
    try
    {
        read_xml(c_XMLFile, c_PropertyTree, boost::property_tree::xml_parser::trim_whitespace);
    }
    catch (std::exception e)
    {
        Erl::debug_cout_e() << "ERROR: " << FILE_STRING << ", Exception: " << e.what() <<", File: "<< _filename.c_str() <<std::endl;
        return nullptr;
    }
    c_XMLFile.close();

    Erl::Transform<Real>  c_InitTrafo = Erl::Transform<Real>::Identity();
    Eigen::Matrix<Real,Eigen::Dynamic,1> c_JVal
            = Eigen::Matrix<Real,Eigen::Dynamic,1>::Zero(1+3*CTR_MAXSECTION,1);

    std::vector<size_t>         c_VectorCC_I;
    std::vector<size_t>         c_VectorVC_I;
    std::vector< std::string >  c_SectionName;
    c_VectorCC_I .reserve(CTR_MAXSECTION);
    c_VectorVC_I .reserve(CTR_MAXSECTION);
    c_SectionName.reserve(CTR_MAXSECTION);

    types::SecCC<Real>* c_VectorCC = Erl::alloc_aligned_init_array_noparameter<types::SecCC<Real> >(CTR_MAXSECTION,alignof(types::SecCC<Real>));
    types::SecVC<Real>* c_VectorVC = Erl::alloc_aligned_init_array_noparameter<types::SecVC<Real> >(CTR_MAXSECTION,alignof(types::SecVC<Real>));

    try
    {
        for( auto v : c_PropertyTree.get_child("robot") )
        {
            std::string c_v_first = v.first;
            std::transform(c_v_first.begin(), c_v_first.end(), c_v_first.begin(), ::tolower);


            if( c_v_first.compare("section")==0 )
            {
                std::array<Real,2> c_stiffness     ;
                std::array<Real,2> c_curvature     ;
                std::array<Real,2> c_radius        ;
                std::array<Real,2> c_possisionRatio;
                std::array<Real,2> c_phi           ;
                std::array<Real,2> c_length        ;

                unsigned c_IdxTube(0);

                for (auto v2 : v.second)
                {


                    std::string c_v2_first = v2.first;
                    std::transform(c_v2_first.begin(), c_v2_first.end(), c_v2_first.begin(), ::tolower);
                    if( c_v2_first.compare("name")==0 )
                    {
                        std::string c_SNameStr = v2.second.data().c_str();
                        c_SectionName.push_back(c_SNameStr);
                    }
                    if( c_v2_first.compare("tube")==0 )
                    {
                        if(c_IdxTube>1){throw std::logic_error("VIOLATED: maximum of 2 tubes per section");}
                        c_stiffness     [c_IdxTube] = v2.second.get<Real>("stiffness"   );
                        c_curvature     [c_IdxTube] = v2.second.get<Real>("curvature"   );
                        c_radius        [c_IdxTube] = v2.second.get<Real>("radius"      );
                        c_possisionRatio[c_IdxTube] = v2.second.get<Real>("poissonRatio");
                        c_phi           [c_IdxTube] = v2.second.get<Real>("phi"         ,Real(0));
                        c_length        [c_IdxTube] = v2.second.get<Real>("length"      );
                        c_IdxTube ++; c_ITube ++;
                    }
                }
                if(c_IdxTube==1)
                {
                    types::SecCC<Real> c_Sc({c_stiffness[0]},{c_curvature[0]},c_length[0],{c_radius[0]},{c_possisionRatio[0]});
                    //c_VectorCC.push_back(std::make_pair(c_Sc,c_ISections++));
                    c_VectorCC[c_VectorCC_I.size()] = c_Sc;
                    c_VectorCC_I.push_back(c_ISections++);
                    c_VectorCC_I.back() = (c_ISections-1);
                    c_JVal(c_IJv++,0) = c_phi[0]; c_IJv+=1;
                }
                else if(c_IdxTube==2)
                {
                    types::SecVC<Real> c_Sv(c_stiffness,c_curvature,c_length[0],c_radius,c_possisionRatio);
                    //c_VectorVC.push_back(std::make_pair(c_Sv,c_ISections++));
                    c_VectorVC[c_VectorVC_I.size()] = c_Sv;
                    c_VectorVC_I.push_back(c_ISections++);
                    c_VectorVC_I.back() = (c_ISections-1);
                    c_JVal(c_IJv++,0) = c_phi[0]; c_IJv+=2;
                }

                if(c_ISections>CTR_MAXSECTION)
                {
                    throw std::invalid_argument(std::string("MAX SECTION: This Programm supports only: ")+std::to_string(int(CTR_MAXSECTION))+" sections, please adapt your defines (CTR_MAXSECTION) at compilation.");
                }
            }
            if( c_v_first.compare("entryframe")==0 )
            {
                c_InitTrafo = Erl::Transform<Real>(v.second.get<Real>("x00"),v.second.get<Real>("x01"),v.second.get<Real>("x02"),v.second.get<Real>("x03"),
                                                   v.second.get<Real>("x10"),v.second.get<Real>("x11"),v.second.get<Real>("x12"),v.second.get<Real>("x13"),
                                                   v.second.get<Real>("x20"),v.second.get<Real>("x21"),v.second.get<Real>("x22"),v.second.get<Real>("x23"));

                if(!Erl::Transform<Real>::verifyOrthogonal(c_InitTrafo,0.01))
                {
                    throw std::invalid_argument("InitTrafo: Not a Transformationmatrix");
                }
                c_InitTrafo = Erl::Transform<Real>::nearestOrthogonal(c_InitTrafo);

            }


        }
        //R = constructRobot<types::RobotA<Real>>(_maxSamples,c_InitTrafo,c_VectorCC,c_VectorVC);

        R = constructRobot2Ptr<types::RobotA<Real>>(_maxSamples,c_InitTrafo,c_VectorCC,c_VectorCC_I,c_VectorVC,c_VectorVC_I);

        size_t c_IterCC = 0; //c_VectorCC.begin();
        size_t c_IterVC = 0; //c_VectorVC.begin();

        if(_verbose)
        {
            std::stringstream c_ss;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            c_ss<<"File: "<<_filename<<std::endl;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            for(unsigned i(0);i<c_ISections;i++)
            {
                if(i<c_SectionName.size())
                { c_ss<<c_SectionName.at(i)<<std::endl; }

                if(c_IterCC != c_VectorCC_I.size())
                { if(c_VectorCC_I[c_IterCC] == i){ c_ss<<c_VectorCC[c_IterCC]; c_IterCC++;} }
                if(c_IterVC != c_VectorVC_I.size())
                { if(c_VectorVC_I[c_IterVC] == i){ c_ss<<c_VectorVC[c_IterVC]; c_IterVC++;} }
                c_ss<<std::endl;
            }
            c_ss<<"InitialFrame: "<<std::endl;
            Erl::debug_set_io_format(c_ss,4,7,std::ios_base::fixed,true);
            c_ss<<c_InitTrafo<<std::flush;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            Erl::debug_cout()<<c_ss.str();
        }
    }
    catch(std::exception& e)
    {
        Erl::debug_cout_e()<<"Error: readConcentricTubefromXMLV1::read, parsing failed: "<<e.what()<<"; File: "<<_filename<<std::endl;
        assert(false && "Error in Construction of Tool");
    }
    if(R != nullptr)
    {
        c_JVal.conservativeResize(c_IJv,1);
        R->calcKinematic(c_JVal, (R->getMaxLength() / R->getMaxSamples()) );
    }
    if(c_VectorCC!=nullptr)
    { Erl::free_aligned_init_array(c_VectorCC,CTR_MAXSECTION); c_VectorCC = nullptr;}
    if(c_VectorVC!=nullptr)
    { Erl::free_aligned_init_array(c_VectorVC,CTR_MAXSECTION); c_VectorVC = nullptr;}

    return R;
}
#else
template<class Real>
types::RobotA<Real>* readConcentricTubefromXMLV1(const std::string& _filename, const size_t& _maxSamples, bool _verbose)
{
    size_t c_ISections(0);
    size_t c_ITube(0);
    size_t c_IJv(1);
    types::RobotA<Real> * R(nullptr);

    std::ifstream c_XMLFile(_filename.c_str(),std::ifstream::in);
    if (!c_XMLFile.good())
    {
        std::cerr<<"ERROR: concentricTubefromXML::XMLParse, could not read file: "<<_filename<<std::endl;
        return nullptr;
    }

    using boost::property_tree::ptree;
    ptree c_PropertyTree;
    try
    {
        read_xml(c_XMLFile, c_PropertyTree, boost::property_tree::xml_parser::trim_whitespace);
    }
    catch (std::exception e)
    {
        Erl::debug_cout_e() << "ERROR: " << FILE_STRING << ", Exception: " << e.what() <<", File: "<< _filename.c_str() <<std::endl;
        return nullptr;
    }
    c_XMLFile.close();

    Erl::Transform<Real>  c_InitTrafo = Erl::Transform<Real>::Identity();
    Eigen::Matrix<Real,Eigen::Dynamic,1> c_JVal
            = Eigen::Matrix<Real,Eigen::Dynamic,1>::Zero(1+3*CTR_MAXSECTION,1);

    std::vector<std::pair<types::SecCC<Real>,size_t > > c_VectorCC;
    std::vector<std::pair<types::SecVC<Real>,size_t > > c_VectorVC;
    std::vector< std::string >                          c_SectionName;
    c_VectorCC   .reserve(CTR_MAXSECTION);
    c_VectorVC   .reserve(CTR_MAXSECTION);
    c_SectionName.reserve(CTR_MAXSECTION);

    try
    {
        for( auto v : c_PropertyTree.get_child("robot") )
        {
            std::string c_v_first = v.first;
            std::transform(c_v_first.begin(), c_v_first.end(), c_v_first.begin(), ::tolower);


            if( c_v_first.compare("section")==0 )
            {
                std::array<Real,2> c_stiffness     ;
                std::array<Real,2> c_curvature     ;
                std::array<Real,2> c_radius        ;
                std::array<Real,2> c_possisionRatio;
                std::array<Real,2> c_phi           ;
                std::array<Real,2> c_length        ;

                unsigned c_IdxTube(0);

                for (auto v2 : v.second)
                {


                    std::string c_v2_first = v2.first;
                    std::transform(c_v2_first.begin(), c_v2_first.end(), c_v2_first.begin(), ::tolower);
                    if( c_v2_first.compare("name")==0 )
                    {
                        c_SectionName.push_back(v2.second.data().c_str());
                    }
                    if( c_v2_first.compare("tube")==0 )
                    {
                        if(c_IdxTube>1){throw std::logic_error("VIOLATED: maximum of 2 tubes per section");}
                        c_stiffness     [c_IdxTube] = v2.second.get<Real>("stiffness"   );
                        c_curvature     [c_IdxTube] = v2.second.get<Real>("curvature"   );
                        c_radius        [c_IdxTube] = v2.second.get<Real>("radius"      );
                        c_possisionRatio[c_IdxTube] = v2.second.get<Real>("poissonRatio");
                        c_phi           [c_IdxTube] = v2.second.get<Real>("phi"         ,Real(0));
                        c_length        [c_IdxTube] = v2.second.get<Real>("length"      );
                        c_IdxTube ++; c_ITube ++;
                    }
                }
                if(c_IdxTube==1)
                {
                    types::SecCC<Real> c_Sc({c_stiffness[0]},{c_curvature[0]},c_length[0],{c_radius[0]},{c_possisionRatio[0]});
                    c_VectorCC.push_back(std::make_pair(c_Sc,c_ISections++));
                    c_JVal(c_IJv++,0) = c_phi[0]; c_IJv+=1;
                }
                else if(c_IdxTube==2)
                {
                    types::SecVC<Real> c_Sv(c_stiffness,c_curvature,c_length[0],c_radius,c_possisionRatio);
                    c_VectorVC.push_back(std::make_pair(c_Sv,c_ISections++));
                    c_JVal(c_IJv++,0) = c_phi[0]; c_IJv+=2;
                }

                if(c_ISections>CTR_MAXSECTION)
                {
                    throw std::invalid_argument(std::string("MAX SECTION: This Programm supports only: ")+std::to_string(int(CTR_MAXSECTION))+" sections, please adapt your defines (CTR_MAXSECTION) at compilation.");
                }
            }
            if( c_v_first.compare("entryframe")==0 )
            {
                c_InitTrafo = Erl::Transform<Real>(v.second.get<Real>("x00"),v.second.get<Real>("x01"),v.second.get<Real>("x02"),v.second.get<Real>("x03"),
                                                   v.second.get<Real>("x10"),v.second.get<Real>("x11"),v.second.get<Real>("x12"),v.second.get<Real>("x13"),
                                                   v.second.get<Real>("x20"),v.second.get<Real>("x21"),v.second.get<Real>("x22"),v.second.get<Real>("x23"));

                if(!Erl::Transform<Real>::verifyOrthogonal(c_InitTrafo,0.01))
                {
                    throw std::invalid_argument("InitTrafo: Not a Transformationmatrix");
                }
                c_InitTrafo = Erl::Transform<Real>::nearestOrthogonal(c_InitTrafo);

            }


        }
        R = constructRobot<types::RobotA<Real>>(_maxSamples,c_InitTrafo,c_VectorCC,c_VectorVC);

        auto c_IterCC = c_VectorCC.begin();
        auto c_IterVC = c_VectorVC.begin();

        if(_verbose)
        {
            std::stringstream c_ss;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            c_ss<<"File: "<<_filename<<std::endl;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            for(unsigned i(0);i<c_ISections;i++)
            {
                if(i<c_SectionName.size())
                { c_ss<<c_SectionName.at(i)<<std::endl; }

                if(c_IterCC != c_VectorCC.end())
                { if(c_IterCC->second == i){ c_ss<<c_IterCC->first; c_IterCC++;} }
                if(c_IterVC != c_VectorVC.end())
                { if(c_IterVC->second == i){ c_ss<<c_IterVC->first; c_IterVC++;} }
                c_ss<<std::endl;
            }
            c_ss<<"InitialFrame: "<<std::endl;
            Erl::debug_set_io_format(c_ss,4,7,std::ios_base::fixed,true);
            c_ss<<c_InitTrafo<<std::flush;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            Erl::debug_cout()<<c_ss.str();

        }

    }
    catch(std::exception& e)
    {
        Erl::debug_cout_e()<<"Error: readConcentricTubefromXMLV1::read, parsing failed: "<<e.what()<<"; File: "<<_filename<<std::endl;
        assert(false && "Error in Construction of Tool");
    }
    if(R != nullptr)
    {
        c_JVal.conservativeResize(c_IJv,1);
        R->calcKinematic(c_JVal, (R->getMaxLength() / R->getMaxSamples()) );
    }

    return R;
}
#endif
template<class Real>
types::RobotA<Real>* readConcentricTubefromXMLV2(const std::string& _filename, const size_t& _maxSamples, bool _verbose)
{
    size_t c_ISections(0);
    size_t c_ITube(0);
    size_t c_IJv(1);
    types::RobotA<Real> * R(nullptr);

    std::ifstream c_XMLFile(_filename.c_str(),std::ifstream::in);
    if (!c_XMLFile.good())
    {
        std::cerr<<"ERROR: concentricTubefromXML::XMLParse, could not read file: "<<_filename<<std::endl;
        return nullptr;
    }

    using boost::property_tree::ptree;
    ptree c_PropertyTree;
    try
    {
        read_xml(c_XMLFile, c_PropertyTree, boost::property_tree::xml_parser::trim_whitespace);
    }
    catch (std::exception e)
    {
        Erl::debug_cout_e() << "ERROR: " << FILE_STRING << ", Exception: " << e.what() <<", File: "<< _filename.c_str() <<std::endl;
        return nullptr;
    }
    c_XMLFile.close();

    Erl::Transform<Real>  c_InitTrafo = Erl::Transform<Real>::Identity();
    Eigen::Matrix<Real,Eigen::Dynamic,1> c_JVal
            = Eigen::Matrix<Real,Eigen::Dynamic,1>::Zero(1+3*CTR_MAXSECTION,1);

    std::vector<std::pair<types::SecCC<Real>,size_t > > c_VectorCC;
    std::vector<std::pair<types::SecVC<Real>,size_t > > c_VectorVC;
    std::vector< std::string >                          c_SectionName;
    c_VectorCC   .reserve(CTR_MAXSECTION);
    c_VectorVC   .reserve(CTR_MAXSECTION);
    c_SectionName.reserve(CTR_MAXSECTION);

    try
    {
        for( auto v : c_PropertyTree.get_child("robot") )
        {
            std::string c_v_first = v.first;
            std::transform(c_v_first.begin(), c_v_first.end(), c_v_first.begin(), ::tolower);


            if( c_v_first.compare("section")==0 )
            {
                Real c_length        ;
                Real c_phi           ;


                std::array<Real,2> c_stiffness     ;
                std::array<Real,2> c_curvature     ;
                std::array<Real,2> c_radius        ;
                std::array<Real,2> c_possisionRatio;


                c_length         = v.second.get<Real>("length"      );
                c_phi            = v.second.get<Real>("phi"         ,Real(0));

                unsigned c_IdxTube(0);

                for (auto v2 : v.second)
                {


                    std::string c_v2_first = v2.first;
                    std::transform(c_v2_first.begin(), c_v2_first.end(), c_v2_first.begin(), ::tolower);
                    if( c_v2_first.compare("name")==0 )
                    {
                        c_SectionName.push_back(v2.second.data().c_str());
                    }
                    if( c_v2_first.compare("tube")==0 )
                    {
                        if(c_IdxTube>1){throw std::logic_error("VIOLATED: maximum of 2 tubes per section");}
                        c_stiffness     [c_IdxTube] = v2.second.get<Real>("stiffness"   );
                        c_curvature     [c_IdxTube] = v2.second.get<Real>("curvature"   );
                        c_radius        [c_IdxTube] = v2.second.get<Real>("radius"      );
                        c_possisionRatio[c_IdxTube] = v2.second.get<Real>("poissonRatio");
                        c_IdxTube ++; c_ITube ++;
                    }
                }
                if(c_IdxTube==1)
                {
                    types::SecCC<Real> c_Sc({c_stiffness[0]},{c_curvature[0]},c_length,{c_radius[0]},{c_possisionRatio[0]});
                    c_VectorCC.push_back(std::make_pair(c_Sc,c_ISections++));
                    c_JVal(c_IJv++,0) = c_phi; c_IJv+=1;
                }
                else if(c_IdxTube==2)
                {
                    types::SecVC<Real> c_Sv(c_stiffness,c_curvature,c_length,c_radius,c_possisionRatio);
                    c_VectorVC.push_back(std::make_pair(c_Sv,c_ISections++));
                    c_JVal(c_IJv++,0) = c_phi; c_IJv+=2;
                }
                if(c_ISections>CTR_MAXSECTION)
                {
                    throw std::invalid_argument(std::string("MAX SECTION: This Programm supports only: ")+std::to_string(int(CTR_MAXSECTION))+" sections, please adapt your defines (CTR_MAXSECTION) at compilation.");
                }
            }
            if( c_v_first.compare("entryframe")==0 )
            {
                c_InitTrafo = Erl::Transform<Real>(v.second.get<Real>("E11"),v.second.get<Real>("E12"),v.second.get<Real>("E13"),v.second.get<Real>("E14"),
                                                   v.second.get<Real>("E21"),v.second.get<Real>("E22"),v.second.get<Real>("E23"),v.second.get<Real>("E24"),
                                                   v.second.get<Real>("E31"),v.second.get<Real>("E32"),v.second.get<Real>("E33"),v.second.get<Real>("E34"));

                if(!Erl::Transform<Real>::verifyOrthogonal(c_InitTrafo,0.01))
                {
                    throw std::invalid_argument("InitTrafo: Not a Transformationmatrix");
                }
                c_InitTrafo = Erl::Transform<Real>::nearestOrthogonal(c_InitTrafo);
            }


        }

        R = constructRobot<types::RobotA<Real>>(_maxSamples,c_InitTrafo,c_VectorCC,c_VectorVC);

        auto c_IterCC = c_VectorCC.begin();
        auto c_IterVC = c_VectorVC.begin();

        if(_verbose)
        {
            std::stringstream c_ss;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            c_ss<<"File: "<<_filename<<std::endl;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            for(unsigned i(0);i<c_ISections;i++)
            {
                if(i<c_SectionName.size())
                { c_ss<<c_SectionName.at(i)<<std::endl; }

                if(c_IterCC != c_VectorCC.end())
                { if(c_IterCC->second == i){ c_ss<<c_IterCC->first; c_IterCC++;} }
                if(c_IterVC != c_VectorVC.end())
                { if(c_IterVC->second == i){ c_ss<<c_IterVC->first; c_IterVC++;} }
                c_ss<<std::endl;
            }
            c_ss<<"InitialFrame: "<<std::endl;
            Erl::debug_set_io_format(c_ss,4,7,std::ios_base::fixed,true);
            c_ss<<c_InitTrafo<<std::flush;
            c_ss<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;
            Erl::debug_cout()<<c_ss.str();

        }

    }
    catch(std::exception& e)
    {
        Erl::debug_cout_e()<<"Error: readConcentricTubefromXMLV1::read, parsing failed: "<<e.what()<<"; File: "<<_filename<<std::endl;
        assert(false && "Error in Construction of Tool");
    }
    if(R != nullptr)
    {
        c_JVal.conservativeResize(c_IJv,1);
        R->calcKinematic(c_JVal, (R->getMaxLength() / R->getMaxSamples()) );
    }

    return R;
}

template<class Real>
types::RobotA<Real>* readConcentricTubefromXML(const std::string& _filename, const size_t& _maxSamples, bool _verbose)
{
    std::ifstream c_XMLFile(_filename.c_str(),std::ifstream::in);
    if (!c_XMLFile.good())
    {
        std::cerr<<"ERROR: concentricTubefromXML::XMLParse, could not read file: "<<_filename<<std::endl;
        return nullptr;
    }
    //****
    // CHECK FOR VERSION
    std::string c_String;
    while(c_XMLFile >> c_String)
    {
        if(c_String.find("version=\"2.0\"")!=std::string::npos)
        {
            c_XMLFile.close();
            return readConcentricTubefromXMLV2<Real>(_filename,_maxSamples,_verbose);

        }
    }
    if(c_XMLFile.is_open()){c_XMLFile.close();}
    return readConcentricTubefromXMLV1<Real>(_filename,_maxSamples,_verbose);
}



}
}
#endif
