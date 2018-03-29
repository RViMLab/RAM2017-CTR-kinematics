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

#include <Erl/Utility/debug.h>
#include "ctr_robot.h"

int main(int _argc, char *_argv[])
{
    typedef double Real;
    typedef CTR::Robot<Real> Robot_t;
    if(_argc>2)
    {
        size_t      c_sample       = std::stoul(_argv[1]);
        std::string c_xml_filename = _argv[2];
        Robot_t c_robot = Robot_t(c_sample, c_xml_filename);
        CT_RNG c_rng; CT_RND<Real> c_rnd;
        Robot_t::VectorJ c_joint_values(c_robot.getNJointValues(),1);
        c_robot.getRandomJointValues(c_rng,c_rnd,c_joint_values);
        Erl::debug_cout_g()<<"Joint values: "<<c_joint_values.transpose()<<std::endl;
        return 0;
    }
    Erl::debug_cout_e()<<"Not enough input arguments."<<std::endl;
    return 1;
}
