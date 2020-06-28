//
//  common.hpp
//  Stream_HPC
//
//  Created by changjiang GOU on 20/04/2020.
//  Copyright © 2020 Changjiang GOU. All rights reserved.
//

#ifndef common_hpp
#define common_hpp

#include <stdio.h>

//configurations
const double SMAX=3700;
const double SMIN=1200;
const unsigned Speed_Options=6;
const array<const double, Speed_Options> DiscSpeed={SMIN,2100,2400,2600,3000,SMAX};
extern vector<double> BandWidth;//bandwidth between processors in the same block or different blocks
const array<const double, 2> Alpha={0.2,0.8};//energy cost of sending a unit of data between cores in the same block or different blocks
const double Static_Power=2.17;
const double ConstantC=1;
extern unsigned int Config_Core;
extern double period;

//const double SMAX=4;
//const double SMIN=1;
//const unsigned Speed_Options=3;
//const array<const double, Speed_Options> DiscSpeed={SMIN,2,SMAX};
//extern vector<double> BandWidth;//bandwidth between processors in the same block or different blocks
//const array<const double, 2> Alpha={1,2};//energy cost of sending a unit of data between cores in the same block or different blocks
//const double Static_Power=0.1;
//const double ConstantC=1;
//extern unsigned int Config_Core;
//extern double period;

#endif /* common_hpp */
