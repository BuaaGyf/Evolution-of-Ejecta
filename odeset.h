//
//  ode.h
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//
#ifndef _ODESET_H_
#define _ODESET_H_

typedef struct ode_option {
    //
    int FunOption, EndStep, StartStep, OutputInterval;
    double StepSize;
    
} ODE_OPTION;

#endif /* defined(____binary__) */
