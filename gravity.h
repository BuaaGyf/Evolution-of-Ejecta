//
//  gravity.h
//  
//
//  Created by Yang Yu and Yunfeng Gao on 2022-5-1.
//
//

#ifndef _GRAVITY_H_
#define _GRAVITY_H_

#include "body.h"

void GravAttraction(fem &alpha, Vector r, Vector F);
void GravBeta(fem &beta, Vector r, Vector F);

#endif /* defined(____gravity__) */
