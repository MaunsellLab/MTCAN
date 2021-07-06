//
//  GRFIntertrialState.h
//  Experiment
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "GRFStateSystem.h"

@interface GRFIntertrialState : LLState {

    NSTimeInterval    expireTime;
}

@property (readonly) BOOL selectTrial;

@end
