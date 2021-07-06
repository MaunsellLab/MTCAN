//
//  MTCIntertrialState.h
//  Experiment
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "MTCStateSystem.h"

@interface MTCIntertrialState : LLState {

    NSTimeInterval    expireTime;
}

@property (readonly) BOOL selectTrial;

@end
