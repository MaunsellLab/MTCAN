//
//  OKIntertrialState.h
//  Experiment
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "OKStateSystem.h"

@interface OKIntertrialState : LLState {

    NSTimeInterval    expireTime;
}

@property (readonly) BOOL selectTrial;

@end
