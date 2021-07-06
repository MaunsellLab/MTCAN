//
//  OKReactState.h
//  Experiment
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "OKStateSystem.h"

@interface OKWaitResponseState : LLState {

	NSTimeInterval expireTime;
}

@end
