//
//  OKEndtrialState.h
//  Experiment
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "OKStateSystem.h"

@interface OKEndtrialState : LLState<NSSpeechSynthesizerDelegate>
{
	NSTimeInterval      expireTime;
}

@property (assign) LLReward   reward;

@end
