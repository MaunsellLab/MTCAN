//
//  FTBlockedState.h
//  Fixate Task
//
//  Created by John Maunsell on Thu May 29 2003.
//  Copyright 2003-2021. All rights reserved.
//

#import "FTStateSystem.h"

@interface FTBlockedState : LLState {

	NSTimeInterval	expireTime;
}

@end
