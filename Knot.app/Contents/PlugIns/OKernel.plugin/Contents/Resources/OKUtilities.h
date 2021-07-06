//
//  OKUtilities.h
//  OKernel
//
//  Copyright (c) 2013-2021 All rights reserved.
//

@interface OKUtilities:NSObject {
    
}

+ (void)postFileEvents;
+ (StimSetting *)getStimSetting;
+ (void)requestReset;
+ (void)reset;
+ (BOOL)setArrayToCount:(StimType)stimType;
+ (void)setArrayValues:(StimType)stimType;
+ (void)updateBlockStatus;

@end

