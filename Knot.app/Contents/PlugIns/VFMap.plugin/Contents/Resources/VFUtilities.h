//
//  VFUtilities.h
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "VF.h"

void			announceEvents(void);
BehaviorSetting *getBehaviorSetting(void);
void			requestReset(void);
void			reset(void);
BOOL			selectTrial(long *pIndex);
float			spikeRateFromStimValue(float normalizedValue);
