//
//  FTStateSystem.h
//  Fixate Task
//
//  Created by John Maunsell on Thu May 29 2003.
//  Copyright 2003-2021. All rights reserved.
//

#import "Fixate.h"

#define		kFixOnSound				@"6C"
#define		kFixateSound			@"7G"
#define		kStimOnSound			@"5C"
#define		kStimOffSound			@"5C"
#define 	kCorrectSound			@"Correct"
#define 	kNotCorrectSound		@"NotCorrect"

#define kFTFixateJitterPCKey        @"FTFixateJitterPC"
#define kFTFixateMSKey              @"FTFixateMS"
#define kFTFixWindowWidthDegKey     @"FTFixWindowWidthDeg"
#define kFTTaskModeKey              @"FTTaskMode"

typedef struct {
	long	stimulusType;
	long	stimulusIndex;
	float   stimulusValue;
	long	stimulusInterval;
	float   respAziDeg;
	float   respEleDeg;
} TrialDesc;

extern long 			eotCode;			// End Of Trial code
extern LLEyeWindow		*fixWindow;

@interface FTStateSystem : LLStateSystem {

}

@end

