/*
OKStimuli.h
*/

#import "OK.h"

@interface OKStimuli : NSObject <LLVisualStimOwner> {

	DisplayParam			display;
    TrialDesc               trial;
}

@property BOOL abortStimuli;
@property (strong) LLGabor *gabor;
@property (strong) LLIntervalMonitor *monitor;
@property float originalContrast;
@property BOOL photodiode;
@property (strong) LLGabor *rampGabor;
@property long stimFrame;
@property BOOL stimSequenceUnderway;
@property BOOL stimulusActive;
@property BOOL stimulusOn;

- (void)backgroundDim;
- (void)backgroundNormal;
- (void)startStimulus:(TrialDesc *)pTrial;
- (void)stopAllStimuli;
- (void)unbindGabor;

@end
