/*
VFStimuli.h
*/

#import "VF.h"
#import "VFMapStimTable.h"

@interface VFStimuli : NSObject  <LLVisualStimOwner> {
    long                currentStimIndex;
    StimDesc            currentStimDescs;
    DisplayParam        display;
    BOOL                doDigitalCodes;
    NSMutableArray      *mapStimList;
}

@property BOOL abortStimuli;
@property long firstFrame;
@property (strong) LLGabor *mappingGabor;
@property (readonly, strong) LLIntervalMonitor *monitor;
@property long stimFrame;
@property BOOL stimSequenceUnderway;
@property BOOL stimulusOn;

- (void)clearStimList;
- (void)makeStimList;
- (LLGabor *)mappingGabor;
- (void)releaseStimuli;
- (void)startStimSequence;
- (void)stopAllStimuli;
- (void)tallyStimList:(long)count;

@end
