/*
MTCStimuli.h
*/

#import "MTC.h"

@interface MTCStimuli : NSObject  <LLVisualStimOwner> {
    DisplayParam        display;
    long                durationMS;
    float               fixSizePix;
    BOOL                fixSpotOn;
    NSArray             *fixTargets;
    BOOL                hideLeftCodes;
    BOOL                hideRightCodes;
    short               selectTable[kMaxDirChanges];
}

@property BOOL abortStimuli;
@property long firstFrame;
@property (readonly, strong) LLFixTarget *fixSpot;
@property (readonly, strong) LLIntervalMonitor *monitor;
@property long stimFrame;
@property (strong) NSMutableArray   *stimList;
@property BOOL stimSequenceUnderway;
@property BOOL stimulusOn;
@property long targetOnFrame;
@property (strong, readonly) NSMutableArray *taskGabors;
@property double trialStartTimeS;

- (void)fixSpotOn;
- (void)makeStimList;
- (void)releaseStimuli;
- (void)shuffleStimListFrom:(short)start count:(short)count;
- (void)shuffleStimSequence;
- (void)startStimuli;
- (void)stopAllStimuli;
- (void)tallyStimuli;

@end
