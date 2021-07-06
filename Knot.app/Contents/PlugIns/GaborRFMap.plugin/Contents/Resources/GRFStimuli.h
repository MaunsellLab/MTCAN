/*
GRFStimuli.h
*/

#import "GRF.h"
#import "GRFMapStimTable.h"

@interface GRFStimuli : NSObject  <LLVisualStimOwner> {
    StimDesc            currentStimDescs[kGabors];
    long                currentStimIndex[kGabors];
    DisplayParam        display;
    long                durationMS;
    float               fixSizePix;
    LLFixTarget         *fixSpot;
    BOOL                fixSpotOn;
    NSArray             *fixTargets;
    NSArray             *gabors;
    NSArray             *gaborStimLists;
    BOOL                hideLeftCodes;
    BOOL                hideRightCodes;
    NSMutableArray      *mapStimList0;
    NSMutableArray      *mapStimList1;
    short               selectTable[kMaxOriChanges];
    NSMutableArray      *taskStimList;
//    TrialDesc           trial;
    LLPlaid             *plaid;
    BOOL                usePlaid;
}

@property BOOL abortStimuli;
@property long firstFrame;
@property long numMap0Completed;
@property long numMap1Completed;
@property (readonly, strong) LLIntervalMonitor *monitor;
@property long stimFrame;
@property BOOL stimSequenceUnderway;
@property BOOL stimulusOn;
@property long targetOnFrame;
@property double trialStartTimeS;

- (void)clearStimLists:(TrialDesc *)pTrial;
- (void)doFixSettings;
- (void)doGabor0Settings;
- (void)dumpStimList:(NSArray *)stimList;
- (void)fixSpotOn;
- (LLGabor *)initializeGabor:(BOOL)bindTemporalFreq;
- (void)loadPlaid:(LLPlaid *)plaid withStimDesc0:(StimDesc *)pSD0 withStimDesc1:(StimDesc *)pSD1;
- (void)makeStimLists;
- (LLGabor *)mappingGabor0;
- (LLGabor *)mappingGabor1;
- (void)releaseStimuli;
- (void)shuffleStimListFrom:(short)start count:(short)count;
- (void)startStimSequence;
- (void)stopAllStimuli;
- (LLGabor *)taskGabor;

@end
