/*
FTStimuli.h
*/

#import <Lablib/Lablib.h>

@interface FTStimuli : NSObject <LLVisualStimOwner> {

    float   fixSizePix;
    BOOL    fixSpotOn;
    long    frameCounter;
}

- (void)drawFixSpot;
- (void)erase;

@property BOOL fixOn;
@property (strong) LLFixTarget *fixSpot;
@property (strong) LLIntervalMonitor *monitor;

@end
