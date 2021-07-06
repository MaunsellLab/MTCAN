//
//  MTCAN.h
//  MTCAN
//
//  Copyright 2006-2021 All rights reserved.
//

#import "MTC.h"
#import "MTCStateSystem.h"
#import "MTCDiagnosticsController.h"
#import "MTCEyeXYController.h"
#import "MTCRoundToStimCycle.h"
#import "MTCStateSystem.h"

@class MTCMapStimTable;

@interface MTCAN:LLTaskPlugIn {

    NSWindowController      *behaviorController;
    NSPoint                 currentEyesUnits[kEyes];
    MTCEyeXYController      *eyeXYController;                // Eye position display
    NSWindowController      *spikeController;
    NSWindowController      *summaryController;
    NSArray                 *topLevelObjects;
}

@property (readonly) BlockStatus blockStatus;
@property (readonly, strong) MTCDiagnosticsController *diagnosticsController;
@property (readonly, strong) LLDigitalOut *digitalOut;
@property (strong, readonly) LLEyeWindow *respWindow;
@property (readonly) BlockStatus *pBS;
@property (strong) MTCStimuli *stimuli;
@property TrialDesc trial;

- (IBAction)doFixSettings:(id)sender;
- (IBAction)doJuice:(id)sender;
- (void)doJuiceOff;
- (IBAction)doReset:(id)sender;
- (IBAction)doTaskGaborSettings:(id)sender;

@end

